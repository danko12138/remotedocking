"""One isolated, resource-bounded worker per compound; receipts are immutable."""
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import time

from .common import disk_check, digest, load, lock, records, save, sha


def check_prepared(row, workspace):
    folder = (Path(workspace) / row["folder"]).resolve(strict=True)
    if not folder.is_relative_to(Path(workspace).resolve()):
        raise ValueError("Prepared folder escaped workspace")
    receipt = load(folder / "PREPARATION.json")
    if sha(folder / "PREPARATION.json") != row["receipt_sha256"]:
        raise ValueError("Preparation receipt changed")
    if receipt["status"] != "READY":
        raise ValueError("Input is not READY")
    for name, expected in receipt["files"].items():
        p = (folder / name).resolve(strict=True)
        if not p.is_relative_to(folder) or sha(p) != expected:
            raise ValueError("Prepared file changed: " + name)
    return folder, receipt


def prepare(c, input_path=None, *, references=False):
    from rdkit import rdBase
    from openbabel import openbabel as ob
    from . import chemistry
    workspace = Path(c["workspace"])
    disk_check(c)
    ref_manifest = Path(c["references_manifest"])
    if references:
        entries = load(ref_manifest)["compounds"]
    else:
        entries = records(input_path)
    ready = failed = 0
    protocol = {"settings": c["preparation"], "chemistry_sha256": sha(chemistry.__file__),
                "preparation_sha256": sha(__file__), "rdkit_version": rdBase.rdkitVersion,
                "openbabel_version": ob.OBReleaseVersion()}
    manifest = workspace / "prepared" / ("references.jsonl" if references else sha(input_path) + ".jsonl")
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with lock(workspace / "preparation.lock"):
        temp = manifest.with_suffix(".tmp")
        with temp.open("w") as index:
            for record in entries:
                disk_check(c)
                identity = digest({"record": record, "protocol": protocol, "references": references, "version": 1})
                folder = workspace / "prepared" / identity
                folder.mkdir(exist_ok=True)
                result_path = folder / "PREPARATION.json"
                if not result_path.exists():
                    try:
                        if references:
                            from .chemistry import native_mapping
                            for name, info in record["files"].items():
                                source = (ref_manifest.parent / info["path"]).resolve(strict=True)
                                if sha(source) != info["sha256"]:
                                    raise ValueError("Reference input changed: " + name)
                                shutil.copyfile(source, folder / name)
                            mapping = native_mapping(folder / "ligand.sdf", folder / "ligand.pdbqt")
                            save(folder / "mapping.json", mapping)
                            result = {"status": "READY", "inchikey": record["inchikey"],
                                      "method": "REUSED_VERIFIED_LOCAL_PREPARATION",
                                      "source_preparation": record["source_preparation"],
                                      "files": {name: sha(folder / name) for name in ("ligand.sdf", "ligand.pdbqt", "mapping.json")}}
                            save(result_path, result)
                        else:
                            save(folder / "INPUT.json", record)
                            env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
                            with (folder / "worker.log").open("w") as log:
                                proc = subprocess.run([sys.executable, "-B", "-m", "zincscreen", "_prepare",
                                                       "--folder", str(folder), "--settings", json.dumps(c["preparation"])],
                                                      stdout=log, stderr=subprocess.STDOUT, env=env,
                                                      timeout=c["preparation"]["timeout_seconds"] + 15)
                            if not result_path.exists():
                                raise RuntimeError("Preparation worker exited without receipt: " + str(proc.returncode))
                    except (ValueError, OSError, RuntimeError, subprocess.TimeoutExpired) as exc:
                        save(result_path, {"status": "PREP_FAILED", "error": str(exc)})
                result = load(result_path)
                row = {"id": record.get("id"), "input": record, "folder": str(folder.relative_to(workspace)),
                       "status": result["status"], "receipt_sha256": sha(result_path)}
                if result["status"] == "READY":
                    check_prepared(row, workspace)
                    ready += 1
                else:
                    row["error"] = result.get("error")
                    failed += 1
                index.write(json.dumps(row, sort_keys=True) + "\n")
                index.flush()
                print(json.dumps({"id": row["id"], "status": row["status"]}), flush=True)
        temp.replace(manifest)
    return {"manifest": str(manifest), "READY": ready, "PREP_FAILED": failed}


def worker(folder, settings):
    import resource
    folder = Path(folder)
    cpu = int(settings["timeout_seconds"])
    memory = int(settings["memory_limit_gib"] * 1024 ** 3)
    resource.setrlimit(resource.RLIMIT_CPU, (cpu, cpu + 5))
    resource.setrlimit(resource.RLIMIT_AS, (memory, memory))
    try:
        from .chemistry import prepare_one
        result = prepare_one(load(folder / "INPUT.json"), folder, settings["seed"])
    except Exception as exc:
        result = {"status": "PREP_FAILED", "error": str(exc)}
    result["finished_unix"] = time.time()
    save(folder / "PREPARATION.json", result)
