"""Serial pilot runner: bounded jobs, shared-GPU admission, immutable attempts."""
from pathlib import Path
import csv
import io
import math
import os
import signal
import subprocess
import time
import uuid

from .common import digest, disk_check, executable, load, lock, records, receptors, save, sha
from .preparation import check_prepared


def gpu_snapshot(c):
    if not c["engine"]["use_gpu"]:
        return {"used": False}
    if os.environ.get("CUDA_VISIBLE_DEVICES"):
        raise ValueError("Unset CUDA_VISIBLE_DEVICES; this pilot uses the configured physical device index")
    query = ["nvidia-smi", "-i", str(c["engine"]["device"]),
             "--query-gpu=uuid,name,memory.total,memory.free,utilization.gpu,driver_version", "--format=csv,noheader,nounits"]
    proc = subprocess.run(query, text=True, capture_output=True, timeout=10, check=True)
    rows = list(csv.reader(io.StringIO(proc.stdout)))
    if len(rows) != 1 or len(rows[0]) != 6:
        raise RuntimeError("Cannot identify exactly one GPU")
    gpu_id, name, total, free, utilization, driver = [v.strip() for v in rows[0]]
    row = {"used": True, "uuid": gpu_id, "name": name, "total_mib": int(total),
           "free_mib": int(free), "utilization_percent": int(utilization), "driver": driver,
           "observed_unix": time.time()}
    row["admitted"] = row["free_mib"] >= c["gpu"]["min_free_mib"] and row["utilization_percent"] <= c["gpu"]["max_utilization_percent"]
    row["limitation"] = "Point-in-time check, not an exclusive reservation; coordinate shared workstation usage"
    return row


def command(c, exe, receptor, ligand, attempt, seed):
    e = c["engine"]
    output = attempt / ("poses.sdf" if e["backend"] == "gnina" else "poses.pdbqt")
    argv = [str(exe), "--receptor", receptor["pdbqt"], "--ligand", str(ligand / "ligand.pdbqt"),
            "--scoring", "vina", "--cpu", str(e["cpu_threads"]), "--seed", str(seed),
            "--exhaustiveness", str(e["exhaustiveness"]), "--num_modes", str(e["num_modes"]),
            "--energy_range", str(e["energy_range"]), "--out", str(output)]
    for key in ("center", "size"):
        for axis, value in zip("xyz", receptor["box"][key]):
            argv += ["--" + key + "_" + axis, str(value)]
    if e["backend"] == "gnina":
        argv += ["--cnn_scoring", "rescore", "--pose_sort_order", "Energy", "--addH", "false"]
        argv += ["--device", str(e["device"])] if e["use_gpu"] else ["--no_gpu"]
    return argv, output


def analyze(c, attempt, ligand, receptor):
    from .chemistry import pdbqt_atoms
    from .poses import vina_poses, gnina_poses, geometry_flags
    if c["engine"]["backend"] == "gnina":
        poses = gnina_poses(attempt / "poses.sdf", ligand / "ligand.sdf")
    else:
        poses = vina_poses(attempt / "poses.pdbqt", pdbqt_atoms(ligand / "ligand.pdbqt"))
    mapping = load(ligand / "mapping.json")
    for pose in poses:
        pose.update(geometry_flags(pose, ligand / "ligand.sdf", mapping, receptor["pdbqt"], receptor["box"]))
    save(attempt / "POSES.json", poses)
    return poses


def stop_owned_process(proc):
    if proc.poll() is None:
        os.killpg(proc.pid, signal.SIGTERM)
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(timeout=10)


def process_identity(pid):
    """Linux process start ticks plus boot ID avoid accepting a recycled PID."""
    try:
        stat = Path("/proc", str(pid), "stat").read_text().rsplit(")", 1)[1].split()
        return {"pid": pid, "start_ticks": stat[19], "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip()}
    except FileNotFoundError:
        return None


def run(c, prepared, receptor_ids, *, seed=42, max_jobs=20, retry_failed=False):
    if max_jobs < 1 or not 0 <= seed < 2147483645:
        raise ValueError("Invalid seed/job limit")
    workspace = Path(c["workspace"])
    disk_check(c)
    recs = receptors(c)
    selected = [recs[name] for name in receptor_ids]
    exe = executable(c)
    engine_sha = sha(exe)
    implementation = digest({p.name: sha(p) for p in sorted(Path(__file__).parent.glob("*.py"))})
    counts = {"executed": 0, "reanalyzed": 0, "reused": 0, "held_failures": 0, "DONE": 0, "failed": 0}
    with lock(workspace / "runner.lock"):
        for row in records(prepared):
            if row["status"] != "READY":
                continue
            ligand, receipt = check_prepared(row, workspace)
            for receptor in selected:
                if counts["executed"] + counts["reanalyzed"] >= max_jobs:
                    return counts
                spec = {"input": row["input"], "preparation_sha256": row["receipt_sha256"],
                        "ligand_files": receipt["files"], "receptor_id": receptor["id"],
                        "receptor_sha256": receptor["pdbqt_sha256"], "box_sha256": receptor["box_file_sha256"],
                        "box": receptor["box"], "engine": c["engine"], "engine_sha256": engine_sha,
                        "implementation_sha256": implementation, "seed": seed}
                key = digest(spec)
                task = workspace / "tasks" / key[:2] / key
                task.mkdir(parents=True, exist_ok=True)
                result_path = task / "RESULT.json"
                previous = load(result_path) if result_path.exists() else None
                if previous and previous["status"] == "DONE":
                    for rel, expected in previous["files"].items():
                        if sha(task / rel) != expected:
                            raise ValueError("Completed result was modified; no silent redocking: " + key)
                    counts["reused"] += 1
                    continue
                if previous and not retry_failed:
                    counts["held_failures"] += 1
                    continue
                # Successful engine output with failed parsing is re-analyzed, never re-docked.
                retained_attempt = task / previous["attempt"] if previous else None
                done_path = retained_attempt / "ENGINE_DONE.json" if retained_attempt else None
                reanalyze = previous and (previous["status"] == "ANALYSIS_FAILED" or done_path.exists())
                if previous and previous["status"] == "RUNNING" and not reanalyze:
                    identity = previous.get("process_identity")
                    if not identity or process_identity(identity["pid"]) == identity:
                        raise RuntimeError("Retained engine may still be alive or its launch is uncertain; inspect before retrying")
                gpu = {"used": False, "analysis_only": True} if reanalyze else gpu_snapshot(c)
                if gpu.get("used") and not gpu["admitted"]:
                    save(workspace / "GPU_WAIT.json", gpu)
                    return {**counts, "waiting_for_gpu": True}
                disk_check(c)
                save(task / "SPEC.json", spec)
                if reanalyze:
                    attempt = retained_attempt
                    retained = load(done_path) if done_path.exists() else previous
                    output = attempt / retained["output"]
                    if sha(output) != retained["output_sha256"]:
                        raise ValueError("Retained native output changed")
                    counts["reanalyzed"] += 1
                else:
                    attempt = task / ("attempt-" + uuid.uuid4().hex)
                    attempt.mkdir()
                    argv, output = command(c, exe, receptor, ligand, attempt, seed)
                    save(attempt / "PLAN.json", {"argv": argv, "gpu": gpu, "spec": spec})
                    started = time.time()
                    save(result_path, {"status": "RUNNING", "attempt": attempt.name, "started_unix": started})
                    env = dict(os.environ)
                    env.update(OMP_NUM_THREADS=str(c["engine"]["cpu_threads"]), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NVIDIA_TF32_OVERRIDE="0")
                    if not c["engine"]["use_gpu"]:
                        env["CUDA_VISIBLE_DEVICES"] = ""
                    counts["executed"] += 1
                    try:
                        with (attempt / "engine.stdout.log").open("w") as stdout, (attempt / "engine.stderr.log").open("w") as stderr:
                            proc = subprocess.Popen(argv, cwd=attempt, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
                            try:
                                save(result_path, {"status": "RUNNING", "attempt": attempt.name, "started_unix": started,
                                                   "process_identity": process_identity(proc.pid)})
                                code = proc.wait(timeout=c["engine"]["timeout_seconds"])
                            except BaseException:
                                stop_owned_process(proc)
                                raise
                        if code:
                            raise RuntimeError("Engine exit code " + str(code))
                        save(attempt / "ENGINE_DONE.json", {"finished_unix": time.time(), "wall_seconds": time.time() - started,
                                                             "output": output.name, "output_sha256": sha(output)})
                    except (OSError, RuntimeError, subprocess.TimeoutExpired, KeyboardInterrupt) as exc:
                        save(result_path, {"status": "INTERRUPTED" if isinstance(exc, KeyboardInterrupt) else "FAILED_RETAINED",
                                           "attempt": attempt.name, "error": str(exc), "finished_unix": time.time()})
                        counts["failed"] += 1
                        if isinstance(exc, KeyboardInterrupt):
                            raise
                        continue
                try:
                    poses = analyze(c, attempt, ligand, receptor)
                    files = {str(p.relative_to(task)): sha(p) for p in attempt.iterdir() if p.is_file()}
                    result = {"status": "DONE", "id": row["id"], "target": receptor["target"], "receptor": receptor["id"],
                              "backend": c["engine"]["backend"], "seed": seed, "attempt": attempt.name,
                              "minimum_energy_kcal_mol": min(p["energy_kcal_mol"] for p in poses),
                              "qc_eligible_pose_count": sum(p["qc_eligible"] for p in poses),
                              "pose_count": len(poses), "files": files, "finished_unix": time.time(),
                              "scientific_admission": False}
                    save(result_path, result)
                    counts["DONE"] += 1
                    print("DONE", row["id"], receptor["id"], result["minimum_energy_kcal_mol"], flush=True)
                except Exception as exc:
                    save(result_path, {"status": "ANALYSIS_FAILED", "attempt": attempt.name, "error": str(exc),
                                       "output": output.name, "output_sha256": sha(output)})
                    counts["failed"] += 1
    return counts


def report(c, output):
    output = Path(output)
    fields = ["id", "target", "receptor", "backend", "method_id", "seed", "status", "minimum_energy_kcal_mol", "minimum_qc_eligible_energy_kcal_mol", "pose_count", "qc_eligible_pose_count", "task_directory", "error"]
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    target_by_id = {r["id"]: r["target"] for r in receptors(c).values()}
    for path in sorted((Path(c["workspace"]) / "tasks").glob("*/*/RESULT.json")):
        result = load(path)
        spec = load(path.parent / "SPEC.json")
        row = {k: result.get(k, "") for k in fields}
        row.update(id=spec["input"]["id"], target=target_by_id[spec["receptor_id"]],
                   receptor=spec["receptor_id"], backend=spec["engine"]["backend"], seed=spec["seed"], task_directory=str(path.parent))
        row["method_id"] = digest({k: spec[k] for k in ("engine", "engine_sha256", "implementation_sha256", "receptor_sha256", "box_sha256")})
        if result["status"] == "DONE":
            for rel, expected in result["files"].items():
                if sha(path.parent / rel) != expected:
                    raise ValueError("Result changed: " + str(path))
            poses = load(path.parent / result["attempt"] / "POSES.json")
            eligible = [p["energy_kcal_mol"] for p in poses if p["qc_eligible"]]
            row["minimum_qc_eligible_energy_kcal_mol"] = min(eligible) if eligible else ""
        rows.append(row)
    rows.sort(key=lambda r: (r["target"], r["backend"], r["receptor"], float(r["minimum_energy_kcal_mol"]) if r["minimum_energy_kcal_mol"] != "" else math.inf))
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return {"rows": len(rows), "output": str(output), "limitation": "Per-receptor task table. Raw scores across targets/backends are not calibrated selectivity."}
