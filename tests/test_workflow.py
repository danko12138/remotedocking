from pathlib import Path
import gzip
import json
import shutil

import pytest

from zincscreen.common import load, save, sha
from zincscreen.library import parse_line, sample
from zincscreen.chemistry import prepare_one, pdbqt_atoms
from zincscreen.preparation import check_prepared
from zincscreen.poses import vina_poses
from zincscreen.runner import run, report


def test_sample_keeps_provenance_and_visits_strata(tmp_path):
    root = tmp_path / "library"
    for name in ("H04", "H29"):
        (root / name).mkdir(parents=True)
        with gzip.open(root / name / "part.smi.gz", "wt") as out:
            out.write("SMILES ZINC_ID\n")
            out.write("BROKEN_RECORD\n")
            for i in range(20):
                out.write("CCO ZINC" + name[1:] + str(i).zfill(3) + "\n")
    a, b = tmp_path / "a.jsonl", tmp_path / "b.jsonl"
    receipt = sample(root, a, limit=4, per_file=2, scan_lines=30)
    sample(root, b, limit=4, per_file=2, scan_lines=30)
    assert a.read_bytes() == b.read_bytes()
    assert receipt["count"] == 4
    assert {r["stratum"] for r in receipt["files_visited"]} == {"H04", "H29"}
    assert len(receipt["invalid_examples"]) == 2
    assert all(json.loads(line)["source_line"] >= 3 for line in a.read_text().splitlines())
    with pytest.raises(FileExistsError):
        sample(root, a)


def test_input_parser_does_not_invent_zinc_identity():
    assert parse_line("SMILES ZINC_ID\n", "a.gz", 1) is None
    assert parse_line("CCO ../../escaped\n", "a.gz", 2)["status"] == "INVALID_RECORD"
    assert parse_line("CCO ZINC000001\n", "a.gz", 3)["id"] == "ZINC000001"


@pytest.fixture
def prepared(tmp_path):
    workspace = tmp_path / "work"
    folder = workspace / "prepared" / "ethanol"
    receipt = prepare_one({"id": "ZINC000001", "smiles": "CCO"}, folder)
    save(folder / "PREPARATION.json", receipt)
    row = {"id": "ZINC000001", "status": "READY", "folder": "prepared/ethanol",
           "input": {"id": "ZINC000001", "smiles": "CCO"},
           "receipt_sha256": sha(folder / "PREPARATION.json")}
    index = workspace / "prepared.jsonl"
    index.write_text(json.dumps(row) + "\n")
    return workspace, folder, row, index


def test_preparation_holds_salts_and_detects_tampering(prepared, tmp_path):
    workspace, folder, row, _ = prepared
    check_prepared(row, workspace)
    with pytest.raises(ValueError, match="MULTICOMPONENT"):
        prepare_one({"id": "ZINC000002", "smiles": "CCO.Cl"}, tmp_path / "salt")
    with (folder / "ligand.pdbqt").open("a") as out:
        out.write("REMARK tampered\n")
    with pytest.raises(ValueError, match="changed"):
        check_prepared(row, workspace)


def test_pose_parser_rejects_truncation_and_atom_loss(prepared, tmp_path):
    _, folder, _, _ = prepared
    atoms = pdbqt_atoms(folder / "ligand.pdbqt")
    content = (folder / "ligand.pdbqt").read_text()
    path = tmp_path / "poses.pdbqt"
    path.write_text("MODEL 1\nREMARK VINA RESULT: -5.5 0 0\n" + content + "ENDMDL\n")
    assert vina_poses(path, atoms)[0]["energy_kcal_mol"] == -5.5
    path.write_text(path.read_text().replace("ENDMDL", ""))
    with pytest.raises(ValueError, match="truncated"):
        vina_poses(path, atoms)
    path.write_text("MODEL 1\nREMARK VINA RESULT: -5.5 0 0\nENDMDL\n")
    with pytest.raises(ValueError):
        vina_poses(path, atoms)


def fake_config(tmp_path, workspace):
    repo = Path(__file__).resolve().parents[1]
    engine = tmp_path / "fake-vina"
    engine.write_text('''#!/usr/bin/env python3
import pathlib,sys
args=sys.argv[1:]
def value(flag): return args[args.index(flag)+1]
text=pathlib.Path(value('--ligand')).read_text()
out=pathlib.Path(value('--out'))
out.write_text('MODEL 1\\nREMARK VINA RESULT: -5.5 0 0\\n'+text+'ENDMDL\\n')
''')
    engine.chmod(0o755)
    c = load(repo / "configs/churchill.example.json")
    c.update(workspace=str(workspace), receptors_manifest=str(repo / "data/receptors/manifest.json"), min_free_disk_gib=0)
    c["engine"].update(backend="vina", executable=str(engine), use_gpu=False, timeout_seconds=10)
    return c, engine


def test_resume_reuses_results_and_refuses_changed_output(prepared, tmp_path):
    workspace, _, _, index = prepared
    c, _ = fake_config(tmp_path, workspace)
    assert run(c, index, ["9CTF"], max_jobs=1)["DONE"] == 1
    again = run(c, index, ["9CTF"], max_jobs=1)
    assert again["executed"] == 0 and again["reused"] == 1
    output = tmp_path / "scores.csv"
    assert report(c, output)["rows"] == 1
    result = next(workspace.glob("tasks/*/*/RESULT.json"))
    value = load(result)
    native = result.parent / value["attempt"] / "poses.pdbqt"
    native.write_text(native.read_text() + "REMARK changed\n")
    with pytest.raises(ValueError, match="modified"):
        run(c, index, ["9CTF"], max_jobs=1)


def test_analysis_retry_does_not_redock(prepared, tmp_path, monkeypatch):
    workspace, _, _, index = prepared
    c, engine = fake_config(tmp_path, workspace)
    from zincscreen import runner
    real = runner.analyze
    def fail(*args):
        raise ValueError("Simulated parser outage")
    monkeypatch.setattr(runner, "analyze", fail)
    first = run(c, index, ["6VXJ"], max_jobs=1)
    assert first["executed"] == 1 and first["failed"] == 1
    assert run(c, index, ["6VXJ"], max_jobs=1)["held_failures"] == 1
    monkeypatch.setattr(runner, "analyze", real)
    def forbid(*args, **kwargs):
        raise AssertionError("A successful engine must not run twice")
    monkeypatch.setattr(runner.subprocess, "Popen", forbid)
    retried = run(c, index, ["6VXJ"], max_jobs=1, retry_failed=True)
    assert retried["reanalyzed"] == 1 and retried["executed"] == 0 and retried["DONE"] == 1


def test_busy_gpu_waits_instead_of_launching(prepared, tmp_path, monkeypatch):
    workspace, _, _, index = prepared
    c, _ = fake_config(tmp_path, workspace)
    from zincscreen import runner
    monkeypatch.setattr(runner, "gpu_snapshot", lambda c: {"used": True, "admitted": False, "utilization_percent": 89})
    result = run(c, index, ["9CTF"], max_jobs=1)
    assert result["waiting_for_gpu"] and result["executed"] == 0
    assert not list(workspace.glob("tasks/*/*/attempt-*"))


def test_positive_energy_is_retained(prepared, tmp_path):
    _, folder, _, _ = prepared
    path = tmp_path / "positive.pdbqt"
    path.write_text("MODEL 1\nREMARK VINA RESULT: 234.5 0 0\n" + (folder / "ligand.pdbqt").read_text() + "ENDMDL\n")
    assert vina_poses(path, pdbqt_atoms(folder / "ligand.pdbqt"))[0]["energy_kcal_mol"] == 234.5


def test_wrong_full_key_is_not_repaired_by_changing_identity(tmp_path):
    with pytest.raises(ValueError, match="IDENTITY_MISMATCH"):
        prepare_one({"id": "ZINC000001", "smiles": "CCO", "inchikey": "WRONG"}, tmp_path / "wrong")


def test_crash_after_engine_done_resumes_analysis_only(prepared, tmp_path, monkeypatch):
    workspace, _, _, index = prepared
    c, _ = fake_config(tmp_path, workspace)
    from zincscreen import runner
    assert run(c, index, ["9CTF"], max_jobs=1)["DONE"] == 1
    path = next(workspace.glob("tasks/*/*/RESULT.json"))
    old = load(path)
    save(path, {"status": "RUNNING", "attempt": old["attempt"]})
    def forbid(*args, **kwargs):
        raise AssertionError("Engine completed before the crash")
    monkeypatch.setattr(runner.subprocess, "Popen", forbid)
    result = run(c, index, ["9CTF"], max_jobs=1, retry_failed=True)
    assert result["reanalyzed"] == 1 and result["executed"] == 0 and result["DONE"] == 1


@pytest.mark.parametrize("smiles,charge", [("CC[NH3+]", 1), ("CCC(=O)[O-]", -1)])
def test_charged_input_preserves_formal_charge(tmp_path, smiles, charge):
    folder = tmp_path / "charged"
    receipt = prepare_one({"id": "ZINC000001", "smiles": smiles}, folder)
    assert receipt["formal_charge"] == charge
    total = sum(float(line[70:76]) for line in (folder / "ligand.pdbqt").read_text().splitlines() if line.startswith(("ATOM  ", "HETATM")))
    assert abs(total - charge) < .1
