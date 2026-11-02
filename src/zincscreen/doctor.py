"""Read-only environment checks, requiring no chemistry imports to start."""
from importlib import import_module, metadata
from pathlib import Path
import platform
import shutil
import subprocess
import sys

from .common import executable, receptors, sha


def doctor(c):
    result = {"python": sys.version, "platform": platform.platform(), "checks": {},
              "jobs_submitted": 0, "scope": "Environment/input integrity only; not scientific validation"}
    checks = result["checks"]
    for name in ("rdkit", "openbabel.openbabel", "numpy", "scipy"):
        try:
            module = import_module(name)
            if name == "openbabel.openbabel":
                version = module.OBReleaseVersion()
            else:
                version = getattr(module, "__version__", "unknown")
            checks[name] = {"ok": True, "version": version}
        except Exception as exc:
            checks[name] = {"ok": False, "error": str(exc)}
    try:
        rec = receptors(c)
        checks["receptors"] = {"ok": True, "ids": sorted(rec)}
    except Exception as exc:
        checks["receptors"] = {"ok": False, "error": str(exc)}
    checks["zinc_root"] = {"ok": Path(c["zinc_root"]).is_dir(), "path": c["zinc_root"]}
    try:
        exe = executable(c)
        proc = subprocess.run([str(exe), "--help"], text=True, capture_output=True, timeout=30)
        help_text = proc.stdout + proc.stderr
        required = ["--receptor", "--ligand", "--seed", "--exhaustiveness", "--num_modes"]
        if c["engine"]["backend"] == "gnina":
            required += ["--cnn_scoring", "--pose_sort_order", "--addH"]
        checks["engine"] = {"ok": proc.returncode == 0 and all(flag in help_text for flag in required),
                            "path": str(exe), "sha256": sha(exe), "help_excerpt": help_text[:3000]}
    except Exception as exc:
        checks["engine"] = {"ok": False, "error": str(exc)}
    if c["engine"]["use_gpu"]:
        try:
            from .runner import gpu_snapshot
            snap = gpu_snapshot(c)
            checks["gpu"] = {"ok": snap["admitted"], **snap}
        except Exception as exc:
            checks["gpu"] = {"ok": False, "error": str(exc)}
    result["ready_for_pilot"] = all(item["ok"] for item in checks.values())
    return result
