from pathlib import Path
import hashlib
import json
import runpy

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "tools/setup_environment.py"
bootstrap = runpy.run_path(str(SCRIPT))


def test_mirror_preserves_package_pins():
    source = (SCRIPT.parents[1] / "environment.yml").read_text()
    rendered = bootstrap["mirror_environment"](source, bootstrap["CHANNELS"]["tuna"])
    assert "  - conda-forge\n" not in rendered
    assert "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge" in rendered
    assert source.split("dependencies:\n")[1] == rendered.split("dependencies:\n")[1]


def test_bad_download_never_executes(tmp_path):
    path = tmp_path / "micromamba"
    path.write_text("partial or unexpected download")
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        bootstrap["check_download"](path)
    bootstrap["check_download"](path, hashlib.sha256(path.read_bytes()).hexdigest())


def test_dry_run_has_no_writes_and_uses_only_mirror(tmp_path, capsys):
    tools, env = tmp_path / "tools", tmp_path / "env"
    bootstrap["main"](["--dry-run", "--channel", "tuna", "--tools-root", str(tools), "--prefix", str(env)])
    plan = json.loads(capsys.readouterr().out)
    assert not tools.exists() and not env.exists()
    assert "conda.anaconda.org" not in json.dumps(plan)
    assert "github.com/mamba-org/micromamba-releases" in plan["micromamba_url"]
    assert "--no-index" in plan["install"]
    assert plan["docking_jobs_submitted"] == 0
