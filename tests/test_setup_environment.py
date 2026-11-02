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
    assert Path(plan["create"][0]).name == "micromamba"


def test_old_verified_download_is_reused_with_supported_basename(tmp_path, monkeypatch):
    cached = tmp_path / "bin" / ("micromamba-" + bootstrap["VERSION"])
    cached.parent.mkdir()
    cached.write_bytes(b"verified fixture executable")
    copied = []
    def verified(path):
        assert Path(path).read_bytes() == cached.read_bytes()
        copied.append(Path(path))
    def forbid(*args, **kwargs):
        raise AssertionError("Cached recovery must not download or run commands")
    monkeypatch.setitem(bootstrap["ensure_micromamba"].__globals__, "check_download", verified)
    monkeypatch.setattr(bootstrap["subprocess"], "run", forbid)
    binary = bootstrap["ensure_micromamba"](tmp_path)
    assert binary.name == "micromamba" and binary.read_bytes() == cached.read_bytes()
    assert binary.stat().st_mode & 0o111
    assert cached.exists() and len(copied) == 2


def test_execution_overrides_stale_mamba_exe(tmp_path, monkeypatch):
    monkeypatch.setenv("MAMBA_EXE", "/old/micromamba-version")
    binary = tmp_path / "bin/micromamba"
    env = bootstrap["execution_environment"](binary)
    assert env["MAMBA_EXE"] == str(binary)
    assert env["PATH"].split(":")[0] == str(binary.parent)


def test_existing_environment_skips_network_and_creation(tmp_path, monkeypatch, capsys):
    tools, prefix = tmp_path / "tools", tmp_path / "environment"
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "conda-meta/history").write_text("existing installation")
    calls = []
    def forbid(*args):
        raise AssertionError("Existing environment recovery must not require DNS")
    def run(argv, **kwargs):
        assert kwargs["env"]["MAMBA_EXE"] == str(tools / "bin/micromamba")
        calls.append(argv)
    globals_ = bootstrap["main"].__globals__
    monkeypatch.setitem(globals_, "ensure_micromamba", lambda tools: tools / "bin/micromamba")
    monkeypatch.setattr(bootstrap["socket"], "getaddrinfo", forbid)
    monkeypatch.setattr(bootstrap["subprocess"], "run", run)
    bootstrap["main"](["--channel", "tuna", "--tools-root", str(tools), "--prefix", str(prefix)])
    assert all("create" not in argv for argv in calls)
    assert len(calls) == 4  # version, prefix check, project install, dependency check
    assert "Wrong Python environment" in calls[1][-2]
    assert (prefix / "conda-meta/history").read_text() == "existing installation"
    assert (tools / "SETUP_RECEIPT.json").exists()
