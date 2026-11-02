#!/usr/bin/env python3
"""Personal-prefix bootstrap using GitHub releases; requires only stdlib + curl."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import shutil
import socket
import subprocess
import sys
from urllib.parse import urlparse

VERSION = "2.10.0-0"
URL = "https://github.com/mamba-org/micromamba-releases/releases/download/" + VERSION + "/micromamba-linux-64"
# Official GitHub release asset digest, checked when adding this bootstrap.
EXPECTED_SHA256 = "701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd"
CHANNELS = {
    "official": "https://conda.anaconda.org/conda-forge",
    "tuna": "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def mirror_environment(text, channel):
    needle = "  - conda-forge\n"
    if text.count(needle) != 1:
        raise ValueError("Expected exactly one conda-forge channel; inspect environment.yml")
    return text.replace(needle, "  - " + channel + "\n", 1)


def check_download(path, expected=EXPECTED_SHA256):
    if sha(path) != expected:
        raise ValueError("Micromamba SHA256 mismatch; downloaded file will not be executed")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channel", choices=CHANNELS, default="official")
    parser.add_argument("--tools-root", type=Path)
    parser.add_argument("--prefix", type=Path)
    parser.add_argument("--dry-run", action="store_true", help="Print plan without downloads or filesystem writes")
    args = parser.parse_args(argv)
    if sys.platform != "linux" or platform.machine() != "x86_64":
        parser.error("This bootstrap currently supports Linux x86_64 / WSL2 only")
    repo = Path(__file__).resolve().parents[1]
    tools = (args.tools_root or repo.parent / ".abc-zinc-tools").expanduser().resolve()
    prefix = (args.prefix or repo.parent / "envs/abc-zinc22").expanduser().resolve()
    binary = tools / "bin" / ("micromamba-" + VERSION)
    root = tools / "root"
    channel = CHANNELS[args.channel]
    text = mirror_environment((repo / "environment.yml").read_text(), channel)
    envfile = tools / ("environment-" + args.channel + ".yml")
    create = [str(binary), "--no-rc", "-r", str(root), "create", "--override-channels",
              "--strict-channel-priority", "-c", channel, "-y", "-p", str(prefix), "-f", str(envfile)]
    base = [str(binary), "--no-rc", "-r", str(root), "run", "-p", str(prefix)]
    install = base + ["python", "-m", "pip", "install", "--no-index", "--no-build-isolation", "--no-deps", "-e", str(repo)]
    validate = base + ["python", "-c", "import sys,rdkit,numpy,scipy; from openbabel import openbabel; print('Python:',sys.executable); print('DEPENDENCIES_OK')"]
    plan = {"micromamba_url": URL, "micromamba_sha256": EXPECTED_SHA256,
            "tools_root": str(tools), "prefix": str(prefix), "channel": channel,
            "create": create, "install": install, "validate": validate,
            "changes_system_configuration": False, "docking_jobs_submitted": 0}
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return
    try:
        if shutil.which("curl") is None:
            raise RuntimeError("curl is required for downloading the verified executable")
        for host in ("github.com", urlparse(channel).hostname):
            try:
                socket.getaddrinfo(host, 443)
            except socket.gaierror as exc:
                raise RuntimeError("DNS resolution failed for " + host + "; no DNS settings were changed") from exc
        tools.mkdir(parents=True, exist_ok=True)
        binary.parent.mkdir(parents=True, exist_ok=True)
        if binary.exists():
            check_download(binary)
        else:
            temp = binary.with_suffix(".download")
            subprocess.run(["curl", "-fL", "--connect-timeout", "15", "--max-time", "300", "--retry", "2",
                            URL, "-o", str(temp)], check=True)
            check_download(temp)
            temp.chmod(0o755)
            temp.replace(binary)
        subprocess.run([str(binary), "--version"], check=True)
        envfile.write_text(text)
        if (prefix / "conda-meta/history").exists():
            print("Existing environment retained; proceeding to import checks/installing this project only.", flush=True)
        elif prefix.exists():
            raise RuntimeError("Incomplete/unrecognized prefix exists: " + str(prefix) + "; inspect it or use --prefix with a new path")
        else:
            subprocess.run(create, check=True)
        subprocess.run(install, check=True)
        subprocess.run(validate, check=True)
        (tools / "SETUP_RECEIPT.json").write_text(json.dumps({**plan, "status": "DEPENDENCIES_OK", "environment_file_sha256": sha(envfile)}, indent=2) + "\n")
        print("Setup complete. Run the project through:", flush=True)
        print(" ".join(base + ["zincscreen", "--help"]))
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        parser.exit(2, "SETUP FAILED: " + str(exc) + "\nExisting jobs and system DNS were not changed.\n")


if __name__ == "__main__":
    main()
