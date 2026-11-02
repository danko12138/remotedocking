"""Small durable records and path-independent configuration."""
from contextlib import contextmanager
from pathlib import Path
import fcntl
import hashlib
import json
import math
import os
import shutil
import time


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp-" + str(os.getpid()))
    with tmp.open("w") as out:
        json.dump(value, out, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    tmp.replace(path)


@contextmanager
def lock(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("Another local process owns " + str(path)) from None
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def config(path):
    path = Path(path).resolve(strict=True)
    c = load(path)
    for key in ("workspace", "zinc_root", "receptors_manifest", "references_manifest"):
        p = Path(c[key]).expanduser()
        c[key] = str((path.parent / p).resolve())
    e = c["engine"]
    if e["backend"] not in ("vina", "gnina"):
        raise ValueError("Supported adapters: vina, gnina")
    if e["backend"] == "vina" and e["use_gpu"]:
        raise ValueError("This Vina adapter is CPU-only; GPU execution uses GNINA")
    for k in ("cpu_threads", "exhaustiveness", "num_modes", "timeout_seconds"):
        if type(e[k]) is not int or e[k] < 1:
            raise ValueError("Positive integer required: engine." + k)
    if not math.isfinite(e["energy_range"]) or e["energy_range"] <= 0:
        raise ValueError("Invalid energy_range")
    return c


def executable(c):
    name = c["engine"]["executable"]
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError("Configured engine not found: " + name)
    return Path(path).resolve()


def receptors(c):
    manifest = Path(c["receptors_manifest"])
    result = {}
    for r in load(manifest)["receptors"]:
        item = dict(r)
        for key in ("pdbqt", "box_file"):
            p = (manifest.parent / r[key]).resolve(strict=True)
            if sha(p) != r[key + "_sha256"]:
                raise ValueError("Receptor input hash changed: " + r["id"] + "/" + key)
            item[key] = str(p)
        box = load(item["box_file"])
        for k in ("center", "size"):
            if len(box[k]) != 3 or not all(math.isfinite(v) for v in box[k]):
                raise ValueError("Invalid box: " + r["id"])
        if min(box["size"]) <= 0:
            raise ValueError("Nonpositive box dimensions")
        item["box"] = {k: box[k] for k in ("center", "size")}
        result[item["id"]] = item
    return result


def records(path):
    with Path(path).open() as stream:
        for number, line in enumerate(stream, 1):
            if line.strip():
                value = json.loads(line)
                if not isinstance(value, dict):
                    raise ValueError("Expected JSON object at line " + str(number))
                yield value


def disk_check(c):
    root = Path(c["workspace"])
    root.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(root).free
    if free < c["min_free_disk_gib"] * 1024 ** 3:
        raise RuntimeError("Workspace disk headroom exhausted; existing results retained")
    return {"free_bytes": free, "observed_unix": time.time()}
