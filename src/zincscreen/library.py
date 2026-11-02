"""Read selected compressed SMILES streams without unpacking the library."""
from contextlib import contextmanager
from pathlib import Path
import bz2
import gzip
import hashlib
import json
import lzma
import random
import re

SUFFIXES = (".smi", ".smiles", ".smi.gz", ".smiles.gz", ".smi.bz2", ".smi.xz", ".gz", ".bz2", ".xz")


def inventory(root):
    root = Path(root).resolve(strict=True)
    if not root.is_dir():
        raise ValueError("zinc_root is not a directory")
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.is_symlink() or not path.name.lower().endswith(SUFFIXES):
            continue
        s = path.stat()
        rel = path.relative_to(root).as_posix()
        match = re.search(r"(?:^|[^A-Za-z0-9])(H\d{2})(?:[^0-9]|$)", rel, re.I)
        yield {"path": rel, "compressed_bytes": s.st_size, "mtime_ns": s.st_mtime_ns,
               "stratum": match.group(1).upper() if match else "unclassified"}


@contextmanager
def text_stream(path):
    path = Path(path)
    opener = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}.get(path.suffix, open)
    with opener(path, "rt", encoding="utf-8", errors="strict") as stream:
        yield stream


def parse_line(line, source, number):
    text = line.strip()
    if not text or text.startswith("#"):
        return None
    fields = text.split()
    if fields[0].lower() in ("smiles", "canonical_smiles"):
        return None
    record = {"source_file": source, "source_line": number,
              "raw_line_sha256": hashlib.sha256(line.encode()).hexdigest()}
    if len(fields) < 2 or not re.fullmatch(r"ZINC[A-Za-z0-9]+", fields[1], flags=re.I):
        return {**record, "status": "INVALID_RECORD", "raw_line": text,
                "error": "Expected SMILES followed by ZINC_ID"}
    return {**record, "id": fields[1], "smiles": fields[0], "status": "INPUT"}


def sample(root, output, *, limit=100, per_file=4, scan_lines=10000, seed=42):
    if min(limit, per_file, scan_lines) < 1:
        raise ValueError("Sample limits must be positive")
    rows = list(inventory(root))
    groups = {}
    rng = random.Random(seed)
    for row in rows:
        groups.setdefault(row["stratum"], []).append(row)
    for values in groups.values():
        rng.shuffle(values)
    # Visit strata in rounds; do not exhaust the first stratum before visiting others.
    ordered = []
    for i in range(max((len(v) for v in groups.values()), default=0)):
        ordered.extend(groups[k][i] for k in sorted(groups) if i < len(groups[k]))
    selected, rejected, visited, seen = [], [], [], set()
    for item in ordered:
        reservoir = []
        invalid = 0
        count = 0
        with text_stream(Path(root) / item["path"]) as stream:
            for number, line in enumerate(stream, 1):
                if number > scan_lines:
                    break
                r = parse_line(line, item["path"], number)
                if r is None:
                    continue
                if r["status"] == "INVALID_RECORD":
                    invalid += 1
                    if len(rejected) < 100:
                        rejected.append(r)
                    continue
                count += 1
                if len(reservoir) < per_file:
                    reservoir.append(r)
                else:
                    j = rng.randrange(count)
                    if j < per_file:
                        reservoir[j] = r
        visited.append({**item, "valid_records_scanned": count, "invalid_records_scanned": invalid})
        for r in reservoir:
            identity = (r["id"], r["smiles"])
            if identity not in seen and len(selected) < limit:
                seen.add(identity)
                selected.append({**r, "source_file_stat": item})
        if len(selected) >= limit:
            break
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError("Sample exists; choose another output to preserve provenance")
    with output.open("x") as stream:
        for r in selected:
            stream.write(json.dumps(r, sort_keys=True) + "\n")
    from .common import save, sha
    result = {"sample_sha256": sha(output), "count": len(selected), "seed": seed,
              "scan_lines_per_file": scan_lines, "reservoir_per_file": per_file,
              "files_visited": visited, "invalid_examples": rejected,
              "limitation": "Bounded prefix reservoirs, stratified by path token, NOT a uniform full-library sample or a molecule-count audit"}
    save(output.with_suffix(output.suffix + ".receipt.json"), result)
    return result
