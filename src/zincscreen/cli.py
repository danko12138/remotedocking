"""Explicit CLI steps; cloning/importing never starts a job."""
import argparse
import json
import sys
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description="ABCB1/ABCG2 ZINC22 portable pilot")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("doctor", "catalog", "sample", "prepare", "run", "report"):
        p = sub.add_parser(name)
        p.add_argument("--config", required=True)
        if name in ("doctor", "catalog", "sample", "report"):
            p.add_argument("--output", required=True)
        if name == "sample":
            p.add_argument("--limit", type=int, default=100)
            p.add_argument("--per-file", type=int, default=4)
            p.add_argument("--scan-lines", type=int, default=10000)
            p.add_argument("--seed", type=int, default=42)
        if name == "prepare":
            group = p.add_mutually_exclusive_group(required=True)
            group.add_argument("--input")
            group.add_argument("--references", action="store_true")
        if name == "run":
            p.add_argument("--prepared", required=True)
            p.add_argument("--receptors", nargs="+", required=True)
            p.add_argument("--seed", type=int, default=42)
            p.add_argument("--max-jobs", type=int, default=20)
            p.add_argument("--retry-failed", action="store_true")
    worker = sub.add_parser("_prepare", help=argparse.SUPPRESS)
    worker.add_argument("--folder", required=True)
    worker.add_argument("--settings", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "_prepare":
            from .preparation import worker
            worker(args.folder, json.loads(args.settings))
            return
        from .common import config, save
        c = config(args.config)
        if args.command == "doctor":
            from .doctor import doctor
            result = doctor(c)
            save(args.output, result)
        elif args.command == "catalog":
            from .library import inventory
            path = Path(args.output)
            path.parent.mkdir(parents=True, exist_ok=True)
            count = total = 0
            with path.open("x") as stream:
                for record in inventory(c["zinc_root"]):
                    stream.write(json.dumps(record, sort_keys=True) + "\n")
                    count += 1
                    total += record["compressed_bytes"]
            result = {"files": count, "compressed_bytes": total, "output": str(path), "molecules_counted": False}
        elif args.command == "sample":
            from .library import sample
            result = sample(c["zinc_root"], args.output, limit=args.limit, per_file=args.per_file, scan_lines=args.scan_lines, seed=args.seed)
        elif args.command == "prepare":
            from .preparation import prepare
            result = prepare(c, args.input, references=args.references)
        elif args.command == "run":
            from .runner import run
            result = run(c, args.prepared, args.receptors, seed=args.seed, max_jobs=args.max_jobs, retry_failed=args.retry_failed)
        elif args.command == "report":
            from .runner import report
            result = report(c, args.output)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        if args.command == "doctor" and not result["ready_for_pilot"]:
            raise SystemExit(2)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(2, "ERROR: " + str(exc) + "\n")
