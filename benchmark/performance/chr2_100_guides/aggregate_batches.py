#!/usr/bin/env python3
import argparse
import csv
import hashlib
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--run-root", required=True)
    p.add_argument("--phase", default="correctness")
    p.add_argument("--iterations", nargs="+", default=["pilot"])
    return p.parse_args()


def main():
    args = parse_args()
    run = Path(args.run_root).resolve()
    guides = read_tsv(SCRIPT_DIR / "guides.tsv")
    batches = read_tsv(SCRIPT_DIR / "batches.tsv")
    expected_guides = [g["guide_id"] for g in guides]
    expected_batches = [b["batch_id"] for b in batches]
    rows = []
    batch_rows = []
    missing = []
    duplicates = []
    seen = set()
    for mode in ["baseline", "columba"]:
        for iteration in args.iterations:
            for batch_id in expected_batches:
                bdir = run / "batches" / args.phase / mode / iteration / batch_id
                if not (bdir / "SUCCESS").exists():
                    missing.append(f"{mode}/{iteration}/{batch_id}")
                    continue
                gr = bdir / "guide_results.tsv"
                if not gr.exists():
                    missing.append(f"{mode}/{iteration}/{batch_id}/guide_results.tsv")
                    continue
                data = read_tsv(gr)
                batch_rows.append({
                    "mode": mode,
                    "iteration": iteration,
                    "batch_id": batch_id,
                    "guide_count": str(len(data)),
                    "exit_failures": str(sum(1 for r in data if r["exit_status"] != "0")),
                    "total_paf_records": str(sum(int(r["paf_records"]) for r in data)),
                    "total_wall_seconds": f"{sum(float(r['wall_seconds']) for r in data):.6f}",
                    "max_peak_rss_kib": str(max(int(r["peak_rss_kib"]) for r in data)),
                    "combined_stdout_sha256": hashlib.sha256("".join(r["stdout_sha256"] for r in data).encode()).hexdigest(),
                })
                for r in data:
                    key = (mode, iteration, r["guide_id"])
                    if key in seen:
                        duplicates.append("/".join(key))
                    seen.add(key)
                    rows.append(r)
    for mode in ["baseline", "columba"]:
        for iteration in args.iterations:
            got = {r["guide_id"] for r in rows if r["mode"] == mode and r["iteration"] == iteration}
            for gid in expected_guides:
                if gid not in got:
                    missing.append(f"{mode}/{iteration}/{gid}")
    summary_dir = run / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    with (summary_dir / "guide_results.tsv").open("w", newline="") as handle:
        fields = ["guide_id", "mode", "phase", "iteration", "exit_status", "wall_seconds", "user_seconds", "system_seconds", "peak_rss_kib", "paf_records", "stdout_sha256", "stderr_sha256"]
        w = csv.DictWriter(handle, fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    with (summary_dir / "batch_summary.tsv").open("w", newline="") as handle:
        fields = ["mode", "iteration", "batch_id", "guide_count", "exit_failures", "total_paf_records", "total_wall_seconds", "max_peak_rss_kib", "combined_stdout_sha256"]
        w = csv.DictWriter(handle, fields, delimiter="\t")
        w.writeheader(); w.writerows(batch_rows)
    with (summary_dir / "aggregate_status.tsv").open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"expected_guides\t{len(expected_guides)}\n")
        handle.write(f"expected_batches\t{len(expected_batches)}\n")
        handle.write(f"missing_items\t{len(missing)}\n")
        handle.write(f"duplicate_items\t{len(duplicates)}\n")
        handle.write(f"all_complete\t{str(not missing and not duplicates).lower()}\n")
    if missing:
        (summary_dir / "missing.txt").write_text("\n".join(missing) + "\n")
    if duplicates:
        (summary_dir / "duplicates.txt").write_text("\n".join(duplicates) + "\n")
    if missing or duplicates:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
