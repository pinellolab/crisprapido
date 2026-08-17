#!/usr/bin/env python3
import argparse
import csv
import hashlib
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description="Aggregate complete whole-genome benchmark batches.")
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--phase", default="correctness")
    parser.add_argument("--iterations", nargs="+", default=["pilot"])
    args = parser.parse_args()
    run_root = Path(args.run_root).resolve()
    guides = read_tsv(SCRIPT_DIR / "guides.tsv")
    batches = read_tsv(SCRIPT_DIR / "batches.tsv")
    expected_guides = [row["guide_id"] for row in guides]
    expected_batches = [row["batch_id"] for row in batches]
    rows = []
    batch_rows = []
    missing = []
    duplicates = []
    seen = set()
    for mode in ["baseline", "columba"]:
        for iteration in args.iterations:
            for batch_id in expected_batches:
                batch_dir = run_root / "batches" / args.phase / mode / iteration / batch_id
                if not (batch_dir / "SUCCESS").exists():
                    missing.append(f"{mode}/{iteration}/{batch_id}")
                    continue
                result_path = batch_dir / "guide_results.tsv"
                if not result_path.exists():
                    missing.append(f"{mode}/{iteration}/{batch_id}/guide_results.tsv")
                    continue
                data = read_tsv(result_path)
                batch_rows.append(
                    {
                        "mode": mode,
                        "iteration": iteration,
                        "batch_id": batch_id,
                        "guide_count": len(data),
                        "exit_failures": sum(row["exit_status"] != "0" for row in data),
                        "total_paf_records": sum(int(row["paf_records"]) for row in data),
                        "total_wall_seconds": f"{sum(float(row['wall_seconds']) for row in data):.6f}",
                        "total_user_seconds": f"{sum(float(row['user_seconds']) for row in data):.6f}",
                        "total_system_seconds": f"{sum(float(row['system_seconds']) for row in data):.6f}",
                        "max_peak_rss_kib": max(int(row["peak_rss_kib"]) for row in data),
                        "combined_stdout_sha256": hashlib.sha256(
                            "".join(row["stdout_sha256"] for row in data).encode()
                        ).hexdigest(),
                    }
                )
                for row in data:
                    key = (mode, iteration, row["guide_id"])
                    if key in seen:
                        duplicates.append("/".join(key))
                    seen.add(key)
                    rows.append(row)
    for mode in ["baseline", "columba"]:
        for iteration in args.iterations:
            observed = {
                row["guide_id"]
                for row in rows
                if row["mode"] == mode and row["iteration"] == iteration
            }
            missing.extend(
                f"{mode}/{iteration}/{guide_id}"
                for guide_id in expected_guides
                if guide_id not in observed
            )

    summary = run_root / "summary"
    summary.mkdir(parents=True, exist_ok=True)
    guide_fields = [
        "guide_id", "mode", "phase", "iteration", "exit_status", "wall_seconds",
        "user_seconds", "system_seconds", "peak_rss_kib", "paf_records",
        "stdout_sha256", "stderr_sha256",
    ]
    with (summary / "guide_results.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, guide_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    batch_fields = list(batch_rows[0]) if batch_rows else ["mode", "iteration", "batch_id"]
    with (summary / "batch_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, batch_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(batch_rows)
    with (summary / "aggregate_status.tsv").open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"expected_guides\t{len(expected_guides)}\n")
        handle.write(f"expected_batches\t{len(expected_batches)}\n")
        handle.write(f"missing_items\t{len(missing)}\n")
        handle.write(f"duplicate_items\t{len(duplicates)}\n")
        handle.write(f"all_complete\t{str(not missing and not duplicates).lower()}\n")
    if missing:
        (summary / "missing.txt").write_text("\n".join(missing) + "\n")
    if duplicates:
        (summary / "duplicates.txt").write_text("\n".join(duplicates) + "\n")
    if missing or duplicates:
        raise SystemExit(2)


if __name__ == "__main__":
    main()

