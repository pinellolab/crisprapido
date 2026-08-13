#!/usr/bin/env python3
import argparse
import csv
import hashlib
import statistics
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
THREADS = [1, 2, 4, 8]
DESIGN = [("baseline", "measured_1"), ("columba", "measured_1"), ("columba", "measured_2"), ("columba", "measured_3")]

def read_tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-root", default=str(SCRIPT_DIR / "raw"), help="Raw directory containing per-thread RUN_ID directories")
    p.add_argument("--run-prefix", required=True, help="Prefix before _t<threads>_timing")
    return p.parse_args()

def main():
    args = parse_args()
    raw = Path(args.raw_root)
    batches = read_tsv(SCRIPT_DIR / "batches.tsv")
    guide_panel = read_tsv(SCRIPT_DIR.parent / "chr22_100_guides" / "chr22_guides.tsv")
    expected_batches = [b["batch_id"] for b in batches]
    expected_guides = [g["guide_id"] for g in guide_panel]
    guide_rows = []
    panel_rows = []
    problems = []
    for threads in THREADS:
        run = raw / f"{args.run_prefix}_t{threads}_timing"
        for mode, iteration in DESIGN:
            rows = []
            batch_rows = []
            seen = set()
            for batch_id in expected_batches:
                bdir = run / f"threads_{threads}" / "batches" / "timing" / mode / iteration / batch_id
                if not (bdir / "SUCCESS").exists():
                    problems.append(f"missing {threads}/{mode}/{iteration}/{batch_id}/SUCCESS")
                    continue
                gr = bdir / "guide_results.tsv"
                if not gr.exists():
                    problems.append(f"missing {threads}/{mode}/{iteration}/{batch_id}/guide_results.tsv")
                    continue
                data = read_tsv(gr)
                for row in data:
                    key = row["guide_id"]
                    if key in seen:
                        problems.append(f"duplicate {threads}/{mode}/{iteration}/{key}")
                    seen.add(key)
                rows.extend(data)
                batch_rows.append({
                    "threads": str(threads), "mode": mode, "iteration": iteration, "batch_id": batch_id,
                    "guide_count": str(len(data)),
                    "exit_failures": str(sum(1 for r in data if r["exit_status"] != "0")),
                    "total_paf_records": str(sum(int(r["paf_records"]) for r in data)),
                    "total_wall_seconds": f"{sum(float(r['wall_seconds']) for r in data):.6f}",
                    "total_user_seconds": f"{sum(float(r['user_seconds']) for r in data):.6f}",
                    "total_system_seconds": f"{sum(float(r['system_seconds']) for r in data):.6f}",
                    "max_peak_rss_kib": str(max(int(r["peak_rss_kib"]) for r in data)),
                    "combined_stdout_sha256": hashlib.sha256("".join(r["stdout_sha256"] for r in data).encode()).hexdigest(),
                    "combined_stderr_sha256": hashlib.sha256("".join(r["stderr_sha256"] for r in data).encode()).hexdigest(),
                })
            if seen != set(expected_guides):
                problems.append(f"guide set mismatch {threads}/{mode}/{iteration}")
            guide_rows.extend(rows)
            panel_rows.append({
                "threads": str(threads), "mode": mode, "iteration": iteration,
                "guide_count": str(sum(int(r["guide_count"]) for r in batch_rows)),
                "exit_failures": str(sum(int(r["exit_failures"]) for r in batch_rows)),
                "total_paf_records": str(sum(int(r["total_paf_records"]) for r in batch_rows)),
                "sum_wall_seconds": f"{sum(float(r['total_wall_seconds']) for r in batch_rows):.6f}",
                "sum_user_seconds": f"{sum(float(r['total_user_seconds']) for r in batch_rows):.6f}",
                "sum_system_seconds": f"{sum(float(r['total_system_seconds']) for r in batch_rows):.6f}",
                "max_peak_rss_kib": str(max(int(r["max_peak_rss_kib"]) for r in batch_rows)),
                "panel_stdout_sha256": hashlib.sha256("".join(r["combined_stdout_sha256"] for r in batch_rows).encode()).hexdigest(),
                "panel_stderr_sha256": hashlib.sha256("".join(r["combined_stderr_sha256"] for r in batch_rows).encode()).hexdigest(),
            })
    with (SCRIPT_DIR / "thread_guide_results.tsv").open("w", newline="") as f:
        fields = ["guide_id","mode","phase","iteration","threads","exit_status","wall_seconds","user_seconds","system_seconds","peak_rss_kib","paf_records","stdout_sha256","stderr_sha256"]
        w = csv.DictWriter(f, fields, delimiter="\t"); w.writeheader(); w.writerows(guide_rows)
    with (SCRIPT_DIR / "thread_run_summary.tsv").open("w", newline="") as f:
        fields = ["threads","mode","iteration","guide_count","exit_failures","total_paf_records","sum_wall_seconds","sum_user_seconds","sum_system_seconds","max_peak_rss_kib","panel_stdout_sha256","panel_stderr_sha256"]
        w = csv.DictWriter(f, fields, delimiter="\t"); w.writeheader(); w.writerows(panel_rows)
    by = {(int(r["threads"]), r["mode"], r["iteration"]): r for r in panel_rows}
    base1 = float(by[(1,"baseline","measured_1")]["sum_wall_seconds"]) if (1,"baseline","measured_1") in by else None
    col1 = statistics.median(float(by[(1,"columba",it)]["sum_wall_seconds"]) for it in ["measured_1","measured_2","measured_3"]) if all((1,"columba",it) in by for it in ["measured_1","measured_2","measured_3"]) else None
    summary = []
    for threads in THREADS:
        b = by.get((threads,"baseline","measured_1"))
        crows = [by.get((threads,"columba",it)) for it in ["measured_1","measured_2","measured_3"]]
        if b:
            wall=float(b["sum_wall_seconds"]); speed=(base1/wall) if base1 else 1.0
            summary.append({"reference":"chr22","guide_count":"100","threads":str(threads),"mode":"baseline","wall_seconds":f"{wall:.6f}","seconds_per_guide":f"{wall/100:.6f}","speedup_vs_1thread":f"{speed:.6f}","parallel_efficiency":f"{speed/threads:.6f}","peak_rss_kib":b["max_peak_rss_kib"],"paf_records":b["total_paf_records"],"replicates":"1","node":"tux05"})
        if all(crows):
            walls=[float(r["sum_wall_seconds"]) for r in crows]; rss=[int(r["max_peak_rss_kib"]) for r in crows]
            wall=statistics.median(walls); speed=(col1/wall) if col1 else 1.0
            summary.append({"reference":"chr22","guide_count":"100","threads":str(threads),"mode":"columba","wall_seconds":f"{wall:.6f}","seconds_per_guide":f"{wall/100:.6f}","speedup_vs_1thread":f"{speed:.6f}","parallel_efficiency":f"{speed/threads:.6f}","peak_rss_kib":str(int(statistics.median(rss))),"paf_records":crows[0]["total_paf_records"],"replicates":"3","node":"tux05"})
    with (SCRIPT_DIR / "thread_scaling_summary.tsv").open("w", newline="") as f:
        fields = ["reference","guide_count","threads","mode","wall_seconds","seconds_per_guide","speedup_vs_1thread","parallel_efficiency","peak_rss_kib","paf_records","replicates","node"]
        w = csv.DictWriter(f, fields, delimiter="\t"); w.writeheader(); w.writerows(summary)
    with (SCRIPT_DIR / "thread_aggregate_status.tsv").open("w") as f:
        f.write("metric\tvalue\n")
        f.write(f"run_prefix\t{args.run_prefix}\n")
        f.write(f"problems\t{len(problems)}\n")
        f.write(f"all_complete\t{str(not problems).lower()}\n")
    if problems:
        (SCRIPT_DIR / "thread_aggregate_problems.txt").write_text("\n".join(problems)+"\n")
        raise SystemExit(2)

if __name__ == "__main__":
    main()
