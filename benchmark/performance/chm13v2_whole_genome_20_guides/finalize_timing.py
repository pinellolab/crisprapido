#!/usr/bin/env python3
"""Finalize the controlled whole-genome timing campaign from raw batches."""

import argparse
import csv
import hashlib
import statistics
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parent
PANEL_SPECS = [
    ("baseline", "measured_1", "baseline_reference_scan"),
    ("columba", "measured_1", "cold_start_candidate"),
    ("columba", "measured_2", "warm_cache_candidate"),
    ("columba", "measured_3", "warm_cache_candidate"),
]
CONSISTENT_ENV_FIELDS = [
    "git_commit", "git_branch", "crisprapido_sha256", "columba_sha256",
    "guide_tsv_sha256", "reference", "columba_index", "pam",
    "max_mismatches", "max_bulges", "max_bulge_size", "min_match_fraction",
    "threads", "candidate_e", "timing_design", "hostname",
    "slurm_partition", "slurm_cpus_per_task",
]
CORRECTNESS = {
    "baseline_valid_loci": 127,
    "shared_baseline_loci": 127,
    "baseline_missing": 0,
    "columba_invalid": 0,
    "columba_valid_loci": 701,
}


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_environment(path):
    result = {}
    for line in path.read_text().splitlines():
        key, value = line.split("\t", 1)
        result[key] = value
    return result


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def combined_digest(rows, field):
    return hashlib.sha256(
        "".join(row[field] for row in rows).encode()
    ).hexdigest()


def fmt(value):
    return f"{value:.6f}"


def main():
    parser = argparse.ArgumentParser(
        description="Finalize one baseline and three Columba timing panels."
    )
    parser.add_argument("--run-root", required=True)
    args = parser.parse_args()
    run_root = Path(args.run_root).resolve()
    run_id = run_root.name

    guides = read_tsv(PACKAGE_DIR / "guides.tsv")
    batches = read_tsv(PACKAGE_DIR / "batches.tsv")
    expected_guides = [row["guide_id"] for row in guides]
    expected_batches = [row["batch_id"] for row in batches]
    expected_by_batch = {
        row["batch_id"]: row["guide_ids"].split(",") for row in batches
    }

    all_rows = []
    batch_rows = []
    provenance_rows = []
    panel_rows = []
    errors = []
    reference_environment = None
    panel_guide_rows = {}

    for mode, iteration, expected_cache_label in PANEL_SPECS:
        panel = []
        seen = set()
        for batch_id in expected_batches:
            batch_dir = (
                run_root / "batches" / "timing" / mode / iteration / batch_id
            )
            if not (batch_dir / "SUCCESS").exists():
                errors.append(f"missing SUCCESS: {mode}/{iteration}/{batch_id}")
                continue
            result_path = batch_dir / "guide_results.tsv"
            environment_path = batch_dir / "environment.tsv"
            if not result_path.exists() or not environment_path.exists():
                errors.append(
                    f"missing result or environment: {mode}/{iteration}/{batch_id}"
                )
                continue

            rows = read_tsv(result_path)
            environment = read_environment(environment_path)
            observed_ids = [row["guide_id"] for row in rows]
            if observed_ids != expected_by_batch[batch_id]:
                errors.append(
                    f"guide order mismatch: {mode}/{iteration}/{batch_id}"
                )
            expected_environment = {
                "run_id": run_id,
                "mode": mode,
                "phase": "timing",
                "iteration": iteration,
                "batch_id": batch_id,
                "cache_state_label": expected_cache_label,
            }
            for field, expected_value in expected_environment.items():
                if environment.get(field) != expected_value:
                    errors.append(
                        f"{field} mismatch: {mode}/{iteration}/{batch_id}"
                    )

            if reference_environment is None:
                reference_environment = environment
            for field in CONSISTENT_ENV_FIELDS:
                if environment.get(field) != reference_environment.get(field):
                    errors.append(
                        f"provenance mismatch for {field}: "
                        f"{mode}/{iteration}/{batch_id}"
                    )

            for row in rows:
                key = (mode, iteration, row["guide_id"])
                if key in seen:
                    errors.append(f"duplicate guide: {'/'.join(key)}")
                seen.add(key)
                if row["mode"] != mode or row["iteration"] != iteration:
                    errors.append(f"result metadata mismatch: {'/'.join(key)}")
                if row["phase"] != "timing" or row["exit_status"] != "0":
                    errors.append(f"failed result: {'/'.join(key)}")
                panel.append(row)
                all_rows.append(row)

            batch_rows.append({
                "mode": mode,
                "iteration": iteration,
                "cache_state_label": expected_cache_label,
                "batch_id": batch_id,
                "guide_count": len(rows),
                "exit_failures": sum(row["exit_status"] != "0" for row in rows),
                "total_paf_records": sum(int(row["paf_records"]) for row in rows),
                "total_wall_seconds": fmt(
                    sum(float(row["wall_seconds"]) for row in rows)
                ),
                "total_user_seconds": fmt(
                    sum(float(row["user_seconds"]) for row in rows)
                ),
                "total_system_seconds": fmt(
                    sum(float(row["system_seconds"]) for row in rows)
                ),
                "max_peak_rss_kib": max(
                    (int(row["peak_rss_kib"]) for row in rows), default=0
                ),
                "combined_stdout_sha256": combined_digest(rows, "stdout_sha256"),
                "combined_stderr_sha256": combined_digest(rows, "stderr_sha256"),
            })
            provenance_rows.append({
                "mode": mode,
                "iteration": iteration,
                "cache_state_label": environment.get("cache_state_label", "NA"),
                "batch_id": batch_id,
                "slurm_job_id": environment.get("slurm_job_id", "NA"),
                "slurm_array_job_id": environment.get(
                    "slurm_array_job_id", "NA"
                ),
                "slurm_array_task_id": environment.get(
                    "slurm_array_task_id", "NA"
                ),
                "slurm_job_dependency": environment.get(
                    "slurm_job_dependency", "NA"
                ),
                "hostname": environment.get("hostname", "NA"),
                "partition": environment.get("slurm_partition", "NA"),
                "date_utc": environment.get("date_utc", "NA"),
                "git_commit": environment.get("git_commit", "NA"),
                "git_branch": environment.get("git_branch", "NA"),
                "crisprapido_sha256": environment.get(
                    "crisprapido_sha256", "NA"
                ),
                "columba_sha256": environment.get("columba_sha256", "NA"),
                "guide_tsv_sha256": environment.get("guide_tsv_sha256", "NA"),
                "reference": environment.get("reference", "NA"),
                "columba_index": environment.get("columba_index", "NA"),
                "pam": environment.get("pam", "NA"),
                "max_mismatches": environment.get("max_mismatches", "NA"),
                "max_bulges": environment.get("max_bulges", "NA"),
                "max_bulge_size": environment.get("max_bulge_size", "NA"),
                "min_match_fraction": environment.get(
                    "min_match_fraction", "NA"
                ),
                "threads": environment.get("threads", "NA"),
                "candidate_e": environment.get("candidate_e", "NA"),
            })

        if [row["guide_id"] for row in panel] != expected_guides:
            errors.append(f"panel guide order mismatch: {mode}/{iteration}")
        panel_guide_rows[(mode, iteration)] = {
            row["guide_id"]: row for row in panel
        }
        panel_rows.append({
            "mode": mode,
            "iteration": iteration,
            "cache_state_label": expected_cache_label,
            "guide_count": len(panel),
            "exit_failures": sum(row["exit_status"] != "0" for row in panel),
            "total_paf_records": sum(int(row["paf_records"]) for row in panel),
            "sum_wall_seconds": fmt(
                sum(float(row["wall_seconds"]) for row in panel)
            ),
            "sum_user_seconds": fmt(
                sum(float(row["user_seconds"]) for row in panel)
            ),
            "sum_system_seconds": fmt(
                sum(float(row["system_seconds"]) for row in panel)
            ),
            "max_peak_rss_kib": max(
                (int(row["peak_rss_kib"]) for row in panel), default=0
            ),
            "panel_stdout_sha256": combined_digest(panel, "stdout_sha256"),
            "panel_stderr_sha256": combined_digest(panel, "stderr_sha256"),
        })

    correctness_rows = read_tsv(PACKAGE_DIR / "guide_results.tsv")
    correctness_by_mode = {
        mode: {
            row["guide_id"]: row
            for row in correctness_rows
            if row["mode"] == mode
        }
        for mode in ("baseline", "columba")
    }
    matches_correctness = {}
    for mode, iteration, _ in PANEL_SPECS:
        timing = panel_guide_rows.get((mode, iteration), {})
        expected = correctness_by_mode[mode]
        matches_correctness[(mode, iteration)] = (
            set(timing) == set(expected)
            and all(
                timing[guide_id]["stdout_sha256"]
                == expected[guide_id]["stdout_sha256"]
                for guide_id in timing
            )
        )

    columba_maps = [
        panel_guide_rows[("columba", iteration)]
        for iteration in ("measured_1", "measured_2", "measured_3")
    ]
    columba_stdout_deterministic = all(
        columba_maps[0][guide_id]["stdout_sha256"]
        == columba_maps[1][guide_id]["stdout_sha256"]
        == columba_maps[2][guide_id]["stdout_sha256"]
        for guide_id in expected_guides
    )
    columba_stderr_deterministic = all(
        columba_maps[0][guide_id]["stderr_sha256"]
        == columba_maps[1][guide_id]["stderr_sha256"]
        == columba_maps[2][guide_id]["stderr_sha256"]
        for guide_id in expected_guides
    )

    panels = {(row["mode"], row["iteration"]): row for row in panel_rows}
    baseline = panels[("baseline", "measured_1")]
    cold = panels[("columba", "measured_1")]
    warm = [
        panels[("columba", "measured_2")],
        panels[("columba", "measured_3")],
    ]
    warm_walls = [float(row["sum_wall_seconds"]) for row in warm]
    warm_users = [float(row["sum_user_seconds"]) for row in warm]
    warm_systems = [float(row["sum_system_seconds"]) for row in warm]
    warm_rss = [int(row["max_peak_rss_kib"]) for row in warm]
    all_columba_walls = [
        float(panels[("columba", iteration)]["sum_wall_seconds"])
        for iteration in ("measured_1", "measured_2", "measured_3")
    ]
    baseline_wall = float(baseline["sum_wall_seconds"])
    cold_wall = float(cold["sum_wall_seconds"])
    warm_wall = statistics.median(warm_walls)
    warm_rss_median = int(statistics.median(warm_rss))

    status_rows = [
        {"metric": "run_id", "value": run_id},
        {"metric": "expected_batches_per_panel", "value": len(expected_batches)},
        {"metric": "expected_guides_per_panel", "value": len(expected_guides)},
        {"metric": "panels", "value": len(PANEL_SPECS)},
        {"metric": "missing_or_duplicate_items", "value": len(errors)},
        {
            "metric": "all_failures_zero",
            "value": str(
                all(row["exit_failures"] == 0 for row in panel_rows)
            ).lower(),
        },
        {
            "metric": "columba_output_deterministic",
            "value": str(columba_stdout_deterministic).lower(),
        },
        {
            "metric": "columba_stderr_deterministic",
            "value": str(columba_stderr_deterministic).lower(),
        },
        {
            "metric": "provenance_consistent",
            "value": str(
                not any("provenance mismatch" in error for error in errors)
            ).lower(),
        },
        {
            "metric": "baseline_array_job_ids",
            "value": ",".join(sorted({
                row["slurm_array_job_id"]
                for row in provenance_rows
                if row["mode"] == "baseline"
            })),
        },
        {"metric": "all_complete", "value": str(not errors).lower()},
    ]

    timing_summary = {
        "run_id": run_id,
        "timing_node": reference_environment["hostname"],
        "timing_partition": reference_environment["slurm_partition"],
        "baseline_replicates": 1,
        "baseline_wall_seconds": baseline["sum_wall_seconds"],
        "baseline_seconds_per_guide": fmt(
            baseline_wall / len(expected_guides)
        ),
        "baseline_user_seconds": baseline["sum_user_seconds"],
        "baseline_system_seconds": baseline["sum_system_seconds"],
        "baseline_peak_rss_kib": baseline["max_peak_rss_kib"],
        "baseline_total_paf_records": baseline["total_paf_records"],
        "baseline_matches_correctness_output": str(
            matches_correctness[("baseline", "measured_1")]
        ).lower(),
        "columba_cold_iteration": "measured_1",
        "columba_cold_wall_seconds": cold["sum_wall_seconds"],
        "columba_cold_seconds_per_guide": fmt(
            cold_wall / len(expected_guides)
        ),
        "columba_cold_user_seconds": cold["sum_user_seconds"],
        "columba_cold_system_seconds": cold["sum_system_seconds"],
        "columba_cold_peak_rss_kib": cold["max_peak_rss_kib"],
        "columba_warm_replicates": 2,
        "columba_warm_iterations": "measured_2,measured_3",
        "columba_warm_wall_seconds": ",".join(
            fmt(value) for value in warm_walls
        ),
        "columba_warm_median_wall_seconds": fmt(warm_wall),
        "columba_warm_seconds_per_guide": fmt(
            warm_wall / len(expected_guides)
        ),
        "columba_warm_min_wall_seconds": fmt(min(warm_walls)),
        "columba_warm_max_wall_seconds": fmt(max(warm_walls)),
        "columba_warm_range_seconds": fmt(
            max(warm_walls) - min(warm_walls)
        ),
        "columba_warm_median_user_seconds": fmt(
            statistics.median(warm_users)
        ),
        "columba_warm_median_system_seconds": fmt(
            statistics.median(warm_systems)
        ),
        "columba_warm_median_peak_rss_kib": warm_rss_median,
        "columba_all_three_median_wall_seconds": fmt(
            statistics.median(all_columba_walls)
        ),
        "columba_total_paf_records": cold["total_paf_records"],
        "columba_output_deterministic": str(
            columba_stdout_deterministic
        ).lower(),
        "columba_stderr_deterministic": str(
            columba_stderr_deterministic
        ).lower(),
        "columba_matches_correctness_output": str(all(
            matches_correctness[("columba", iteration)]
            for iteration in ("measured_1", "measured_2", "measured_3")
        )).lower(),
        "cold_start_speedup": fmt(baseline_wall / cold_wall),
        "steady_state_speedup": fmt(baseline_wall / warm_wall),
        "warm_memory_ratio": fmt(
            warm_rss_median / int(baseline["max_peak_rss_kib"])
        ),
        **CORRECTNESS,
        "timing_complete": str(not errors).lower(),
    }

    guide_fields = [
        "guide_id", "mode", "phase", "iteration", "exit_status",
        "wall_seconds", "user_seconds", "system_seconds", "peak_rss_kib",
        "paf_records", "stdout_sha256", "stderr_sha256",
    ]
    write_tsv(PACKAGE_DIR / "timing_runs.tsv", all_rows, guide_fields)
    write_tsv(
        PACKAGE_DIR / "batch_timing_summary.tsv",
        batch_rows,
        list(batch_rows[0]),
    )
    write_tsv(PACKAGE_DIR / "run_summary.tsv", panel_rows, list(panel_rows[0]))
    write_tsv(
        PACKAGE_DIR / "timing_provenance.tsv",
        provenance_rows,
        list(provenance_rows[0]),
    )
    write_tsv(
        PACKAGE_DIR / "timing_aggregate_status.tsv",
        status_rows,
        ["metric", "value"],
    )
    write_tsv(
        PACKAGE_DIR / "timing_summary.tsv",
        [timing_summary],
        list(timing_summary),
    )

    scaling_path = PACKAGE_DIR / "scaling_summary.tsv"
    scaling_rows = read_tsv(scaling_path)
    whole_row = {
        "reference": "CHM13v2_whole_genome",
        "reference_length_bp": 3117292070,
        "guide_count": len(expected_guides),
        "timing_node": reference_environment["hostname"],
        "baseline_wall_seconds": baseline["sum_wall_seconds"],
        "columba_wall_seconds": fmt(warm_wall),
        "baseline_seconds_per_guide": fmt(
            baseline_wall / len(expected_guides)
        ),
        "columba_seconds_per_guide": fmt(
            warm_wall / len(expected_guides)
        ),
        "observed_speedup": fmt(baseline_wall / warm_wall),
        "baseline_peak_rss_kib": baseline["max_peak_rss_kib"],
        "columba_peak_rss_kib": warm_rss_median,
        "memory_ratio": fmt(
            warm_rss_median / int(baseline["max_peak_rss_kib"])
        ),
        "baseline_paf_records": baseline["total_paf_records"],
        "columba_paf_records": cold["total_paf_records"],
        "baseline_valid_loci": CORRECTNESS["baseline_valid_loci"],
        "columba_valid_loci": CORRECTNESS["columba_valid_loci"],
        "baseline_missing": CORRECTNESS["baseline_missing"],
        "columba_invalid": CORRECTNESS["columba_invalid"],
        "baseline_replicates": 1,
        "columba_replicates": 2,
    }
    scaling_rows = [
        row
        for row in scaling_rows
        if row["reference"] != "CHM13v2_whole_genome"
    ]
    scaling_rows.append(whole_row)
    write_tsv(scaling_path, scaling_rows, list(scaling_rows[0]))

    manifest_rows = []
    timing_root = run_root / "batches" / "timing"
    for path in sorted(
        item for item in timing_root.rglob("*") if item.is_file()
    ):
        manifest_rows.append({
            "relative_path": path.relative_to(run_root),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    write_tsv(
        PACKAGE_DIR / "raw_timing_manifest.tsv",
        manifest_rows,
        ["relative_path", "size_bytes", "sha256"],
    )

    if errors:
        for error in errors:
            print(error)
        raise SystemExit(2)
    print(f"finalized {run_id}")
    print(f"baseline_wall_seconds={baseline['sum_wall_seconds']}")
    print(f"columba_cold_wall_seconds={cold['sum_wall_seconds']}")
    print(f"columba_warm_median_wall_seconds={fmt(warm_wall)}")
    print(f"steady_state_speedup={fmt(baseline_wall / warm_wall)}")


if __name__ == "__main__":
    main()
