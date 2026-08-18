# CHM13 chr22 100-Guide Performance Benchmark

This package extends the 20-guide chr22 pilot to a deterministic 100-guide panel. It compares original CRISPRapido sliding-window candidate generation with automatic Columba candidate generation followed by anchored WFA2 verification, filtering, PAM extraction, CFD scoring, and PAF reporting.

The recorded result uses the committed bounded-fallback implementation:

- full-guide-aware WFA2 verification
- coordinate-anchored imported Columba candidate verification
- Columba candidate bound `candidate_e = m + b*z`
- bounded fallback expansion for Columba-suppressed insertion/bulge candidates

## Guide Panel

Reference: `../data/real_reference/chm13v2_chr22.fa`  
Columba index prefix: `../results/chm13_chr22_index/chm13v2_chr22`  
Chromosome/header: `22`  
Reference length: `51324926` bp

Selection rules:

- Forward-strand 20-nt protospacers only.
- The two reference bases immediately following the protospacer must be `GG`.
- Guide and PAM must contain only `A/C/G/T`.
- Guides with homopolymer runs longer than 5 bases are excluded.
- Exact duplicate guide sequences are excluded before interval selection.
- The chromosome is divided into 100 approximately equal intervals.
- One candidate nearest the interval center is selected per interval.
- Ties are broken by lowest zero-based coordinate.

Preparation result:

```text
candidate_positions=3264489
unique_candidate_guides=2368140
selected_guides=100
unique_guides=74
low_copy_guides=21
repetitive_guides=5
no_exact_hit_guides=0
all_columba_validated=True
all_crisprapido_validated=True
```

Direct validation at edit distance 0 found `2322` total Columba `20M`/`NM:i:0` mapped records across the 100 guides, and CRISPRapido automatic Columba exact-match validation produced `2322` total PAF records. All 100 guide validation rows passed.

## Benchmark Parameters

- PAM: `GG`
- maximum mismatches: `0`
- maximum bulges: `1`
- maximum bulge size: `2`
- Columba candidate edit-distance bound: `candidate_e = m + b*z = 2`
- minimum match fraction: `0.75`
- threads: `1`

Correctness uses conservative biological-locus equivalence, not raw PAF byte identity. This is necessary because baseline sliding-window mode reports one WFA2-selected result per search window, while Columba mode enumerates candidate loci up to the configured edit-distance bound and verifies them at the reported coordinate.

## Recorded Run

Recorded run: `raw/pilot_100_fallback_20260801T180205Z/`

All baseline and Columba commands exited successfully. Correctness eligibility passed, so the harness ran one warm-up panel per mode and three measured panels per mode.

Correctness summary:

```text
baseline_raw_records: 2892
columba_raw_records: 12598
baseline_valid_loci: 690
columba_valid_loci: 8896
shared_baseline_loci: 690
baseline_missing_from_columba: 0
columba_only_valid_loci: 8206
columba_invalid_records: 0
columba_non_gg_pam_records: 3702
timing_eligible: yes
```

The bounded fallback recovered the known Columba-suppressed loci exposed by the initial 100-guide comparison. Fallback workload counters are not exposed by the CLI, and automatic-mode temporary SAM files are deleted, so exact fallback-trigger and expansion counts are not recorded in this package.

## Timing Result

Measured panel runs, one thread:

```text
baseline wall seconds: 3424.420954, 3425.267860, 3419.766600
Columba wall seconds: 64.896516, 64.244734, 63.991992
baseline median wall: 3424.420954 s
Columba median wall: 64.244734 s
observed speedup in this pilot: 53.302749x
baseline median peak RSS: 125116 KiB
Columba median peak RSS: 293164 KiB
memory ratio: 2.343138
```

Deterministic output checks passed within each mode: the combined stdout checksum was identical across all three measured baseline panels and identical across all three measured Columba panels.

## Scaling From 20 Guides

The 20-guide chr22 pilot reported:

```text
baseline median wall: 664.124025 s
Columba median wall: 9.099324 s
baseline median RSS: 124276 KiB
Columba median RSS: 290436 KiB
baseline raw records: 310
Columba raw records: 3095
```

Per-guide median wall time:

```text
20-guide baseline: 33.206201 s/guide
100-guide baseline: 34.244210 s/guide
20-guide Columba: 0.454966 s/guide
100-guide Columba: 0.642447 s/guide
```

Baseline scaling is close to linear with guide count. Columba remains much faster, but per-guide time increases on the 100-guide panel because the panel contains higher-output guides and the bounded fallback verifies additional nearby spans for suppressed insertion candidates. Output volume also increased from `3095` Columba PAF rows in the 20-guide pilot to `12598` rows here.

## Raw Outputs

Raw performance outputs are written under `raw/<RUN_ID>/` and excluded from Git. Preparation validation raw outputs are written under `preparation_raw/` and excluded from Git. Compact summaries and checksum manifests are retained:

- `raw_manifest.tsv`
- `preparation_raw_manifest.tsv`

The raw manifest records relative path, file size, and SHA-256 checksum for the retained raw run.

## Limitations

- This is a 100-guide chr22 pilot, not a whole-genome benchmark.
- Runs used one thread only.
- The Linux filesystem cache was not cleared.
- `/proc` RSS sampling can miss brief memory peaks between samples.
- Candidate and fallback workload counters are not exposed by the CLI.
- Columba output volume is larger because candidate generation is more exhaustive than baseline sliding windows.
- `3702` independently valid Columba records had non-`GG` PAM and PAM filtering remains a separate issue from candidate-generation performance.

## Recommendation

The 100-guide chr22 benchmark is timing-eligible and shows stable baseline-locus recovery with zero independently invalid Columba alignments. The next useful step is multithreaded chr22 testing, followed by chr1. Whole-genome benchmarking should wait until multithreaded behavior and output-volume management are characterized.

## Matched-Node tux05 Timing Rerun Plan

A timing-only Slurm harness is available to rerun this same validated 100-guide chr22 panel on partition `tux`, node `tux05`, matching the chr2_100 timing hardware. This rerun does not regenerate guides and does not repeat correctness unless a timing run reveals a provenance or execution problem.

The rerun uses the same biological parameters as the finalized chr2_100 benchmark: `PAM=GG`, `m=0`, `b=1`, `z=2`, `f=0.75`, and `threads=1`. The panel remains `chr22_guides.tsv`, the reference remains `../data/real_reference/chm13v2_chr22.fa`, and the Columba index prefix remains `../results/chm13_chr22_index/chm13v2_chr22`.

Slurm timing scripts use `SLURM_SUBMIT_DIR` as the repository root and write raw outputs under ignored timestamped directories in `raw/`. The prepared design is one measured baseline panel and three measured Columba panels. Because only one baseline replicate is planned, report its wall time as a single measured pass, not a median.

Use `submission_commands.txt` or `DRY_RUN=1 ./submit_plan.sh` to print the smoke and full timing commands. Every prepared command includes `-p tux -w tux05`.

## Matched-Node tux05 Timing Result

The timing-only rerun on Slurm partition `tux`, node `tux05`, completed as RUN_ID `chr22_100_guides_tux05_20260813T160700Z_timing`. It reused the same 100-guide panel, CHM13v2 chr22 reference, Columba index, production release binary, and biological parameters as the finalized chr2_100 benchmark: `PAM=GG`, `m=0`, `b=1`, `z=2`, `f=0.75`, and `threads=1`.

Completeness checks passed: baseline `measured_1` completed 10/10 batches, and Columba `measured_1`, `measured_2`, and `measured_3` each completed 10/10 batches. Each full panel contained 100 guides, with no missing or duplicate guide assignments and no failed guide runs.

Matched-node timing summary:

```text
baseline measured wall time: 1932.584993 s
baseline user CPU time: 1921.051972 s
baseline system CPU time: 2.355807 s
baseline peak RSS: 127584 KiB
baseline PAF records: 2892
Columba measured wall times: 84.294101, 84.993735, 83.573807 s
Columba median wall time: 84.294101 s
Columba wall-time min/max: 83.573807 / 84.993735 s
Columba median peak RSS: 292428 KiB
Columba PAF records: 12608
Columba output deterministic across measured runs: true
observed speedup using one baseline replicate and the median of three Columba replicates: 22.926693x
memory ratio, Columba median peak RSS over baseline measured peak RSS: 2.292043x
```

The old non-tux05 chr22_100 timing reported baseline `3424.420954` s, Columba `64.244734` s, and speedup `53.302749x`. On `tux05`, baseline was `1.771938x` faster than the old-node baseline, while Columba was `1.312078x` slower than the old-node Columba run. This is why the chr2 comparison should use the matched-node values rather than the older chr22_100 timing.

Matched-node comparison with chr2_100 on `tux05`:

```text
chr22 length: 51,324,926 bp
chr2 length: 242,696,752 bp
reference length ratio chr2/chr22: 4.728633x
baseline runtime ratio chr2/chr22: 4.721868x
Columba runtime ratio chr2/chr22: 3.304718x
baseline RSS ratio chr2/chr22: 4.515331x
Columba RSS ratio chr2/chr22: 4.056000x
```

With hardware matched, baseline runtime scales almost exactly with reference length for this 100-guide panel. Columba runtime also increases with reference size, but less than linearly relative to chromosome length in this comparison. RSS increases substantially with the larger chr2 reference and index; the observed memory ratios are lower than the reference length ratio but still close enough to show reference-size dependence.
