# CHM13 chr22 Pilot Performance Benchmark

This benchmark compares CRISPRapido's original sliding-window candidate generation against automatic Columba candidate generation on CHM13 chromosome 22. Both modes use the same release binary after candidate generation: WFA2 verification, hit filtering, PAM extraction, CFD scoring, and PAF reporting are shared.

The current CRISPRapido CLI accepts one guide sequence per invocation, so the harness runs the deterministic 20-guide chr22 panel one guide at a time and aggregates per-guide measurements.

Default inputs, relative to the repository root:

- Reference: `../data/real_reference/chm13v2_chr22.fa`
- Columba index prefix: `../results/chm13_chr22_index/chm13v2_chr22`
- Guide panel: `../results/chm13_chr22_guide_preparation/chr22_guides.tsv`
- CRISPRapido binary: `target/release/crisprapido`
- Columba binary: `../columba/build_Vanilla/columba`

Default biological and computational parameters:

- PAM: `GG`
- maximum mismatches: `0`
- maximum bulges: `1`
- maximum bulge size: `2`
- Columba candidate edit-distance bound: `candidate_e = m + b*z = 2`
- minimum match fraction: `0.75`
- threads: `1`

The baseline command scans the full chr22 FASTA with CRISPRapido's normal sliding-window implementation. The Columba command invokes Columba with the existing chr22 index and then verifies each imported candidate at its SAM-reported coordinate. Columba index construction time is excluded from query runtime.

Manual byte-identical output is not expected for this benchmark. Baseline mode reports one WFA2-selected result per sliding window/strand and may omit valid loci because of windowing or traceback selection. Columba mode exhaustively generates candidates up to `candidate_e=2`, then coordinate-anchored WFA2 verifies each locus. Correctness eligibility is based on a conservative biological-locus comparison: every independently valid baseline locus must be recovered by Columba, allowing small coordinate/indel-placement equivalence, and Columba must not emit independently invalid alignments.

Memory measurement is implemented by `collect_peak_rss.py`, which launches each process, samples aggregate `VmRSS` for the process and descendants through `/proc`, and records the peak. The sampler also records wall time with a monotonic clock and child user/system CPU time through Python's `resource` module. Very short-lived memory peaks between samples may be missed.

Run:

```bash
benchmark/performance/chr22/run_chr22_pilot.sh
```

Useful overrides:

```bash
THREADS=1 MAX_MISMATCHES=0 MAX_BULGES=1 MAX_BULGE_SIZE=2 RUN_ID=my_run benchmark/performance/chr22/run_chr22_pilot.sh
```

Raw outputs are written under `raw/<RUN_ID>/`. Top-level TSV files mirror the latest completed run for convenient inspection.

## Final Pilot Run

Recorded run: `raw/pilot_final_20260731T163824Z/`

Correctness eligibility passed:

```text
baseline_raw_records: 310
columba_raw_records: 3095
baseline_valid_loci: 291
columba_valid_loci: 3036
shared_baseline_loci: 291
baseline_missing_from_columba: 0
columba_only_valid_loci: 2745
columba_invalid_records: 0
columba_non_gg_pam_records: 59
timing_eligible: yes
```

The 59 non-`GG` PAM Columba rows are alignment-valid under the full-guide oracle but are not counted as PAM-valid biological loci for the eligibility comparison.

Timing used one warm-up panel per mode and three measured panels per mode. Panels contain all 20 guides. The measured panels were run in a balanced alternating order: baseline, Columba, Columba, baseline, baseline, Columba.

```text
mode      measured wall seconds                 median wall  median peak RSS KiB  total PAF records  deterministic output
baseline  670.746058,664.124025,663.699258      664.124025   124276               310,310,310        true
columba   9.244436,8.691813,9.099324            9.099324     290436               3095,3095,3095     true
```

Derived ratios:

```text
speedup = baseline median / Columba median = 72.986084
memory_ratio = Columba median peak RSS / baseline median peak RSS = 2.337024
```

Interpretation limits:

- Output volume differs by design: Columba reports an exhaustive verified candidate set, while baseline reports one WFA2-selected hit per sliding window/strand.
- Candidate counts are not exposed by the current CRISPRapido CLI and are recorded as `not_exposed`.
- The benchmark includes output writing time consistently for both modes.
- Filesystem cache was not cleared between runs.
- The `/proc` RSS sampler can miss very short-lived child-process memory peaks.

Recommendations:

- Larger chr22 guide panels are appropriate next because this pilot passed locus-recall eligibility and deterministic-output checks.
- Multithreaded testing is appropriate after adding candidate-count instrumentation or log preservation for automatic Columba mode.
- chr1 benchmarking is a reasonable next scale-up before whole-genome runs.
- Whole-genome benchmarking should wait until runtime, memory, candidate volume, and output-size behavior are characterized on larger single-chromosome panels.
