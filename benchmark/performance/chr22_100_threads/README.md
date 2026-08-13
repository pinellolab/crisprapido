# chr22_100 Thread-Scaling Benchmark

This package contains the timing-only thread-scaling benchmark for the finalized CRISPRapido + Columba implementation. It reuses the validated `chr22_100_guides` panel and does not regenerate guides or repeat correctness unless a timing run reveals a provenance or output-invariance problem.

Reference: `../data/real_reference/chm13v2_chr22.fa`  
Columba index prefix: `../results/chm13_chr22_index/chm13v2_chr22`  
Guide panel: `benchmark/performance/chr22_100_guides/chr22_guides.tsv`  
Batch layout: 100 guides, 10 batches, 10 guides per batch

Biological parameters are fixed to match the tux05 chr22_100 and chr2_100 benchmarks:

- `PAM=GG`
- `m=0`
- `b=1`
- `z=2`
- `f=0.75`
- Columba candidate bound from CRISPRapido: `candidate_e = m + b*z = 2`

## Threading Behavior Audit

CRISPRapido initializes a global Rayon pool when `-t/--threads` is supplied. In baseline sliding-window mode, windows are processed through `into_par_iter()`, so the reference-window scan is parallelized by Rayon. WFA2 aligners are created per Rayon worker through `map_init`, so WFA2 verification of baseline windows is parallel across windows.

Automatic Columba mode passes the same CLI thread value to Columba as `-t`. Columba help documents this as the number of threads used by Columba candidate generation. After Columba returns SAM candidates, CRISPRapido imported-candidate verification currently iterates candidates serially in Rust. Anchored WFA2 verification, fallback interval generation, fallback WFA2 calls, hit deduplication, and final reporting/filtering are therefore not parallelized by Rayon in imported mode.

Expected scaling implication: baseline can benefit from CRISPRapido threading across scan windows. Automatic Columba can benefit from Columba's internal candidate-generation threading, but imported verification and fallback may become serial bottlenecks as candidate volume grows.

## Design

Run on Slurm partition `tux`, node `tux05`, with thread counts:

```text
1 2 4 8
```

For each thread count:

- baseline: one measured full-panel replicate
- Columba: three measured full-panel replicates

Each Slurm command must set both `-c THREADS` and `THREADS=THREADS`; `run_batch.sh` refuses to run if `THREADS` differs from `SLURM_CPUS_PER_TASK`. CRISPRapido receives `-t THREADS`, and automatic Columba mode forwards that value to Columba as `-t THREADS`.

Raw outputs are written under ignored `raw/<RUN_ID>/threads_<N>/...` directories. Slurm logs are ignored under `slurm_logs/`.

## Existing 1-Thread Baseline for Cost Estimates

The matched-node chr22_100 timing on tux05 reported:

```text
baseline measured wall: 1932.584993 s
Columba median wall: 84.294101 s
baseline peak RSS: 127584 KiB
Columba median RSS: 292428 KiB
```

A full four-thread-count run with one baseline replicate at each thread count and three Columba replicates at each thread count is expected to require roughly 4 baseline panel-equivalents plus 12 Columba panel-equivalents. If baseline threading scales well, the total summed compute should be below four times the 1-thread baseline; if it does not, it may approach that cost. Three baseline replicates per thread count would add substantial compute and should be considered only after the one-replicate scaling shape is known.

## Completed Run

The completed tux05 timing campaign used run prefix:

```text
chr22_100_threads_tux05_20260813T182822Z
```

All thread counts completed the planned timing design:

- baseline: one measured full-panel replicate
- Columba: three measured full-panel replicates
- 10/10 batches per complete panel
- 100/100 guides represented exactly once per complete panel
- zero command failures

The stdout and stderr checksums were invariant across thread counts within each mode, confirming that reported PAF output was deterministic with respect to thread count for this benchmark.

Final wide summary:

```text
threads  baseline_wall_seconds  columba_wall_seconds  baseline_speedup_vs_1thread  columba_speedup_vs_1thread  baseline_parallel_efficiency  columba_parallel_efficiency  baseline_peak_rss_kib  columba_peak_rss_kib  baseline_paf_records  columba_paf_records
1        1947.143865            89.529341             1.000000                      1.000000                    1.000000                        1.000000                      127484                 292548                2892                  12608
2        990.520169             85.237319             1.965779                      1.050354                    0.982890                        0.525177                      131152                 292748                2892                  12608
4        559.864686             78.495753             3.477883                      1.140563                    0.869471                        0.285141                      137348                 294540                2892                  12608
8        294.387895             80.831806             6.614212                      1.107600                    0.826776                        0.138450                      151384                 294488                2892                  12608
```

Interpretation:

- Baseline sliding-window mode scales strongly through 8 threads, reaching 6.61x speedup over the one-thread run with 82.7% parallel efficiency at 8 threads.
- Automatic Columba mode improves only modestly from one to four threads and is slightly slower at eight threads than at four threads. The best observed Columba time is at four threads, with 1.14x speedup over one thread.
- The weak Columba-side scaling is consistent with the implementation audit: Columba candidate generation receives the thread count, but imported-candidate verification, fallback expansion, anchored WFA2 verification, deduplication, and reporting are currently serial in CRISPRapido.
- Baseline memory rises with threads from 127,484 KiB to 151,384 KiB. Columba memory is nearly flat around 292-295 MiB across thread counts.
- Output volume is invariant across thread counts: 2,892 baseline records and 12,608 Columba records.

## Commands

Print planned commands without submitting:

```bash
DRY_RUN=1 ./benchmark/performance/chr22_100_threads/submit_plan.sh
```

A focused smoke test for threads=2 is:

```bash
RUN_ID=chr22_100_threads_tux05_smoke_t2_baseline THREADS=2 ITERATION=smoke \
  sbatch -p tux -w tux05 -c 2 --array=1-1 \
  benchmark/performance/chr22_100_threads/slurm_thread_baseline.sbatch

RUN_ID=chr22_100_threads_tux05_smoke_t2_columba THREADS=2 ITERATION=smoke \
  sbatch -p tux -w tux05 -c 2 --array=1-1 \
  benchmark/performance/chr22_100_threads/slurm_thread_columba.sbatch
```

After all timing jobs complete, aggregate with:

```bash
python3 benchmark/performance/chr22_100_threads/aggregate_thread_scaling.py \
  --run-prefix <RUN_ID_PREFIX>
```

The aggregation writes `thread_scaling_summary.tsv` with columns:

```text
reference guide_count threads mode wall_seconds seconds_per_guide speedup_vs_1thread parallel_efficiency peak_rss_kib paf_records replicates node
```

## Validation Checks

Before submission:

```bash
bash -n benchmark/performance/chr22_100_threads/*.sh benchmark/performance/chr22_100_threads/*.sbatch
python3 -m py_compile benchmark/performance/chr22_100_threads/*.py
benchmark/performance/chr22_100_threads/test_slurm_paths.sh
```

Correctness/output should remain invariant across thread counts. Aggregation reports output checksums and failures so thread-dependent nondeterminism can be detected.
