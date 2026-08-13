# chr2_100_guides CRISPRapido + Columba Benchmark

This package prepares a deterministic guide panel and Slurm-compatible batched
benchmark for CRISPRapido baseline sliding-window candidate generation versus
automatic Columba candidate generation with WFA2 verification and CFD reporting.

Reference: `/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/data/real_reference/chm13v2_chr2.fa`  
Columba index prefix: `/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/results/chm13_chr2_index/chm13v2_chr2`  
Chromosome/header: `2`  
Reference length: `242696752` bp

Biological parameters for correctness and timing:

- PAM: `GG`
- max mismatches: `0`
- max bulges: `1`
- max bulge size: `2`
- minimum match fraction: `0.75`
- threads: `1`
- automatic Columba candidate bound: `candidate_e = m + b*z = 2`

Guide selection is deterministic: the reference is divided into
100 intervals, and one unique forward-strand 20 nt guide with
an immediately adjacent `GG` PAM is selected nearest each interval center. If a nominal interval has no valid guide after filters, the nearest unused valid guide globally is selected and marked with `inside_interval=false`.
Guides with non-ACGT bases or homopolymer runs longer than five bases are
excluded before interval selection.

Raw outputs are excluded by `.gitignore`. Use `prepare_guides.py` to regenerate
the panel and exact-match validation summaries. Use the Slurm scripts for
batched correctness and measured runs, then `aggregate_batches.py` to combine
completed batches.

The Slurm scripts intentionally do not set a partition or account. Provide
site-specific values with `SBATCH_PARTITION`, `SBATCH_ACCOUNT`, and related
environment variables when submitting.

## Batch and Slurm Design

Batch layout: 100 guides; 10 guides per batch; 10 batches.

Each batch has a stable `batch_###` ID in `batches.tsv`, an ordered guide list,
and a SHA-256 checksum over that guide order. `run_batch.sh` writes each batch
to an independent raw directory and creates a `SUCCESS` marker only after all
guides in the batch exit with status 0. Existing successful batches are skipped
unless `FORCE=1` is set.

Slurm array scripts are provided for correctness baseline, correctness Columba,
measured baseline, and measured Columba runs. They record commit, branch,
hostname, date, Slurm task ID, command lines, exit status, checksums, wall time,
CPU time, and peak RSS through `collect_peak_rss.py`.

Recommended resources are in `resource_recommendations.tsv`. No partition or
account is hard-coded; use site-specific `sbatch` options or environment values
when submitting. Exact future submission commands are listed in
`submission_commands.txt`.

Measurement design for this stage:

- full correctness comparison for all guides
- one complete measured baseline pass
- three measured Columba passes
- single thread
- no filesystem cache clearing
- Columba index construction excluded
- no baseline median should be claimed from one baseline replicate

Aggregation is performed with `aggregate_batches.py`. It verifies that every
expected guide appears exactly once for each mode/iteration, all expected batches
are present, no duplicate guide assignments appear, all successful batches have
completion markers, and batch ordering remains deterministic. Biological-locus
correctness aggregation should be run after both correctness modes complete,
using the same independent oracle method as the existing chr22 packages.

## Final Correctness Result

Correctness was aggregated from `raw/chr2_100_guides_correctness_20260812T142104Z` using the same independent full-guide oracle and conservative biological-locus equivalence used for the finalized chr22_500 benchmark. The ordinary biological-locus rule is applied first; the narrow D/I indel-equivalence rule is applied only for validated equivalent indel representations.

Final correctness summary:

- baseline raw records: `848`
- Columba raw records: `6658`
- baseline valid loci: `222`
- Columba valid loci: `5514`
- shared baseline loci: `222`
- baseline missing from Columba: `0`
- Columba-only valid loci: `5292`
- invalid Columba records: `0`
- Columba non-`GG` PAM records: `1144`
- timing eligible: `yes`

## Final Timing Result

Timing was run on Slurm partition `tux` with all measured tasks on node `tux05`. This placement is recorded explicitly because the earlier chr22 timing packages were run on a different node class, so cross-reference runtime comparisons are hardware-confounded.

Timing was run from commit `63da0e8` (`Add 500-guide chr22 scaling benchmark`). Production source is unchanged from implementation commit `ab4a3a4` (`Recover Columba candidates across equivalent indel representations`). The timing RUN_ID was `chr2_100_guides_timing_20260813T133456Z`.

The measured jobs used the same CHM13v2 chr2 reference, Columba index, 100-guide panel, and biological parameters as the correctness campaign: `PAM=GG`, `m=0`, `b=1`, `z=2`, `f=0.75`, and `threads=1`.

Timing used deterministic Slurm batches of 10 guides each. Per-batch wall times were summed to estimate single-thread panel compute time. This is distinct from concurrent Slurm campaign elapsed time when array tasks run in parallel. The Columba index build time is excluded; output writing time is included for both modes. Peak RSS was measured by the `/proc` sampling wrapper and can miss very brief peaks.

The baseline sliding-window mode has one complete measured full-panel replicate, so no baseline median is reported for this 100-guide chr2 benchmark. Automatic Columba mode has three measured full-panel replicates. The observed speedup below is calculated as baseline single measured wall time divided by the median of the three Columba measured wall times.

Final timing summary:

- baseline measured wall time: `9125.410422` s
- baseline user CPU time: `9055.472011` s
- baseline system CPU time: `9.231481` s
- baseline peak RSS: `576084` KiB
- baseline PAF records: `848`
- Columba measured wall times: `278.720803`, `278.568223`, `275.215601` s
- Columba median wall time: `278.568223` s
- Columba wall-time min/max: `275.215601` / `278.720803` s
- Columba peak RSS values: `1186080`, `1186088`, `1186688` KiB
- Columba median peak RSS: `1186088` KiB
- Columba PAF records: `6658`
- Columba output deterministic across measured runs: `true`
- observed speedup using one baseline replicate and the median of three Columba replicates: `32.758261x`
- memory ratio, Columba median peak RSS over baseline measured peak RSS: `2.058880x`

## chr22 Comparison Caveat

The chr2 reference length is `242696752` bp, while chr22 is `51324926` bp, giving a chr2/chr22 reference length ratio of `4.728633x`. However, chr2 timing was run entirely on `tux05`, while the earlier chr22 timing packages were not run on `tux05`. The combined scaling table therefore provides descriptive cross-reference ratios only; it should not be interpreted as pure genome-size scaling.

Compared with chr22_100 measured values, chr2_100 has lower output volume but longer runtime:

- baseline wall ratio chr2/chr22_100: `2.664804x`
- Columba wall ratio chr2/chr22_100: `4.336048x`
- chr22_100 baseline records / valid loci: `2892` / `690`
- chr2_100 baseline records / valid loci: `848` / `222`
- chr22_100 Columba records / valid loci: `12598` / `8896`
- chr2_100 Columba records / valid loci: `6658` / `5514`

Limitations: this is a 100-guide, chr2-only, single-thread benchmark. The filesystem cache was not cleared. Candidate and fallback counts are not exposed by the CLI. Columba output is more exhaustive than baseline sliding-window output, so raw PAF equality is not a correctness criterion. Non-`GG` PAM reporting remains a separate issue tracked by the correctness summaries.

