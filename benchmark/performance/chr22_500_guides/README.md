# chr22_500_guides CRISPRapido + Columba Benchmark

This package prepares a deterministic guide panel and Slurm-compatible batched
benchmark for CRISPRapido baseline sliding-window candidate generation versus
automatic Columba candidate generation with WFA2 verification and CFD reporting.

Reference: `/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/data/real_reference/chm13v2_chr22.fa`  
Columba index prefix: `/moosefs/raid5/salehi/crispr2/crispr-progress/columba-crispr-benchmark/results/chm13_chr22_index/chm13v2_chr22`  
Chromosome/header: `22`  
Reference length: `51324926` bp

Biological parameters for correctness and timing:

- PAM: `GG`
- max mismatches: `0`
- max bulges: `1`
- max bulge size: `2`
- minimum match fraction: `0.75`
- threads: `1`
- automatic Columba candidate bound: `candidate_e = m + b*z = 2`

Guide selection is deterministic: the reference is divided into
500 intervals, and one unique forward-strand 20 nt guide with
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

Batch layout: 500 guides; 25 guides per batch; 20 batches.

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
completion markers, and batch ordering remains deterministic.

Biological-locus correctness is evaluated with `correctness_oracle.py`, not by
raw PAF identity. Raw PAF equality is inappropriate for this benchmark because
baseline sliding-window mode and Columba-anchored mode can report different
valid representations of the same local off-target event. In repetitive
sequence, a reference-deletion representation such as `18=1D2=` and a shorter
guide-insertion representation such as `2I18=` can describe the same
PAM-relative neighborhood.

The ordinary comparison first requires same guide, reference, strand, PAM, and
overlapping intervals whose start/end differ by at most the configured maximum
bulge size. If that does not match, a conservative indel-equivalence rule may
merge a pair only when all of the following hold:

- same guide
- same reference
- same strand
- same observed requested PAM
- both alignments are independently valid
- intervals overlap
- same start coordinate
- one representation contains `D` and the other contains `I`
- zero mismatches
- no more than one gap group in either row
- max gap in either row is no larger than configured `z`
- endpoint difference is no larger than the sum of the two gap sizes
- oriented target strings are prefix/suffix-compatible in the local sequence
  context

This rule is intentionally not a generic overlap rule. In the final correctness
campaign, it changed exactly 14 apparent missing baseline rows into recovered
biological loci. Negative controls for same-start different-PAM rows,
different-strand overlaps, overly large endpoint differences, and unrelated
repetitive overlaps produced zero false merges. No independently invalid
Columba alignments were found.

Final correctness summary for
`raw/chr22_500_guides_20260810T040146_correctness_final`:

- baseline raw records: `17719`
- Columba raw records: `54518`
- baseline valid loci: `1748`
- Columba valid loci: `19765`
- shared baseline loci: `1748`
- baseline missing from Columba: `0`
- Columba-only valid loci: `18017`
- invalid Columba records: `0`
- Columba non-`GG` PAM records: `34753`
- timing eligible: `yes`

## Final Timing Result

Timing was run from source commit `ab4a3a4` (`Recover Columba candidates across equivalent indel representations`) using the release binary from that implementation. The timing campaign RUN_ID was `chr22_500_guides_20260810T163619_final_timing`.

The RUN_ID is restart-safe and contains a mixed Slurm provenance: batches `001`-`019` use earlier successful task outputs under the same RUN_ID, while batch `020` records the final array IDs supplied for this campaign (`2218825`, `2218827`, `2218828`, and `2218829`). All recorded batches used commit `ab4a3a4`, branch `columba-wfa2-cfd`, and the same biological parameters.

The measured jobs used the same CHM13v2 chr22 reference, Columba index, 500-guide panel, and biological parameters as the correctness campaign: `PAM=GG`, `m=0`, `b=1`, `z=2`, `f=0.75`, and `threads=1`.

Timing used deterministic Slurm batches of 25 guides each. Per-batch wall times were summed to estimate single-thread panel compute time. This is distinct from concurrent Slurm campaign elapsed time when array tasks run in parallel. The Columba index build time is excluded; output writing time is included for both modes. Peak RSS was measured by the `/proc` sampling wrapper and can miss very brief peaks.

The baseline sliding-window mode has one complete measured full-panel replicate, so no baseline median is reported for this 500-guide benchmark. Automatic Columba mode has three measured full-panel replicates. The observed speedup below is calculated as baseline single measured wall time divided by the median of the three Columba measured wall times.

Final timing summary:

- baseline measured wall time: `19918.452960` s
- baseline user CPU time: `19797.981050` s
- baseline system CPU time: `50.141836` s
- baseline peak RSS: `127608` KiB
- baseline PAF records: `17719`
- Columba measured wall times: `697.473811`, `659.807203`, `653.744126` s
- Columba median wall time: `659.807203` s
- Columba wall-time min/max: `653.744126` / `697.473811` s
- Columba peak RSS values: `308272`, `301824`, `312796` KiB
- Columba median peak RSS: `308272` KiB
- Columba PAF records: `54518`
- Columba output deterministic across measured runs: `true`
- observed speedup using one baseline replicate and the median of three Columba replicates: `30.188293x`
- memory ratio, Columba median peak RSS over baseline measured peak RSS: `2.415773x`

Limitations: this is a 500-guide, chr22-only, single-thread benchmark. The filesystem cache was not cleared. Candidate and fallback counts are not exposed by the CLI. Columba output is more exhaustive than baseline sliding-window output, so raw PAF equality is not a correctness criterion. Non-`GG` PAM reporting remains a separate issue tracked by the correctness summaries.
