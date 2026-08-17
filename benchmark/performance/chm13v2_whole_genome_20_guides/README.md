# CHM13v2 whole-genome 20-guide pilot

This package prepares a correctness-first pilot comparing CRISPRapido's
original sliding-window candidate generation with automatic Columba candidate
generation followed by the same anchored WFA2 verification, CFD scoring, and
PAF reporting.

## Reference and index

The completed pilot uses the full CHM13v2 reference:

- reference: `../data/real_reference/chm13v2.fa`
- reference records: `25`
- reference length: `3,117,292,070` bp
- reference SHA-256: `15a4ba1246f6021a89699bf5083da7f2bad3f79c86acd7bc1eb0ca3a13164e85`
- index prefix: `../results/chm13_whole_genome_index/chm13v2`
- Columba mapper: `../columba/build_Vanilla/columba`
- Columba index builder: `../columba/build_Vanilla/columba_build`

All paths can be overridden with `CHM13_FASTA`, `COLUMBA_INDEX`,
`COLUMBA_BIN`, and `COLUMBA_BUILD`.

## Panel design

`panel_targets.tsv` defines 20 deterministic target slots distributed across
20 chromosomes: 12 unique, 6 low-copy (2-10 exact whole-genome occurrences),
and 2 repetitive (>10 exact occurrences). For each slot,
`prepare_guides.py` examines valid forward-strand 20 nt protospacers nearest
the chromosome center in deterministic coordinate order. Candidates must use
only A/C/G/T, have an immediately adjacent `GG` PAM, contain no homopolymer
longer than five bases, and not duplicate another selected guide.

Candidate pools are classified by one direct whole-genome Columba `-e 0`
invocation per pool. The completed panel is then independently revalidated by
`validate_exact.py`: every guide must have at least one `20M`/`NM:i:0` Columba
record and at least one automatic CRISPRapido exact-match PAF record. No
correctness or timing job should run before all 20 rows pass.

The finalized panel contains 20 guides on chromosomes 1-18, 22, and X: 12
unique, 6 low-copy, and 2 repetitive guides. Exact validation passed for all
20 guides.

## Biological parameters

- PAM: `GG`
- max mismatches: `0`
- max bulges: `1`
- max bulge size: `2`
- minimum match fraction: `0.75`
- threads: `1`
- automatic Columba candidate bound: `candidate_e = m + b*z = 2`

Correctness uses the independent full-guide oracle and applies ordinary
biological-locus equivalence first. The validated narrow D/I equivalence rule
from the chr22_500 and chr2_100 benchmarks is considered only when ordinary
equivalence fails. Timing is eligible only when all jobs succeed, every valid
baseline locus is recovered, and no independently invalid Columba record is
found.

## Final correctness result

The completed correctness campaign is stored outside Git under
`raw/chm13_wg20_correctness_20260814T164635Z`. It used 10 deterministic
two-guide batches per mode on `tux05`. All 20 baseline and all 20 Columba
guide-level commands exited successfully.

- baseline raw records: `583`
- Columba raw records: `4,249`
- baseline valid biological loci: `127`
- Columba valid biological loci: `701`
- shared baseline loci: `127`
- baseline valid loci missing from Columba: `0`
- Columba-only valid loci: `574`
- independently invalid Columba records: `0`
- Columba records with independently valid alignments but non-`GG` PAMs: `3,548`
- timing eligible: `yes`

The correctness pass accumulated `27,525.006023` baseline wall seconds and
`11,594.921666` Columba wall seconds. Peak sampled RSS was `588,932` KiB for
baseline and `14,622,992` KiB for Columba. These are diagnostic correctness-run
measurements, not final timing claims.

Columba batches 1-4 each took about 2,867-2,871 seconds, while batches 5-10
took about 17-24 seconds. The early tasks ran concurrently against a cold
whole-genome index on MooseFS and spent nearly all elapsed time in I/O/index
loading: their aggregate user CPU time was only 5-6 seconds per batch. Later
tasks ran after the index was cached. Candidate/output volume has a smaller
secondary effect (for example, batch 8 emitted 2,740 records in 23.65 seconds),
so repetitiveness does not explain the early 48-minute runs. Final timing must
therefore use complete, consistently scheduled panel iterations and report the
cache limitation.

## Previous contention-affected timing result

Timing RUN_ID `chm13_wg20_timing_20260815T201449Z` used Slurm partition
`tux`, node `tux05`, and one thread. It contains one complete measured
baseline panel and three complete measured Columba panels. All 40 batches
completed successfully, and every panel contains each of the 20 guides exactly
once.

- baseline measured wall time: `23,641.369151` seconds
- baseline seconds per guide: `1,182.068458`
- baseline user/system CPU: `23,205.777754` / `28.490806` seconds
- baseline peak RSS: `590,888` KiB
- baseline PAF records: `583`
- Columba measured wall times: `8,737.818451`, `8,736.328265`, and
  `8,735.355034` seconds
- Columba median wall time: `8,736.328265` seconds
- Columba median seconds per guide: `436.816413`
- Columba median user/system CPU: `58.935054` / `184.094035` seconds
- Columba median peak RSS: `14,611,636` KiB
- Columba PAF records: `4,249`
- Columba stdout and stderr deterministic across replicates: `true`
- observed speedup using one baseline replicate and the median of three
  Columba replicates: `2.706099x`
- memory ratio, Columba median over baseline: `24.728267x`

Timing outputs are byte-identical to their corresponding validated correctness
outputs. The three Columba panels were launched concurrently, and the first
four batches in every panel spent about 36 minutes loading the cold 11.8 GB
index from MooseFS while later batches completed in seconds. Consequently,
these replicates are highly repeatable measurements of this specific
cold-index contention schedule, but they are not independent isolated or
steady-state Columba measurements.

On the same `tux05` node, whole-genome baseline seconds per guide are
`12.953592x` chr2_100 and `61.165147x` chr22_100, closely tracking reference
length ratios of `12.844391x` and `60.736416x`. Columba seconds per guide are
`156.807709x` and `518.205204x`, respectively, because cold whole-index loading
dominates this campaign. Columba peak RSS is `12.319184x` chr2_100 and scales
approximately with index/reference size. These comparisons are descriptive:
the whole-genome panel has 20 guides with a different copy-number composition,
whereas the chromosome panels have 100 guides.

## Final controlled timing result

The final timing campaign is stored outside Git under
`raw/chm13_wg20_controlled_timing_20260816T225040Z`. All batches ran on
partition `tux`, node `tux05`, with one thread and array concurrency
`1-10%1`. The baseline has one measured panel. Columba `measured_1` is
the first index-using panel and is labeled `cold_start_candidate`;
`measured_2` and `measured_3` are labeled
`warm_cache_candidate` and ran sequentially through Slurm `afterok`
dependencies.

Baseline task 1 initially encountered a transient MooseFS `EIO`. Its
incomplete batch was replaced by successful retry job `2695881`; successful
batches 2-10 from array `2657851` were retained. The finalizer verifies that
this restart-safe panel contains every guide exactly once with consistent
binary, reference, index, panel, parameter, and output checksums.

- baseline wall time: `32,036.953851` seconds
- baseline seconds per guide: `1,601.847693`
- baseline user/system CPU: `31,817.832266` / `23.471411` seconds
- baseline peak RSS: `589,000` KiB
- baseline PAF records: `583`
- Columba cold-start-candidate wall time: `388.085890` seconds
- Columba cold-start-candidate seconds per guide: `19.404294`
- Columba cold-start-candidate user/system CPU:
  `77.132347` / `165.211901` seconds
- Columba cold-start-candidate peak RSS: `14,623,928` KiB
- Columba warm-cache wall times: `242.092775` and `236.552530` seconds
- Columba warm-cache median wall time: `239.322653` seconds
- Columba warm-cache range: `236.552530`-`242.092775` seconds
- Columba warm-cache median seconds per guide: `11.966133`
- Columba warm-cache median user/system CPU:
  `74.337520` / `154.651235` seconds
- Columba warm-cache median peak RSS: `14,631,642` KiB
- Columba PAF records per panel: `4,249`
- cold-start speedup: `82.551195x`
- steady-state speedup, one baseline panel over two-replicate warm median:
  `133.865113x`
- warm-cache memory ratio over baseline: `24.841497x`

All 20 guides and ten batches are complete in every panel, with zero failures
or duplicates. Columba stdout and stderr checksums are identical across all
three panels, and both modes match the validated correctness outputs. The
carried-forward correctness result is 127 valid baseline loci, all 127
recovered, zero baseline loci missing, and zero independently invalid Columba
alignments.

The first Columba label records run order, not a guaranteed physically cold
cache: Linux and MooseFS caches were not cleared. The steady-state number uses
only `measured_2` and `measured_3`; the all-three Columba median is
`242.092775` seconds and is retained in `timing_summary.tsv`.

The earlier contention-affected campaign reported a Columba median of
`8,736.328265` seconds and a `2.706099x` speedup because three panels
and four batches per panel loaded the 11.8 GB index concurrently. Sequential
batches reduce the warm median by `36.504393x`. The controlled baseline is
`1.355123x` slower than the earlier single baseline panel, reflecting
single-replicate scheduling/runtime variability and the restart. Performance
reporting therefore uses the baseline and Columba measurements from this same
controlled campaign, reports cold and warm Columba behavior separately, and
does not use the previous contention-affected speedup.

Regenerate compact artifacts from the ignored raw run with:

`python3 benchmark/performance/chm13v2_whole_genome_20_guides/finalize_timing.py --run-root benchmark/performance/chm13v2_whole_genome_20_guides/raw/chm13_wg20_controlled_timing_20260816T225040Z`

## Workflow

1. Audit/build the index when preparing the assets:
   `BUILD_INDEX=1 benchmark/performance/chm13v2_whole_genome_20_guides/prepare_reference.sh`
2. Prepare and validate the panel:
   `python3 benchmark/performance/chm13v2_whole_genome_20_guides/prepare_guides.py`
3. Submit correctness arrays using the commands printed by `submit_plan.sh`.
4. Aggregate completed arrays with `aggregate_batches.py`, then run
   `correctness_oracle.py` against the same run root.
5. Submit timing only if the oracle reports `timing_eligible=yes`.
6. Finalize a completed controlled timing run with `finalize_timing.py`.

The full panel is split into ten stable two-guide batches. Successful batches
are restart-safe and are skipped unless `FORCE=1`. Raw data and Slurm logs are
ignored. Index construction is always excluded from query timing.

## Measurement limitations

Peak aggregate RSS is sampled from `/proc`; short-lived peaks may be missed.
The filesystem cache is not cleared. Per-batch wall times are summed as
single-core panel compute time and must not be confused with concurrent Slurm
campaign elapsed time. Imported verification and fallback are largely serial,
so this initial pilot stays at one thread. Summing concurrently executed batch
times can overstate sequential panel cost when the tasks contend for cold
shared index pages, as observed in the previous contention-affected campaign.

