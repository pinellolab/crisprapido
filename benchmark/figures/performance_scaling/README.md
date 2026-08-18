# Figure 2: Performance scaling

This package builds the manuscript performance figure from finalized benchmark
TSVs. It does not read raw benchmark output or run CRISPRapido or Columba.

## Outputs

- `figure2_performance_scaling.pdf`: vector manuscript figure.
- `figure2_performance_scaling.png`: 320 dpi raster figure.
- `figure2_source_data.tsv`: one row per plotted mode and experiment, including
  source-file and source-row provenance.
- `make_figure2.py`: deterministic source-data extraction, validation, and
  plotting script. One invocation produces both the combined figure and all
  standalone panels.
- `panels/panelA_guide_count_scaling.pdf` and `.png`: standalone chr22
  guide-count scaling.
- `panels/panelB_reference_size_scaling.pdf` and `.png`: standalone matched-node
  reference-size scaling.
- `panels/panelC_memory_scaling.pdf` and `.png`: standalone matched-node peak
  memory scaling.

## Source artifacts

Panel A reads rows 2-4 of:

`benchmark/performance/chr22_500_guides/scaling_summary.tsv`

These are the finalized 20-, 100-, and 500-guide chr22 measurements. The
script cross-checks them against the corresponding `timing_summary.tsv` files.
The 20- and 100-guide runs were measured on `octopus01`; the 500-guide result
was assembled from Slurm batches on `octopus07` through `octopus11`. Panel A
therefore summarizes the established guide-count series but is not a
strictly single-node scaling experiment.

Panels B and C read rows 2-4 of:

`benchmark/performance/chm13v2_whole_genome_20_guides/scaling_summary.tsv`

All three rows use `tux05`. The script cross-checks the chr22 row against
`chr22_100_guides/timing_tux05_summary.tsv`, the chr2 row against
`chr2_100_guides/timing_summary.tsv`, and the whole-genome row against
`chm13v2_whole_genome_20_guides/timing_summary.tsv`.

The whole-genome value is the controlled result: one measured baseline
replicate versus the median of two sequential warm-cache Columba replicates.
The previous contention-affected whole-genome timing is excluded.

## Figure semantics

- **A, guide-count scaling:** total wall time for 20, 100, and 500 guides on
  chr22. Both axes use logarithmic spacing. The 20- and 100-guide values are
  medians of three runs per mode. At 500 guides, baseline has one measured
  panel and Columba is the median of three panels.
- **B, reference-size scaling:** seconds per guide against reference length
  for chr22, chr2, and CHM13v2 on matched `tux05` hardware. Both axes are
  logarithmic.
- **C, memory scaling:** peak aggregate RSS for the same matched-node
  experiments, converted from KiB to GiB using 1 GiB = 1,048,576 KiB.

Numbers between paired series are the observed baseline/Columba wall-time
ratios in panels A and B and Columba/baseline memory ratios in panel C. No
confidence intervals are shown because several baseline measurements have a
single replicate. Values are reported as observations from these benchmark
runs, not general performance guarantees.

## Rebuild

The script requires Python 3.10 or newer and Matplotlib. From the repository
root:

```bash
python3 benchmark/figures/performance_scaling/make_figure2.py
```

The generated source-data TSV records every plotted value, its replicate
interpretation, and its originating file, row, and fields. The script exits
with an error if a consolidated scaling value disagrees with its finalized
timing summary.

## Figure legend draft

**Figure 2. Columba accelerates CRISPRapido candidate generation across guide
counts and reference sizes.** (A) Total wall time for the original
sliding-window workflow (baseline) and automatic Columba candidate generation
followed by anchored WFA2 verification on CHM13v2 chr22 for 20, 100, and 500
guides. (B) Wall time per guide as reference size increases from chr22 to chr2
and the complete CHM13v2 assembly on matched `tux05` hardware. (C) Peak
aggregate resident memory for the matched-node experiments. The controlled
whole-genome result compares one baseline replicate with the median of two
sequential warm-cache Columba replicates. Ratio labels show observed speedup
in A-B and Columba/baseline memory ratio in C. No confidence intervals are
implied; replicate counts are shown beneath the panels.


## Standalone panel captions

**Panel A. Guide-count scaling on chr22.** Total wall time for baseline
CRISPRapido and Columba-enabled CRISPRapido for 20, 100, and 500 guides. Ratio
labels show the observed baseline/Columba speedup. The 20- and 100-guide values
are three-run medians; the 500-guide comparison uses one baseline replicate
and the median of three Columba replicates.

**Panel B. Matched-node reference-size scaling.** Wall time per guide for chr22,
chr2, and the complete CHM13v2 assembly on `tux05`. The controlled whole-genome
comparison uses one baseline replicate and the median of two sequential
warm-cache Columba replicates.

**Panel C. Matched-node peak memory scaling.** Peak aggregate resident memory
for the same `tux05` experiments. Ratio labels show Columba/baseline memory
use. The whole-genome Columba value is the median of two sequential warm-cache
replicates.
