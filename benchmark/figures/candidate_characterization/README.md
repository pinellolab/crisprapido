# Figure 4: Characterization of Columba-reported loci

This package characterizes final Columba-enabled CRISPRapido PAF records from
the finalized 20-guide CHM13v2 whole-genome correctness campaign. It does not
run CRISPRapido, Columba, or a benchmark.

## Files

- `figure4_source_data.tsv`: one row for each of the 4,249 final Columba PAF
  records, including reconstructed PAM, independent-oracle classification,
  canonical edit event, shared/Columba-only status, reported CIGAR, and CFD.
- `build_figure4_source_data.py`: reproducibly rebuilds the compact TSV from
  retained correctness PAFs and the CHM13v2 reference.
- `make_figure4.py`: validates the compact TSV against the finalized
  correctness summary and renders the figure.
- `figure4_candidate_characterization.pdf`: vector manuscript figure.
- `figure4_candidate_characterization.png`: 320 dpi raster figure.

## Source and classification

The source campaign is:

`benchmark/performance/chm13v2_whole_genome_20_guides/raw/chm13_wg20_correctness_20260814T164635Z`

The builder reads the 20 guide sequences and both baseline and Columba PAFs,
streams the existing CHM13v2 reference to reconstruct strand-aware targets and
adjacent PAMs, and applies the same independent oracle used by the finalized
whole-genome benchmark. Shared loci use the same deterministic one-to-one
matching order, ordinary locus equivalence, and narrowly validated D/I
equivalence from:

- `benchmark/performance/chm13v2_whole_genome_20_guides/correctness_oracle.py`
- `benchmark/performance/chr22_500_guides/correctness_oracle.py`

Canonical edit classes come from the oracle CIGAR rather than the reported
WFA2 CIGAR. This makes event labels deterministic when indels have equivalent
placements in repetitive sequence. The whole-genome benchmark used `m=0`, so
Figure 4 characterizes exact and insertion/deletion events; it does not measure
a biological mismatch distribution.

Only independently alignment-valid records with the requested `GG` PAM enter
the shared versus Columba-only and CFD comparisons. Non-`GG` rows remain in
the source table and Panel B because they were present in final PAF output.
The 574 Columba-only records are described as additional oracle-valid Columba
candidate loci, not as true biological off-targets of equal relevance.

## Final values

- final Columba PAF records: 4,249
- independently alignment-valid records: 4,249
- requested-`GG` and oracle-valid loci: 701
- shared baseline loci: 127
- additional Columba-only oracle-valid loci: 574
- independently invalid Columba records: 0
- non-`GG` records: 3,548

The plotting script cross-checks these totals against:

`benchmark/performance/chm13v2_whole_genome_20_guides/correctness_summary.tsv`

Panel D begins with final reported PAF rows. It is a post-output
classification, not a candidate-generation or filtering funnel. The original
Columba SAM counts and internal verification/fallback counters were not
retained by this campaign.

## Rebuild

From the CRISPRapido repository root:

```bash
python3 benchmark/figures/candidate_characterization/build_figure4_source_data.py
python3 benchmark/figures/candidate_characterization/make_figure4.py
```

The source-data builder requires the retained ignored raw campaign and the
reference at `../data/real_reference/chm13v2.fa`. Both paths can be overridden:

```bash
python3 benchmark/figures/candidate_characterization/build_figure4_source_data.py \
  --run-root /path/to/chm13_wg20_correctness_20260814T164635Z \
  --reference /path/to/chm13v2.fa
```

## Figure legend draft

**Figure 4. Characterization of candidate loci reported by Columba-enabled
CRISPRapido.** (A) Canonical edit-event composition of independently valid
`GG` loci shared with baseline CRISPRapido and loci reported only by the
Columba path. The benchmark allowed no mismatches (`m=0`). (B) Dinucleotide
PAMs reconstructed from the reference for all final Columba PAF records; `GG`
is highlighted. (C) Reported CFD-score distributions for shared and
Columba-only oracle-valid `GG` loci. (D) Post-output classification of final
Columba records by independent alignment validity, requested PAM, and overlap
with baseline-valid loci. The 574 additional records are oracle-valid Columba
candidate loci and are not assumed to have equal biological relevance.
