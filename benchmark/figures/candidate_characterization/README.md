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
  correctness summary and generates the combined figure and all standalone
  panels with one command.
- `figure4_candidate_characterization.pdf` and `.png`: combined four-panel
  manuscript figure.
- `panels/`: independently sized vector PDF and 320 dpi PNG exports for each
  panel.

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
a mismatch distribution.

`I` denotes a guide insertion relative to the reference. `D` denotes a
reference deletion, meaning a reference-consuming gap in the alignment.

## Panel values

Panel A uses oracle-valid requested-`GG` loci and plots counts:

| Group | Exact | 1-nt I | 2-nt I | 1-nt D | 2-nt D | Total |
|---|---:|---:|---:|---:|---:|---:|
| Shared with baseline | 43 | 67 | 15 | 1 | 1 | 127 |
| Columba-only | 0 | 532 | 28 | 9 | 5 | 574 |

Panel B includes all 4,249 final Columba PAF records:

| PAM | AA | AC | AG | AT | CA | CC | CG | CT | GA | GC | GG | GT | TA | TC | TG | TT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Records | 73 | 23 | 73 | 142 | 114 | 25 | 23 | 43 | 46 | 63 | 701 | 59 | 33 | 14 | 2,754 | 63 |

The requested `GG` PAM accounts for 701 records; 3,548 records have another
PAM dinucleotide.

Panel C uses canonical oracle edit distance:

| Group | ED0 | ED1 | ED2 | Total |
|---|---:|---:|---:|---:|
| Shared with baseline | 43 | 68 | 16 | 127 |
| Columba-only | 0 | 541 | 33 | 574 |

Panel D is a post-output classification, not a candidate-generation or
filtering funnel:

- final Columba PAF records: 4,249
- independently alignment-valid records: 4,249
- requested-`GG` and oracle-valid loci: 701
- shared baseline loci: 127
- additional oracle-valid Columba candidate loci: 574
- independently invalid Columba records: 0
- non-`GG` records: 3,548

The 574 additional loci are not described as experimentally confirmed
off-targets and are not assumed to have equal biological relevance.

## CFD audit and limitation

CFD is retained in the per-locus TSV for provenance but is not plotted in the
main figure. All 43 exact shared loci have canonical and reported CIGAR `20=`,
PAM `GG`, and reported CFD `1.0000`. However, some terminal indel alignments
also receive CFD `1.0000`. The current CFD input conversion can omit a terminal
reference-consuming deletion or represent a terminal guide insertion as an
unpenalized gap-to-gap comparison. Therefore CFD `1` for these cases does not
mean exact biological identity, and CFD scores for indel-containing loci are
not used as a primary biological comparison in Figure 4. Production CFD code
is unchanged by this figure package.

## Validation

The plotting script derives every plotted count from `figure4_source_data.tsv`
and checks final totals against:

`benchmark/performance/chm13v2_whole_genome_20_guides/correctness_summary.tsv`

It also verifies that Panel A and Panel C each contain 127 shared and 574
Columba-only loci, that `GG + non-GG = 4,249`, that `127 + 574 = 701`, that all
canonical mismatch counts are zero, and that only edit distances 0, 1, and 2
occur in the plotted oracle-valid loci.

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

## Figure legend

**Figure 4. Characterization of loci reported by Columba-enabled
CRISPRapido.** (A) Canonical edit-event counts for independently valid `GG`
loci shared with baseline CRISPRapido and loci reported only by the Columba
path. `I` is a guide insertion relative to the reference, and `D` is a
reference-consuming deletion. (B) PAM dinucleotides reconstructed from the
reference for all final Columba PAF records; the requested `GG` PAM is
highlighted. (C) Canonical edit-distance counts for shared and Columba-only
oracle-valid `GG` loci. The benchmark used `m=0`, so edit distance is caused by
indel events rather than mismatches. (D) Post-output classification of final
Columba records by independent alignment validity, requested PAM, and overlap
with baseline-valid loci. This panel is not a candidate-generation or
filtering funnel. The 574 additional loci are oracle-valid Columba candidate
loci, not experimentally confirmed off-targets.

## Standalone captions

**Panel A. Edit/event composition.** Counts of canonical oracle edit events
among valid requested-`GG` loci shared with baseline CRISPRapido and loci
reported only by the Columba path. The benchmark used `m=0`; `I` denotes a
guide insertion and `D` a reference-consuming deletion.

**Panel B. PAM composition.** Reference-reconstructed PAM dinucleotides for all
4,249 final Columba PAF records. The requested `GG` PAM is highlighted; 701
records have `GG` and 3,548 have another PAM.

**Panel C. Edit-distance composition.** Canonical edit-distance counts among
shared and Columba-only oracle-valid `GG` loci. Because the benchmark used
`m=0`, edit distances 1 and 2 reflect insertion/deletion events, not
mismatches.

**Panel D. Post-output classification.** Classification of 4,249 final Columba
PAF records after output: all were independently alignment-valid, 701 had the
requested `GG` PAM, 127 were shared with baseline-valid loci, and 574 were
additional oracle-valid Columba candidate loci. This is not a
candidate-generation or filtering funnel.
