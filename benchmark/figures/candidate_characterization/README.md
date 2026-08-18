# Figure 4: Candidate characterization and mismatch recovery

This package combines two validated experiments without treating them as one
dataset:

- Panels A, B, and D characterize final Columba-enabled CRISPRapido PAF
  records from the 20-guide CHM13v2 whole-genome correctness campaign
  (`m=0,b=1,z=2`).
- Panel C reports the separate real-reference chr22 mismatch validation using
  ten unique guides with `b=0,z=0`, `PAM=GG`, and one thread.

The package does not run CRISPRapido, Columba, or a benchmark.

## Files

- `figure4_source_data.tsv`: one row for each of the 4,249 whole-genome final
  Columba PAF records, including reconstructed PAM, independent-oracle
  classification, canonical edit event, shared/Columba-only status, reported
  CIGAR, and CFD.
- `figure4_mismatch_source_data.tsv`: compact summary, mismatch-count, and
  strand rows from the finalized chr22 mismatch oracle outputs.
- `build_figure4_source_data.py`: rebuilds the whole-genome per-locus TSV from
  retained correctness outputs and the CHM13v2 reference.
- `build_figure4_mismatch_source_data.py`: rebuilds the compact mismatch TSV
  from finalized `summary.tsv`, `by_mismatch.tsv`, and `by_strand.tsv` files.
- `make_figure4.py`: validates both source-data files and generates the
  combined figure and all standalone panels.
- `figure4_candidate_characterization.pdf` and `.png`: combined four-panel
  manuscript figure.
- `panels/`: independently sized vector PDF and 320 dpi PNG exports.

## Whole-genome source and classification

Panels A, B, and D use the finalized campaign:

`benchmark/performance/chm13v2_whole_genome_20_guides/raw/chm13_wg20_correctness_20260814T164635Z`

The per-locus builder reads the 20 guide sequences and both baseline and
Columba PAFs, streams the CHM13v2 reference to reconstruct strand-aware targets
and PAMs, and applies the same independent oracle used by the finalized
whole-genome benchmark. Shared loci use ordinary biological-locus equivalence
first and the narrowly validated D/I equivalence second. Canonical event labels
come from the oracle CIGAR rather than the reported WFA2 CIGAR.

The whole-genome experiment used `m=0`; Panels A, B, and D therefore do not
measure mismatch recovery. `I` denotes a guide insertion relative to the
reference. `D` denotes a reference deletion, meaning a reference-consuming gap
in the alignment.

## Mismatch-validation source

Panel C uses the independent-oracle output produced by:

`benchmark/performance/mismatch_validation/chr22_10_guides/run_chr22_mismatch_validation.sh`

The experiment used ten exact-copy-class unique chr22 guides. Baseline and
automatic Columba modes were each run twice at `m=1,b=0,z=0` and
`m=2,b=0,z=0`. The strict oracle requires a 20-base target, the requested `GG`
PAM, exact guide/reference/strand/coordinates, matching mismatch count, no gap
operations, and canonical or exactly reversed operation order on reverse-
strand records. It does not use the indel-equivalence rule.

`figure4_mismatch_source_data.tsv` was generated from finalized oracle files
with these SHA-256 values:

- `summary.tsv`: `83b3a4129abb4349ab7a76b3dec0abd6f3f3fcfec00b09437f2a8b315d3fafb8`
- `by_mismatch.tsv`: `1939f991a1a0e8efc93161e1d3c82a79e62f37160c5a842e142537eefb6a9b50`
- `by_strand.tsv`: `c7bf540a1f1f315a323960c27832a9e5b37218c1c4ef8075118456fa83efa232`

## Panel values

### Panel A: whole-genome edit/event composition

| Group | Exact | 1-nt I | 2-nt I | 1-nt D | 2-nt D | Total |
|---|---:|---:|---:|---:|---:|---:|
| Shared with baseline | 43 | 67 | 15 | 1 | 1 | 127 |
| Columba-only | 0 | 532 | 28 | 9 | 5 | 574 |

### Panel B: whole-genome PAM composition

| PAM | AA | AC | AG | AT | CA | CC | CG | CT | GA | GC | GG | GT | TA | TC | TG | TT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Records | 73 | 23 | 73 | 142 | 114 | 25 | 23 | 43 | 46 | 63 | 701 | 59 | 33 | 14 | 2,754 | 63 |

The requested `GG` PAM accounts for 701 records; 3,548 records have another
PAM dinucleotide. The total is 4,249 final Columba PAF records.

### Panel C: real chr22 mismatch recovery

The main grouped bars use the `m=2,b=0,z=0` experiment:

| Mismatch count | Baseline-valid loci | Recovered by Columba |
|---:|---:|---:|
| 0 | 10 | 10 |
| 1 | 21 | 21 |
| 2 | 1,007 | 1,007 |
| **Total** | **1,038** | **1,038** |

There were zero missing baseline-valid loci, zero additional requested-`GG`
Columba loci, and zero independently invalid Columba loci. Forward recovery was
955/955 and reverse recovery was 83/83. The `m=1` check recovered 31/31 loci
(25/25 forward and 6/6 reverse), also with zero missing or invalid loci. Both
modes were deterministic across repeated runs.

### Panel D: whole-genome post-output classification

Panel D is a post-output classification, not a candidate-generation or
filtering funnel:

- final Columba PAF records: 4,249
- independently alignment-valid records: 4,249
- requested-`GG` and oracle-valid loci: 701
- shared baseline loci: 127
- additional oracle-valid Columba candidate loci: 574
- independently invalid Columba records: 0
- non-`GG` records: 3,548

The 574 additional loci are not experimentally confirmed off-targets and are
not assumed to have equal biological relevance.

## CFD audit and limitation

CFD remains in the whole-genome per-locus TSV for provenance but is not plotted.
All 43 exact shared loci have canonical and reported CIGAR `20=`, PAM `GG`, and
reported CFD `1.0000`. Some terminal indel alignments also receive CFD `1.0000`
because the current CFD input conversion can omit or neutralize terminal gap
differences. CFD `1` for these cases does not mean exact biological identity,
so CFD distributions for indel-containing loci are not used as a primary
comparison. Production CFD code is unchanged by this figure package.

## Validation

`make_figure4.py` checks Panels A, B, and D against
`benchmark/performance/chm13v2_whole_genome_20_guides/correctness_summary.tsv`.
It verifies 4,249 total records, 701 requested-`GG` loci, 127 shared loci, 574
additional oracle-valid loci, 3,548 non-`GG` records, and zero independently
invalid records.

For Panel C, it checks both mismatch configurations, requires `b=0,z=0`, zero
missing/Columba-only/invalid loci, deterministic output, mismatch-count totals
that equal the configuration summary, and forward/reverse totals that equal
the same summary. All plotted values are read from
`figure4_mismatch_source_data.tsv`.

## Rebuild

Generate the compact mismatch source data from a finalized validation run:

```bash
python3 benchmark/figures/candidate_characterization/build_figure4_mismatch_source_data.py \
  --run-root /path/to/finalized-chr22-mismatch-run
```

Then generate the combined figure and all standalone panels:

```bash
python3 benchmark/figures/candidate_characterization/make_figure4.py
```

Rebuilding the whole-genome per-locus TSV is optional and requires the retained
ignored whole-genome correctness outputs plus `../data/real_reference/chm13v2.fa`:

```bash
python3 benchmark/figures/candidate_characterization/build_figure4_source_data.py
```

## Figure legend

**Figure 4. Characterization and correctness validation of loci reported by
Columba-enabled CRISPRapido.** (A) Canonical indel-event counts for
independently valid whole-genome `GG` loci shared with baseline CRISPRapido and
loci reported only by the Columba path. `I` is a guide insertion relative to
the reference, and `D` is a reference-consuming deletion. (B) PAM
dinucleotides reconstructed from the whole-genome reference for all final
Columba PAF records; the requested `GG` PAM is highlighted. (C) Independent
validation of mismatch-containing loci for ten real chr22 guides with
`m=2,b=0,z=0`. Columba recovered all 1,038 baseline-valid loci, including all
1,007 two-mismatch loci; the separate `m=1` check recovered 31/31 loci. (D)
Post-output classification of whole-genome Columba records by alignment
validity, requested PAM, and overlap with baseline-valid loci. This is not a
candidate-generation or filtering funnel. The 574 additional loci are
oracle-valid Columba candidate loci, not experimentally confirmed off-targets.

## Standalone captions

**Panel A. Whole-genome edit/event composition.** Counts of canonical oracle
indel events among valid requested-`GG` loci shared with baseline CRISPRapido
and loci reported only by the Columba path. The experiment used `m=0`; `I`
denotes a guide insertion and `D` a reference-consuming deletion.

**Panel B. Whole-genome PAM composition.** Reference-reconstructed PAM
dinucleotides for all 4,249 final Columba PAF records. The requested `GG` PAM
is highlighted; 701 records have `GG` and 3,548 have another PAM.

**Panel C. Real-reference mismatch validation.** Baseline-valid and
Columba-recovered loci by mismatch count for ten unique chr22 guides with
`m=2,b=0,z=0`. Columba recovered all 1,038 baseline-valid loci, including all
1,007 loci with two mismatches; no baseline locus was missing and no Columba
locus was independently invalid.

**Panel D. Whole-genome post-output classification.** Classification of 4,249
final Columba PAF records: all were independently alignment-valid, 701 had the
requested `GG` PAM, 127 were shared with baseline-valid loci, and 574 were
additional oracle-valid Columba candidate loci. This is not a
candidate-generation or filtering funnel.
