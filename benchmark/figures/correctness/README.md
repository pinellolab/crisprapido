# Figure 1: Columba integration and correctness

This package builds the manuscript workflow/correctness figure from the
current production implementation and finalized correctness summaries. It
does not run CRISPRapido, Columba, a correctness oracle, or any benchmark.

## Outputs

- `figure1_correctness.pdf`: vector manuscript figure.
- `figure1_correctness.png`: 320 dpi raster figure.
- `figure1_source_data.tsv`: workflow provenance and every correctness value
  used by the figure.
- `make_figure1.py`: deterministic source validation, data extraction, and
  plotting script. One invocation produces both the combined figure and all
  standalone panels.
- `panels/panelA_workflow.pdf` and `.png`: standalone integration workflow.
- `panels/panelB_whole_genome_recovery.pdf` and `.png`: standalone
  whole-genome recovery result.
- `panels/panelC_recovery_across_benchmarks.pdf` and `.png`: standalone
  cross-benchmark recovery result.

## Panel A: implementation workflow

The workflow was audited against the current source before drawing:

- CLI inputs and Columba/baseline branch selection:
  `src/main.rs:43-185`.
- Original FASTA streaming, overlapping window generation, and parallel scan:
  `src/main.rs:188-270`.
- Original ends-free WFA2 window verification and mismatch/bulge/match-fraction
  filtering: `src/verification.rs:559-699`.
- Columba candidate-distance bound and process execution:
  `src/columba.rs:120-229`.
- Candidate union, SAM parsing, strand/POS/CIGAR span recovery:
  `src/columba.rs:232-435`.
- Coordinate-anchored WFA2 verification and narrowly bounded fallback:
  `src/verification.rs:314-556` and `src/verification.rs:591-627`.
- Aligned CFD target and adjacent PAM extraction:
  `src/verification.rs:160-242`.
- PAM-aware overlap filtering, CFD invocation, and deterministic PAF-style
  output: `src/reporting.rs:39-285` and `src/cfd_score.rs:147-180`.

The current imported-candidate path does not perform an unconstrained flanked
realignment. It verifies the SAM-implied reference interval end-to-end and
uses bounded alternative intervals when fallback is triggered. The schematic
therefore labels this stage as coordinate/span recovery followed by anchored
WFA2 verification. Columba changes candidate generation; WFA2 and the existing
CRISPR/PAM, CFD, overlap-filtering, and PAF-reporting code remain downstream.

## Panels B-C: finalized correctness results

Panel B uses row 2 of:

`benchmark/performance/chm13v2_whole_genome_20_guides/correctness_summary.tsv`

It displays independently validated loci rather than raw PAF records:

- baseline-valid loci: 127
- Columba-valid loci: 701
- baseline-valid loci recovered: 127/127 (100%)
- baseline-valid loci missing: 0
- additional Columba-only oracle-valid candidate loci: 574
- independently invalid Columba records: 0

The 574 additional loci passed the benchmark's independent alignment oracle.
They demonstrate broader candidate-locus recovery by the Columba path; the
figure does not claim that they are biological off-targets of equal relevance
or priority.

Panel C reads row 2 of each finalized correctness summary:

- `benchmark/performance/chr22/correctness_summary.tsv`
- `benchmark/performance/chr22_100_guides/correctness_summary.tsv`
- `benchmark/performance/chr22_500_guides/correctness_summary.tsv`
- `benchmark/performance/chr2_100_guides/correctness_summary.tsv`
- `benchmark/performance/chm13v2_whole_genome_20_guides/correctness_summary.tsv`

All five oracle-backed experiments recovered every baseline-valid locus and
reported zero independently invalid Columba records. Comparisons use ordinary
biological-locus equivalence first and the narrow validated D/I equivalence
rule second; they do not merge loci using broad overlap alone.

The controlled synthetic suite is not mixed into the oracle-locus bars because
`benchmark/controlled/summary.tsv` reports candidate-matched raw PAF identity,
not independent biological-locus counts. Its supported result is shown
separately: all six configurations exited successfully and produced
byte-identical manual/automatic PAF output.

## Rebuild

The script requires Python 3.10 or newer and Matplotlib. From the repository
root:

```bash
python3 benchmark/figures/correctness/make_figure1.py
```

The script validates source-code anchors and correctness-accounting identities
before writing outputs. It exits if a workflow anchor disappears, a benchmark
summary is empty, locus totals are inconsistent, or a controlled configuration
is not successful and byte-identical.

## Figure legend draft

**Figure 1. Integration of Columba candidate generation preserves validated
CRISPRapido locus recovery.** (A) The original workflow applies ends-free WFA2
verification to overlapping reference windows, whereas the integrated workflow
uses Columba edit-distance candidates, SAM-derived coordinates, and anchored
WFA2 verification with bounded fallback. Both paths retain CRISPRapido's
downstream PAM/CRISPR filtering, CFD scoring, overlap filtering, and PAF-style
reporting. (B) In the 20-guide CHM13v2 whole-genome benchmark, Columba recovered
all 127 baseline-valid loci and 574 additional oracle-valid candidate loci,
with no independently invalid Columba records. (C) All baseline-valid loci were
recovered across five oracle-backed chromosome and whole-genome experiments;
labels report recovered/baseline-valid counts, missing baseline loci, and
invalid Columba records. The controlled synthetic suite separately produced
byte-identical PAF output in all six candidate-matched configurations.


## Standalone panel captions

**Panel A. CRISPRapido integration workflow.** The baseline workflow searches
overlapping reference windows before ends-free WFA2 verification. The
Columba-enabled workflow obtains edit-distance candidates and SAM-derived
coordinates before anchored WFA2 verification with bounded fallback. WFA2,
PAM/CRISPR filtering, CFD scoring, and PAF-style reporting remain downstream.

**Panel B. Whole-genome validated locus recovery.** In the 20-guide CHM13v2
benchmark, Columba recovered all 127 baseline-valid loci and reported 574
additional oracle-valid candidate loci, with no independently invalid Columba
records.

**Panel C. Correctness across benchmark scales.** Columba recovered every
baseline-valid locus in five oracle-backed chromosome and whole-genome
experiments, with no independently invalid Columba records. The controlled
synthetic benchmark separately produced byte-identical PAF output for all six
candidate-matched configurations.
