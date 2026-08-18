# Targeted Mismatch Validation

This benchmark isolates substitution handling in the baseline sliding-window
workflow and the automatic Columba candidate-generation workflow. It does not
measure performance and does not mix mismatches with DNA/RNA bulges.

## Existing coverage

The repository's controlled benchmark already includes `m=1,b=0,z=0` and
shows byte-identical manual-SAM and automatic Columba PAF output. Its fixture
contains one forward one-substitution locus. A two-substitution locus is
accepted by a mixed `m=2,b=1,z=2` configuration. Rust unit tests also cover a
single WFA2 mismatch and one imported one-substitution candidate.

That coverage does not provide an independent, locus-level oracle for
mismatch-only `m=2` or `m=3`, reverse-strand mismatch loci, or negative controls
immediately beyond each mismatch threshold.

## Fixture

The deterministic fixture uses one 20-nt guide:

```text
GAGTCCGAGCAGAAGAAGAA
```

Each locus is placed on a separate short reference record with an adjacent
oriented `GG` PAM. The fixture contains:

- two exact loci, one per strand;
- three one-mismatch loci at guide positions 3, 10, and 17;
- three two-mismatch loci with substitutions at distal, internal, proximal,
  and reverse-strand positions;
- two three-mismatch loci, one per strand;
- two four-mismatch negative controls.

`prepare_fixture.py` chooses deterministic flanks and proves by exhaustive
Hamming-distance enumeration that no other 20-nt substring on either strand is
within three mismatches of the guide. `expected_hits.tsv` records coordinates,
strand, mismatch count, compact canonical CIGAR, PAM, and threshold membership.

## Configurations

| Configuration | Mismatches | Bulges | Bulge size | Columba candidate e | Expected loci |
|---|---:|---:|---:|---:|---:|
| A | 1 | 0 | 0 | 1 | 5 |
| B | 2 | 0 | 0 | 2 | 8 |
| C | 3 | 0 | 0 | 3 | 10 |

The four-mismatch records are negative controls for all configurations. The
three-mismatch records are also negative controls for A and B, and the
two-mismatch records are negative controls for A.

## Independent oracle

`correctness_oracle.py` does not call CRISPRapido alignment or acceptance
helpers. For every PAF row it reconstructs the oriented 20-nt reference target
and adjacent PAM directly from `reference.fa`, computes Hamming distance,
constructs a canonical `=`/`X` CIGAR, and checks `nm:i`, `ng:i`, `bs:i`, and
`cg:Z`. It rejects non-20-nt spans, gaps, non-`GG` PAMs, unexpected loci,
duplicates, and rows above the configured mismatch limit.

For reverse-strand substitutions, baseline scanning may emit compact CIGAR
operations in forward-reference order while anchored imported verification
emits them in oriented-guide order. The oracle accepts only these two exact
operation-order representations after independently confirming the same
target, coordinates, strand, mismatch count, and PAM.

`expected_summary.tsv` records the expected accepted-locus totals before any
workflow is run. Measured summaries are generated under the ignored `raw/`
run directory and are not pre-populated.

Success requires both baseline and automatic Columba mode to recover exactly
the expected locus set with no missing or independently invalid rows. Each mode
is run twice, and its two PAF files must be byte-identical.

## Run locally

Build the current release binary first, then run:

```bash
benchmark/performance/mismatch_validation/run_mismatch_validation.sh
```

To choose an output directory:

```bash
benchmark/performance/mismatch_validation/run_mismatch_validation.sh \
  /tmp/crisprapido-mismatch-validation
```

Environment overrides are supported for `CRISPRAPIDO_BIN`, `COLUMBA_BIN`, and
`COLUMBA_BUILD`. The tiny index is generated inside the ignored run directory.

## Optional real-reference check

A chr22 check is not required to establish mismatch-threshold correctness: the
synthetic fixture exercises the complete baseline and automatic integration
with exact ground truth. After it passes, an optional external-validity check
can use ten evenly spaced guides from the existing deterministic chr22
100-guide panel at `m=1` and `m=2`, with `b=0,z=0`. It should remain a
correctness check, not a performance result, and needs a separate locus oracle
because the true set of chr22 mismatch off-targets is not known in advance.
