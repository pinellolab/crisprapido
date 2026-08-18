# chr22 10-guide mismatch validation

This small correctness validation compares baseline CRISPRapido with automatic
Columba candidate generation on ten existing exact-copy-class unique guides
from `benchmark/performance/chr22_100_guides/chr22_guides.tsv`.

The mismatch-only configurations are `m=1,b=0,z=0,f=0.75` and
`m=2,b=0,z=0,f=0.75`. Both use `PAM=GG` and one thread. Automatic mode uses
Columba candidate edit-distance bounds `e=1` and `e=2`, respectively. Each
workflow runs twice per guide and configuration, and raw PAF output must be
byte-identical between replicates.

## Independent oracle

`chr22_mismatch_oracle.py` reconstructs each reported target and oriented PAM
directly from the chr22 FASTA. Requested-PAM valid loci must span exactly 20
bases, have Hamming distance no greater than `m`, report matching `nm:i`, have
`ng:i:0` and `bs:i:0`, contain only `=`/`X` CIGAR operations, and have an
adjacent oriented `GG` PAM.

Loci are compared by exact guide, reference, strand, start, and end. The
indel-equivalence rule used by bulge benchmarks is deliberately disabled.
For reverse-strand substitution alignments, only the canonical CIGAR operation
order or its exact reverse is accepted after sequence, coordinate, mismatch,
and PAM reconstruction succeeds.

Non-`GG` rows are excluded from valid loci and reported separately from
alignment-invalid rows. The oracle independently validates reported loci but
does not exhaustively enumerate every possible 20-mer in chr22.

## Run

Build the intended release binary, then run:

```bash
benchmark/performance/mismatch_validation/chr22_10_guides/run_chr22_mismatch_validation.sh \
  /tmp/crisprapido-chr22-mismatch-validation
```

Raw outputs remain outside Git or under the ignored `raw/` directory. This is
a correctness validation and must not be interpreted as a performance result.
