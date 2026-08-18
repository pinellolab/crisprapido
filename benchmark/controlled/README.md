# Controlled Correctness Benchmark

This benchmark validates that CRISPRapido's manual Columba-SAM import mode and automatic Columba execution mode produce byte-identical PAF output when they use the same Columba candidate-generation threshold.

For each configuration, the benchmark calculates:

```text
candidate_e = max_mismatches + max_bulges * max_bulge_size
```

The benchmark then runs Columba directly with `-e candidate_e` to generate a fresh `manual.sam`, runs CRISPRapido manual mode with that SAM, and runs CRISPRapido automatic mode with the same downstream verification parameters. The archived `../results/controlled_k*.sam` files are not used for integration equivalence because their `k` value is not necessarily equal to `candidate_e` when bulges are allowed.

The two thresholds have different roles:

- Columba `candidate_e` is a candidate-generation edit-distance bound. It should be a superset large enough to include alignments that downstream CRISPRapido may accept.
- CRISPRapido `-m`, `-b`, `-z`, and `-f` are downstream WFA2 verification limits: maximum mismatches, maximum gap groups, maximum gap-group size, and minimum match fraction.

WFA2 verification remains authoritative. Columba may return candidates that satisfy total edit distance but are rejected because they violate the separate mismatch, gap-group, gap-size, PAM, or match-fraction checks.

Both modes use the same controlled reference, guide, PAM, WFA2 verification path, CFD scoring path, hit filtering, and deterministic PAF reporting. Success means the manual and automatic PAF files are byte-for-byte identical without sorting.

Observed results on `columba-wfa2-cfd` after full-guide accounting, lowercase FASTA normalization, flanked-window ends-free alignment, and coordinate-anchored imported verification:

```text
config	m	b	z	f	candidate_e	manual_candidates	automatic_candidates	manual_records	automatic_records	manual_exit	automatic_exit	paf_byte_identical	stderr_diff
A	0	0	0	0.75	0	12	not_observed	12	12	0	0	yes	none
B	1	0	0	0.75	1	17	not_observed	13	13	0	0	yes	none
C	0	1	1	0.75	1	17	not_observed	16	16	0	0	yes	none
D	0	1	2	0.75	2	22	not_observed	18	18	0	0	yes	none
E	1	1	2	0.75	3	22	not_observed	21	21	0	0	yes	none
F	2	1	2	0.75	4	22	not_observed	22	22	0	0	yes	none
```

`automatic_candidates` is recorded as `not_observed` because the current CRISPRapido CLI deletes the generated temporary SAM by default and does not expose a candidate counter. The manual candidate count is taken from the freshly generated `manual.sam` for the same `candidate_e`, so it represents the candidate set that both modes should use.

To regenerate this benchmark:

```bash
benchmark/scripts/run_controlled_benchmark.sh benchmark/controlled
```

To validate without overwriting checked-in outputs:

```bash
benchmark/scripts/run_controlled_benchmark.sh /tmp/crisprapido-controlled-benchmark-check
```
