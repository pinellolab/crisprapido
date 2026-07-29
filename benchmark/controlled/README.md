# Controlled Correctness Benchmark

This benchmark validates that automatic Columba execution produces the same CRISPRapido output as importing an existing Columba SAM file.

For each edit-distance threshold `k`, the benchmark runs:

- Manual mode: `--columba-sam ../results/controlled_k{k}.sam`
- Automatic mode: `--columba-bin ... --columba-index ... -m k`

Both modes use the same controlled reference, guide, PAM, WFA2 verification path, CFD scoring path, hit filtering, and PAF reporting. The value `k` is passed through CRISPRapido's existing `--max-mismatches` option and through automatic Columba execution as the Columba edit-distance threshold.

Success means the manual and automatic PAF files are byte-for-byte identical without sorting.

Observed results for `columba-v0.1`:

```text
k	manual_records	automatic_records	paf_byte_identical	manual_exit	automatic_exit	manual_candidates	automatic_candidates	stderr_diff
0	12	12	yes	0	0	12	not_observed	none
1	16	16	yes	0	0	17	not_observed	qname_or_other
2	21	21	yes	0	0	22	not_observed	qname_or_other
3	21	21	yes	0	0	22	not_observed	qname_or_other
4	21	21	yes	0	0	22	not_observed	qname_or_other
```

The expected stderr difference for `k=1..4` is only the query name in the WFA2-failed candidate warning: archived SAM files use `guide_20bp`, while automatic mode creates a temporary query FASTA with header `guide`.

Both modes also emit the same near-end PAM warning where the adjacent PAM cannot be extracted:

```text
Warning: unable to extract adjacent PAM for near_end:10..30 on strand +
```

Automatic Columba candidate counts are recorded as `not_observed` because the current CRISPRapido CLI deletes the generated temporary SAM by default and does not expose a candidate counter.

To regenerate this benchmark:

```bash
benchmark/scripts/run_controlled_benchmark.sh benchmark/controlled
```

To validate without overwriting checked-in outputs:

```bash
benchmark/scripts/run_controlled_benchmark.sh /tmp/crisprapido-controlled-benchmark-check
```
