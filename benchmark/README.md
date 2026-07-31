# CRISPRapido Benchmarks

This directory contains reproducible benchmarks for the Columba-enabled CRISPRapido implementation.

The `controlled/` benchmark checks correctness on a small reference by comparing manual Columba SAM import against automatic Columba execution. For every verification configuration, the benchmark computes `candidate_e = m + b*z`, generates a fresh manual SAM with that exact Columba edit-distance threshold, and verifies that manual and automatic modes produce byte-identical PAF output after the shared WFA2 verification, CFD scoring, filtering, and reporting pipeline.

This benchmark distinguishes Columba candidate generation from downstream CRISPRapido acceptance. Columba `candidate_e` controls which candidate alignments are generated; CRISPRapido `-m`, `-b`, `-z`, and `-f` control which candidates are accepted after WFA2 verification.

Future chromosome-scale and genome-scale benchmarks should live alongside this directory and focus on runtime, memory use, and scalability.

The checked-in controlled results were originally introduced with implementation tag `columba-v0.1` and updated on `columba-wfa2-cfd` after lowercase FASTA normalization, flanked-window WFA2 verification fixes, full-guide accounting, coordinate-anchored imported verification, and candidate-threshold-aligned benchmark generation.
