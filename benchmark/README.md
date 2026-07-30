# CRISPRapido Benchmarks

This directory contains reproducible benchmarks for the Columba-enabled CRISPRapido implementation.

The `controlled/` benchmark checks correctness on a small reference with known Columba SAM outputs. It compares manual Columba SAM import against automatic Columba execution, then verifies that both modes produce byte-identical PAF output after the shared WFA2 verification, CFD scoring, and reporting pipeline.

Future chromosome-scale and genome-scale benchmarks should live alongside this directory and focus on runtime, memory use, and scalability.

The checked-in controlled results were originally introduced with implementation tag `columba-v0.1` and updated after the lowercase FASTA and flanked-window WFA2 verification fixes on `columba-wfa2-cfd`.
