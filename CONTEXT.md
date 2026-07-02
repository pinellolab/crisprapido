# crisprapido — Refactor Context for opencode

## What this tool does
Scans a reference genome (FASTA, optionally gzipped) for CRISPR off-target sites.
- Input: a guide RNA sequence (20bp) + PAM (NGG), given as 22bp FASTA entries
- It searches both strands of every chromosome/contig
- Outputs: hit positions, strand, CIGAR, CFD score, mismatch count

## Current search backend: sassy
- Crate: `sassy` (from github.com/RagnarGrootKoerkamp/sassy)
- Used via a thread-local `SEARCHER` (see src/main.rs ~line 623) to avoid re-allocation per call
- Called at ~line 713: `SEARCHER.with(|searcher| searcher.borrow_mut().search(guide, &contig, max_errors))`
- The genome is read record-by-record from a fasta::Reader at runtime (line 975)
- Each contig is searched as it is streamed — genome is NOT pre-loaded into memory

## Goal 1: Replace sassy with ropebwt3
- ropebwt3 is already compiled at: ./ropebwt3/ropebwt3 (binary)
- Source is at: ./ropebwt3/
- Integration approach: call ropebwt3 as a subprocess OR write a Rust FFI wrapper against its C library
- ropebwt3 uses an FM-index (BWT-based) for fast approximate string matching
- The index must be built from the genome FASTA first, then queried per guide RNA

## Goal 2: Load genome index once
- Currently the genome FASTA is streamed per run (acceptable), but if ropebwt3 needs
  to build/load an index, that index must be built ONCE at startup and reused for all
  guide RNA queries — not rebuilt per query
- The index file can be saved to disk (e.g. genome.fa.rbwt) and reloaded on next run

## Goal 3: Configurable mismatches and indels
- Add CLI flag: --mismatches <0|1|2|3>  (default: 3)
  - 0 = perfect match only
  - 1, 2, 3 = allow up to N mismatches
- Add CLI flag: --indels <bool> or --allow-indels (default: false)
  - Controls whether insertions/deletions are permitted in hits
- SNPs: the tool already handles mismatches; make sure SNP-like single-base differences
  are counted within the --mismatches budget

## Existing CLI args to preserve
- --reference / -r : path to genome FASTA (or .fa.gz)
- --threads / -t : number of threads (default: num logical CPUs)
- Keep all existing output columns and format

## Test files available
- test_genome.fa : small synthetic genome (5 lines, 1 chromosome)
- test.fa : guide RNA inputs including Ns
- test_20bp_mismatches.fa : 5 guide RNAs with 0, 1, 1, 1, and 2 mismatches vs genome
- simple_test.fa : basic sequences including NNN ambiguous bases

## Key source files
- src/main.rs : CLI, threading, genome reading, search dispatch
- src/cfd_score.rs : CFD scoring (do not change)
- src/lib.rs : library interface
- Cargo.toml : dependencies (sassy will be replaced or kept as fallback)

## Language & constraints
- Rust (edition 2021)
- ropebwt3 is written in C — use FFI (bindgen) or subprocess (std::process::Command)
- Subprocess is simpler and acceptable if performance is not degraded
- Do not break existing output format or CFD scoring
