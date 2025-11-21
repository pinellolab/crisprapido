# CRISPRapido
![CRISPRapido Logo](crisprapido.png)

CRISPRapido is a reference-free tool for comprehensive detection of CRISPR off-target sites using complete genome assemblies. Unlike traditional approaches that rely on reference genomes and variant files, CRISPRapido directly analyzes haplotype-resolved assemblies to identify potential off-targets arising from any form of genetic variation. By leveraging efficient approximate string matching algorithms and parallel processing, CRISPRapido enables fast scanning of whole genomes while considering both mismatches and DNA/RNA bulges. The tool is particularly valuable for therapeutic applications, where comprehensive off-target analysis is critical for safety assessment. CRISPRapido can process both complete assemblies and raw sequencing data, providing flexibility for different analysis scenarios while maintaining high computational efficiency through its robust Rust implementation.

## Features
- Fast parallel scanning of genomic sequences
- Support for both gzipped and plain FASTA files
- Configurable mismatch and bulge tolerances
- Automatic reverse complement scanning
- PAF-format output compatible with downstream analysis tools
- Multi-threaded processing for improved performance
- CFD (Cutting Frequency Determination) scoring for off-targets

## Installation

### Prerequisites
- Rust toolchain (install from https://rustup.rs/)

### Simple Installation

```bash
# Clone the repository
git clone https://github.com/FarnazSalehi94/crisprapido.git
cd crisprapido

# Build with Cargo
cargo build --release

# The binary will be at:
./target/release/crisprapido
```

### Install System-wide

```bash
# Install from local directory
cargo install --path .

# Or install directly from GitHub
cargo install --git https://github.com/FarnazSalehi94/crisprapido.git
```

## Usage
```bash
./target/release/crisprapido -r <reference.fa> -g <guide_sequence> -p <pam_sequence> [OPTIONS]
```

### Required Arguments
- `-r, --reference <FILE>`: Input reference FASTA file (supports .fa and .fa.gz)
- `-g, --guide <SEQUENCE>`: Guide RNA sequence (without PAM)
- `-p, --pam <SEQUENCE>` : PAM sequence for CFD

### Optional Arguments
- `-m, --max-mismatches <NUM>`: Maximum number of mismatches allowed (default: 4)
- `-b, --max-bulges <NUM>`: Maximum number of bulges allowed (default: 1)
- `-z, --max-bulge-size <NUM>`: Maximum size of each bulge in bp (default: 2)
- `-w, --window-size <NUM>`: Size of sequence window to scan (default: 4x guide length)
- `-t, --threads <NUM>`: Number of threads to use (default: number of logical CPUs)
- `--no-filter`: Disable all filtering (report every alignment)

## Output Format
CRISPRapido outputs results in the Pairwise Alignment Format (PAF), which is widely used for representing genomic alignments. Each line represents a potential off-target site with the following tab-separated fields:

| Column | Field | Description |
|--------|-------|-------------|
| 1 | Query name | "Guide" (the guide RNA sequence) |
| 2 | Query length | Length of the guide RNA |
| 3 | Query start | 0-based start position in the guide sequence |
| 4 | Query end | 0-based end position in the guide sequence |
| 5 | Strand | '+' (forward) or '-' (reverse complement) |
| 6 | Target name | Reference sequence name (e.g., chromosome) |
| 7 | Target length | Length of the target reference sequence |
| 8 | Target start | 0-based start position in reference |
| 9 | Target end | 0-based end position in reference |
| 10 | Matches | Number of matching bases |
| 11 | Block length | Total alignment block length |
| 12 | Mapping quality | Always 255 for CRISPRapido |

Additionally, CRISPRapido includes these custom tags:

| Tag | Description |
|-----|-------------|
| `as:i` | Alignment score (lower is better) |
| `nm:i` | Number of mismatches |
| `ng:i` | Number of gaps (indels) |
| `bs:i` | Biggest gap size in bases |
| `cg:Z` | CIGAR string representing alignment details |
| `cf:f` | CFD score |

### CFD Score
The Cutting Frequency Determination (CFD) score estimates the likelihood of a guide RNA cutting at an off-target site.
The score ranges from 0.0 to 1.0,  taking into account:
- Position-specific mismatch penalties
- PAM sequence efficiency
- Bulge and gap effects

This implementation requires two data files:
- `mismatch_scores.txt` : Position-specific mismatch penalties
- `pam_scores.txt` : Efficiency scores for different PAM sequences

### Example Output
```
Guide   20      0       20      +       chr1    248956422       10050   10070   19      21      255     as:i:6  nm:i:1  ng:i:0  bs:i:0  cg:Z:19=1X  cf:f:0.0549
```
This indicates:
- A 20bp guide RNA aligned to chromosome 1
- Position 10050-10070 on the forward strand
- 19 bases match with 1 mismatch (nm:i:1)
- No gaps (ng:i:0)
- Alignment score of 6 (as:i:6)
- CIGAR string shows 19 matches followed by 1 mismatch

### PAF Format Specification
For more details on the PAF format, see the [official specification](https://github.com/lh3/miniasm/blob/master/PAF.md) from the developers of miniasm.

## Example

```bash
# Basic usage
./target/release/crisprapido -r genome.fa -g ATCGATCGATCG -p GG -m 3 -b 1 -z 2

# Quick test with a small file
echo ">test_seq" > test.fa
echo "AAATCGATCGATCGAAATCG" >> test.fa
./target/release/crisprapido -r test.fa -g ATCGATCGATCG -p GG -m 1
```

## Testing

Run the test suite:
```bash
cargo test
```

Enable debug output during development:
```bash
cargo run --features debug
```

## License
See LICENSE file

## Citation
Stay tuned!
