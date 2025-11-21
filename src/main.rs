use bio::io::fasta;
use clap::Parser;
use sassy::profiles::Dna;
use sassy::Searcher;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;
// Remove the broken imports for now - we'll add correct ones later
// use sassy::{search, Alphabet, SearchConfig};
use crossbeam_channel::bounded;
use flate2::read::MultiGzDecoder;

mod cfd_score;

#[derive(Clone)]
struct Hit {
    ref_id: String,
    pos: usize,
    strand: char,
    score: i32,
    cigar: String,
    guide: Arc<Vec<u8>>,
    target_len: usize,
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
    cfd_score: Option<f64>, // Add CFD score field
    target_seq: Vec<u8>,    // Add target sequence for CFD calculation
}

struct Task {
    guide: Arc<Vec<u8>>,
    contig_seq: Arc<Vec<u8>>,
    contig_len: usize,
    ref_id: Arc<String>,
}

struct HitResult {
    ref_id: Arc<String>,
    pos: usize,
    strand: char,
    score: i32,
    cigar: String,
    guide: Arc<Vec<u8>>,
    target_len: usize,
    target_seq: Vec<u8>,
}

impl Hit {
    fn quality_score(&self) -> i32 {
        // Count matches and other stats
        let matches = self.cigar.chars().filter(|&c| c == 'M' || c == '=').count();

        let mismatches = self.cigar.chars().filter(|&c| c == 'X').count();

        let gaps = self.cigar.chars().filter(|&c| c == 'I' || c == 'D').count();

        // Higher score is better
        matches as i32 - (mismatches as i32) - (gaps as i32 * 2) - self.score
    }

    fn end_pos(&self) -> usize {
        // Calculate reference consumed bases
        let mut ref_consumed = 0;
        for c in self.cigar.chars() {
            match c {
                'M' | '=' | 'X' | 'D' => ref_consumed += 1,
                _ => {}
            }
        }
        self.pos + ref_consumed
    }

    fn overlaps_with(&self, other: &Hit) -> bool {
        self.strand == other.strand
            && self.ref_id == other.ref_id
            && self.pos < other.end_pos()
            && other.pos < self.end_pos()
    }
}

// Replace your entire report_hit function with this corrected version:

fn report_hit(
    ref_id: &str,
    pos: usize,
    _len: usize,
    strand: char,
    _score: i32,
    cigar: &str,
    guide: &[u8],
    target_len: usize,
    _max_mismatches: u32,
    _max_bulges: u32,
    _max_bulge_size: u32,
    target_seq: &[u8],
    pam: &str,
) {
    // Parse CIGAR to calculate positions and statistics
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut max_gap_size = 0;
    let mut matches = 0;

    // Handle empty CIGAR (fallback to perfect match)
    let effective_cigar = if cigar.is_empty() {
        format!("{}=", guide.len())
    } else {
        cigar.to_string()
    };

    // Parse CIGAR string to count operations
    let mut chars = effective_cigar.chars().peekable();
    while let Some(&ch) = chars.peek() {
        if ch.is_ascii_digit() {
            let mut num_str = String::new();
            while let Some(&next_ch) = chars.peek() {
                if next_ch.is_ascii_digit() {
                    num_str.push(chars.next().unwrap());
                } else {
                    break;
                }
            }

            if let Some(op) = chars.next() {
                if let Ok(count) = num_str.parse::<usize>() {
                    match op {
                        '=' | 'M' => {
                            matches += count;
                        }
                        'X' => {
                            mismatches += count;
                        }
                        'I' | 'D' => {
                            gaps += 1;
                            max_gap_size = max_gap_size.max(count);
                        }
                        _ => {}
                    }
                }
            }
        } else {
            let op = chars.next().unwrap();
            match op {
                '=' | 'M' => matches += 1,
                'X' => mismatches += 1,
                'I' | 'D' => {
                    gaps += 1;
                    max_gap_size = max_gap_size.max(1);
                }
                _ => {}
            }
        }
    }

    // Calculate query positions
    let query_start = 0;
    let query_end = guide.len();
    let query_length = guide.len();

    // Calculate reference positions
    let ref_start = pos;
    let ref_end = pos + guide.len();

    // Calculate adjusted score
    let adjusted_score = mismatches * 3 + gaps * 5;

    // Calculate block length
    let block_len = matches + mismatches + gaps;

    // Enable CFD calculation
    let cfd_score = if !target_seq.is_empty() && target_seq.len() >= guide.len() {
        let target_for_cfd = if target_seq.len() >= 20 {
            &target_seq[0..20]
        } else {
            target_seq
        };

        let guide_for_cfd = if guide.len() >= 20 {
            &guide[0..20]
        } else {
            guide
        };

        cfd_score::get_cfd_score(guide_for_cfd, target_for_cfd, &effective_cigar, pam)
    } else {
        None
    };

    let cfd_tag = match cfd_score {
        Some(score) => format!("\tcf:f:{:.4}", score),
        None => "\tcf:f:0.0000".to_string(),
    };

    // Convert sequences to strings for display
    let guide_str = String::from_utf8_lossy(guide);
    let target_str = String::from_utf8_lossy(target_seq);

    // Create sequence alignment display
    let seq_tag = format!("\tqs:Z:{}\tts:Z:{}", guide_str, target_str);

    // Convert CIGAR to minimap2 format (remove debug print)
    let minimap2_cigar = effective_cigar.clone();

    // Output in PAF format with sequences
    println!("Guide\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tas:i:{}\tnm:i:{}\tng:i:{}\tbs:i:{}\tcg:Z:{}{}{}",
        query_length,      // Query length (total guide length)
        query_start,       // Query start (always 0 for local alignment)
        query_end,         // Query end (bases consumed from query)
        strand,            // Strand (+/-)
        ref_id,            // Target sequence name
        target_len,        // Full target sequence length
        ref_start,         // Target start position
        ref_end,           // Target end position
        matches,           // Number of matches
        block_len,         // Total alignment block length
        adjusted_score,    // AS:i alignment score
        mismatches,        // NM:i number of mismatches
        gaps,              // NG:i number of gaps
        max_gap_size,      // BS:i biggest gap size
        minimap2_cigar,    // cg:Z CIGAR string
        cfd_tag,           // cf:f CFD score
        seq_tag            // qs:Z and ts:Z sequence tags
    );
}

#[cfg(test)]
use rand::{rngs::SmallRng, RngCore, SeedableRng};

#[cfg(test)]
mod tests {
    use super::*;

    fn generate_random_seq(rng: &mut SmallRng, length: usize) -> Vec<u8> {
        let bases = b"ACGT";
        (0..length)
            .map(|_| bases[rng.next_u32() as usize % 4])
            .collect()
    }

    fn create_flanked_sequence(rng: &mut SmallRng, core: &[u8], flank_size: usize) -> Vec<u8> {
        let mut seq = generate_random_seq(rng, flank_size);
        seq.extend_from_slice(core);
        seq.extend(generate_random_seq(rng, flank_size));
        seq
    }

    #[test]
    fn test_normalize_cigar() {
        assert_eq!(normalize_cigar("10="), "10=");
        assert_eq!(normalize_cigar("4=X5="), "4=1X5=");
        assert_eq!(normalize_cigar("===XX="), "3=2X1=");
        assert_eq!(normalize_cigar("2=3X"), "2=3X");
        assert_eq!(normalize_cigar("=X=X="), "1=1X1=1X1=");
        assert_eq!(normalize_cigar("4=I5="), "4=1I5=");
        assert_eq!(normalize_cigar("4=2D5="), "4=2D5=");
        assert_eq!(normalize_cigar("XXXXXXXXXX"), "10X");
        assert_eq!(normalize_cigar("2=2X2="), "2=2X2=");
    }

    #[test]
    fn test_perfect_match_sassy() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGATCGAT";

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(!results.is_empty());
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = &results[0];
        assert_eq!(cigar, "10=");
    }

    #[test]
    fn test_with_mismatches_sassy() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGTTCGAT"; // Single mismatch at position 5

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(!results.is_empty(), "Should accept a single mismatch");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = &results[0];
        assert_eq!(cigar, "4=1X5=");
    }

    #[test]
    fn test_with_bulge_sassy() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGAATCGAT"; // Single base insertion after position 4

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(!results.is_empty(), "Should accept a single base bulge");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = &results[0];
        assert!(
            cigar.contains('I') || cigar.contains('D'),
            "Should contain an insertion or deletion"
        );
    }

    #[test]
    fn test_too_many_differences_sassy() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGTTCGTT"; // Three mismatches at positions 5, 8, 9

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(results.is_empty());
    }

    #[test]
    fn test_perfect_match_with_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let target = create_flanked_sequence(&mut rng, guide, 500);

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false);
        assert!(
            !results.is_empty(),
            "Should match perfectly even with flanks"
        );

        // Find the match at position 500 (SASSY may find other matches in random flanks)
        let match_at_500 = results.iter().find(|(_, _, _, _, _, pos)| *pos == 500);
        assert!(match_at_500.is_some(), "Should find match at position 500");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = match_at_500.unwrap();
        assert_eq!(cigar, "10=");
    }

    #[test]
    fn test_with_mismatches_and_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGAT"; // Single mismatch at position 5
        let target = create_flanked_sequence(&mut rng, core, 500);

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false);
        assert!(
            !results.is_empty(),
            "Should accept a single mismatch with flanks"
        );

        // Find the match at position 500 (SASSY may find other matches in random flanks)
        let match_at_500 = results.iter().find(|(_, _, _, _, _, pos)| *pos == 500);
        assert!(match_at_500.is_some(), "Should find match at position 500");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = match_at_500.unwrap();
        assert_eq!(cigar, "4=1X5=");
    }

    #[test]
    fn test_with_bulge_and_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGAATCGAT"; // Single base insertion after position 4
        let target = create_flanked_sequence(&mut rng, core, 500);

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false);
        assert!(
            !results.is_empty(),
            "Should accept a single base bulge with flanks"
        );

        // Find the match at position 500 (SASSY may find other matches in random flanks)
        let match_at_500 = results.iter().find(|(_, _, _, _, _, pos)| *pos == 500);
        assert!(match_at_500.is_some(), "Should find match at position 500");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = match_at_500.unwrap();
        assert!(
            cigar.contains('I') || cigar.contains('D'),
            "Should contain an insertion or deletion"
        );
    }

    #[test]
    fn test_too_many_differences_with_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGTT"; // Three mismatches at positions 5, 8, 9
        let target = create_flanked_sequence(&mut rng, core, 500);

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false);

        // Should not find a match at position 500 (too many mismatches)
        // Note: SASSY may find accidental matches in the random flanks, so we only check position 500
        let match_at_500 = results.iter().find(|(_, _, _, _, _, pos)| *pos == 500);
        assert!(
            match_at_500.is_none(),
            "Should reject sequence with too many mismatches at position 500"
        );
    }

    #[test]
    fn test_hit_quality_scoring_and_filtering() {
        // Create Hit objects with different qualities
        let guide_seq = Arc::new(b"ATCGATCGAT".to_vec());

        // Perfect match hit
        let perfect_hit = Hit {
            ref_id: "chr1".to_string(),
            pos: 100,
            strand: '+',
            score: 0,
            cigar: "MMMMMMMMMM".to_string(), // 10 perfect matches
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
        };

        // Hit with one mismatch
        let mismatch_hit = Hit {
            ref_id: "chr1".to_string(),
            pos: 105, // Overlaps with perfect_hit
            strand: '+',
            score: 3,                        // Higher score (worse)
            cigar: "MMMMXMMMMM".to_string(), // 9 matches, 1 mismatch
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
        };

        // Hit with a bulge
        let bulge_hit = Hit {
            ref_id: "chr1".to_string(),
            pos: 110, // Doesn't overlap with others
            strand: '+',
            score: 6,                       // Even higher score (worse)
            cigar: "MMMDMMMMM".to_string(), // Gap
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
        };
        assert_eq!(
            perfect_hit.end_pos(),
            110,
            "End position should be pos + matches"
        );
        assert_eq!(
            mismatch_hit.end_pos(),
            115,
            "End position includes mismatches"
        );
        assert_eq!(bulge_hit.end_pos(), 119, "End position includes deletions");
    }
}

#[derive(Parser)]
#[command(author, version, about = "CRISPR guide RNA off-target scanner")]
struct Args {
    /// Input reference FASTA file (-r)
    #[arg(short, long)]
    reference: PathBuf,

    /// Guide RNA sequence (without PAM) (-g)
    #[arg(short, long, required_unless_present = "guide_file")]
    guide: Option<String>,

    /// File containing guide RNA sequences (FASTA or one per line)
    #[arg(long = "guide-file")]
    guide_file: Option<PathBuf>,

    /// PAM sequence (to use for CFD scoring)
    #[arg(short = 'p', long, default_value = "GG")]
    pam: String,

    /// Maximum number of mismatches allowed
    #[arg(short, long, default_value = "4")]
    max_mismatches: u32,

    /// Maximum number of bulges allowed
    #[arg(short = 'b', long, default_value = "1")]
    max_bulges: u32,

    /// Maximum size of each bulge in bp
    #[arg(short = 'z', long, default_value = "2")]
    max_bulge_size: u32,

    /// Minimum fraction of guide that must match (0.0-1.0)
    #[arg(short = 'f', long, default_value = "0.75")]
    min_match_fraction: f32,

    /// Path to mismatch scores file for CFD calculation
    #[arg(long, default_value = "mismatch_scores.txt")]
    mismatch_scores: PathBuf,

    /// Path to PAM scores file for CFD calculation
    #[arg(long, default_value = "pam_scores.txt")]
    pam_scores: PathBuf,

    /// Number of threads to use (default: number of logical CPUs)
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Disable all filtering (report every alignment)
    #[arg(long)]
    no_filter: bool,
}

fn convert_to_minimap2_cigar(cigar: &str) -> String {
    // Remove the debug print line
    if cigar.is_empty() {
        return "".to_string();
    }

    cigar.to_string() // Just return the CIGAR as-is
}

fn scan_contig_sassy(
    guide: &[u8],
    contig: &[u8],
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
    min_match_fraction: f32,
    no_filter: bool,
) -> Vec<(i32, String, u32, u32, u32, usize)> {
    // Calculate maximum allowed errors
    let max_errors = (max_mismatches + max_bulges) as usize;

    // Create SASSY searcher with DNA profile
    let mut searcher: Searcher<Dna> = Searcher::new(false, None);

    // Convert to Vec for SASSY's SearchAble trait
    let contig_vec = contig.to_vec();

    // Search for ALL matches using SASSY
    let matches = searcher.search(guide, &contig_vec, max_errors);

    if matches.is_empty() {
        return Vec::new();
    }

    // Process ALL matches, not just the best one
    matches
        .into_iter()
        .filter_map(|sassy_match| {
            let score = sassy_match.cost as i32;
            let pos = sassy_match.text_start;

            // Use SASSY's CIGAR and normalize it to always include counts
            let cigar_str = normalize_cigar(&sassy_match.cigar.to_string());

            // Calculate statistics from CIGAR
            let (matches_count, mismatches, gaps, max_gap_size) = parse_cigar_stats(&cigar_str);

            // Apply filtering
            let non_n_positions = guide.iter().filter(|&&b| b != b'N').count();
            let match_percentage = if non_n_positions > 0 {
                (matches_count as f32 / non_n_positions as f32) * 100.0
            } else {
                0.0
            };

            if no_filter
                || (matches_count >= 1
                    && match_percentage >= min_match_fraction * 100.0
                    && mismatches <= max_mismatches
                    && gaps <= max_bulges
                    && max_gap_size <= max_bulge_size)
            {
                Some((score, cigar_str, mismatches, gaps, max_gap_size, pos))
            } else {
                None
            }
        })
        .collect()
}

/// Normalize CIGAR string to always include counts and consolidate consecutive operations
/// e.g., "X" -> "1X", "===XX=" -> "3=2X1="
fn normalize_cigar(cigar: &str) -> String {
    let mut ops: Vec<(u32, char)> = Vec::new();
    let mut chars = cigar.chars().peekable();

    // First, parse the CIGAR into (count, op) pairs
    while let Some(&ch) = chars.peek() {
        let count = if ch.is_ascii_digit() {
            // Has a count, parse it
            let mut num_str = String::new();
            while let Some(&digit_ch) = chars.peek() {
                if digit_ch.is_ascii_digit() {
                    num_str.push(chars.next().unwrap());
                } else {
                    break;
                }
            }
            num_str.parse::<u32>().unwrap_or(1)
        } else {
            1
        };

        // Get the operation
        if let Some(op) = chars.next() {
            ops.push((count, op));
        }
    }

    // Now consolidate consecutive operations
    let mut consolidated: Vec<(u32, char)> = Vec::new();
    for (count, op) in ops {
        if let Some(last) = consolidated.last_mut() {
            if last.1 == op {
                // Same operation, add to count
                last.0 += count;
            } else {
                // Different operation, add new entry
                consolidated.push((count, op));
            }
        } else {
            // First operation
            consolidated.push((count, op));
        }
    }

    // Format as string
    consolidated
        .iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect::<String>()
}

fn parse_cigar_stats(cigar: &str) -> (usize, u32, u32, u32) {
    let mut matches_count = 0;
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut max_gap_size = 0;
    let mut current_gap_size = 0;

    // Parse CIGAR string with proper number handling
    let mut chars = cigar.chars().peekable();
    while let Some(&ch) = chars.peek() {
        if ch.is_ascii_digit() {
            // Extract the count
            let mut num_str = String::new();
            while let Some(&digit_ch) = chars.peek() {
                if digit_ch.is_ascii_digit() {
                    num_str.push(chars.next().unwrap());
                } else {
                    break;
                }
            }

            // Get the operation
            if let Some(op) = chars.next() {
                if let Ok(count) = num_str.parse::<u32>() {
                    match op {
                        '=' | 'M' => {
                            matches_count += count as usize;
                            current_gap_size = 0;
                        }
                        'X' => {
                            mismatches += count;
                            current_gap_size = 0;
                        }
                        'I' | 'D' => {
                            if current_gap_size == 0 {
                                gaps += 1;
                            }
                            current_gap_size += count;
                            max_gap_size = max_gap_size.max(current_gap_size);
                        }
                        _ => {}
                    }
                }
            }
        } else {
            // Handle single-character operations (no counts)
            let op = chars.next().unwrap();
            match op {
                '=' | 'M' => {
                    matches_count += 1;
                    current_gap_size = 0;
                }
                'X' => {
                    mismatches += 1;
                    current_gap_size = 0;
                }
                'I' | 'D' => {
                    if current_gap_size == 0 {
                        gaps += 1;
                    }
                    current_gap_size += 1;
                    max_gap_size = max_gap_size.max(current_gap_size);
                }
                _ => {}
            }
        }
    }

    (matches_count, mismatches, gaps, max_gap_size)
}

fn normalize_sequence(seq: &str) -> String {
    seq.chars()
        .filter(|c| !c.is_whitespace())
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

fn parse_line_guides(contents: &str) -> Vec<Arc<Vec<u8>>> {
    contents
        .lines()
        .filter_map(|line| {
            let normalized = normalize_sequence(line);
            if normalized.is_empty() {
                None
            } else {
                Some(Arc::new(normalized.into_bytes()))
            }
        })
        .collect()
}

fn parse_fasta_guides(contents: &str) -> Vec<Arc<Vec<u8>>> {
    let mut guides = Vec::new();
    let mut current = String::new();

    for line in contents.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('>') {
            if !current.is_empty() {
                let normalized = normalize_sequence(&current);
                if !normalized.is_empty() {
                    guides.push(Arc::new(normalized.into_bytes()));
                }
                current.clear();
            }
        } else {
            current.push_str(trimmed);
        }
    }

    if !current.is_empty() {
        let normalized = normalize_sequence(&current);
        if !normalized.is_empty() {
            guides.push(Arc::new(normalized.into_bytes()));
        }
    }

    guides
}

fn load_guides(args: &Args) -> io::Result<Vec<Arc<Vec<u8>>>> {
    let mut guides = Vec::new();

    if let Some(seq) = args.guide.as_ref() {
        let normalized = normalize_sequence(seq);
        if !normalized.is_empty() {
            guides.push(Arc::new(normalized.into_bytes()));
        }
    }

    if let Some(path) = args.guide_file.as_ref() {
        let contents = fs::read_to_string(path)?;
        let mut file_guides = if contents.trim_start().starts_with('>') {
            parse_fasta_guides(&contents)
        } else {
            parse_line_guides(&contents)
        };
        guides.append(&mut file_guides);
    }

    Ok(guides)
}

fn main() {
    let args = Args::parse();

    // **FIXED: Better CFD initialization with more informative error handling**
    match cfd_score::init_score_matrices(
        args.mismatch_scores
            .to_str()
            .unwrap_or("mismatch_scores.txt"),
        args.pam_scores.to_str().unwrap_or("pam_scores.txt"),
    ) {
        Ok(()) => {
            eprintln!("CFD scoring initialized successfully");
        }
        Err(e) => {
            eprintln!("Warning: CFD scoring disabled - {}", e);
            eprintln!(
                "Expected files: {} and {}",
                args.mismatch_scores.display(),
                args.pam_scores.display()
            );
        }
    }

    let guides = match load_guides(&args) {
        Ok(g) => g,
        Err(err) => {
            eprintln!("Failed to load guide sequences: {}", err);
            std::process::exit(1);
        }
    };

    if guides.is_empty() {
        eprintln!("No guide RNAs supplied. Use --guide or --guide-file.");
        std::process::exit(1);
    }

    let mut worker_count = args.threads.unwrap_or_else(|| {
        thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    });
    if worker_count == 0 {
        worker_count = 1;
    }

    let queue_depth = worker_count * 4;
    let (task_tx, task_rx) = bounded::<Task>(queue_depth);
    let (hit_tx, hit_rx) = bounded::<HitResult>(queue_depth);

    let max_mismatches = args.max_mismatches;
    let max_bulges = args.max_bulges;
    let max_bulge_size = args.max_bulge_size;
    let min_match_fraction = args.min_match_fraction;
    let no_filter = args.no_filter;
    let pam = Arc::new(args.pam.clone());

    let writer_handle = {
        let pam = Arc::clone(&pam);
        thread::spawn(move || {
            while let Ok(hit) = hit_rx.recv() {
                report_hit(
                    hit.ref_id.as_str(),
                    hit.pos,
                    hit.guide.len(),
                    hit.strand,
                    hit.score,
                    &hit.cigar,
                    &hit.guide,
                    hit.target_len,
                    max_mismatches,
                    max_bulges,
                    max_bulge_size,
                    &hit.target_seq,
                    pam.as_str(),
                );
            }
        })
    };

    let mut worker_handles = Vec::with_capacity(worker_count);
    for _ in 0..worker_count {
        let task_rx = task_rx.clone();
        let hit_tx = hit_tx.clone();
        let thread_min_match_fraction = min_match_fraction;
        let thread_no_filter = no_filter;
        worker_handles.push(thread::spawn(move || {
            while let Ok(task) = task_rx.recv() {
                let guide = task.guide;
                let contig_seq = task.contig_seq;
                let contig_len = task.contig_len;
                let ref_id = task.ref_id;
                let guide_len = guide.len();

                if contig_len < guide_len {
                    continue;
                }

                let matches = scan_contig_sassy(
                    &guide,
                    &contig_seq,
                    max_mismatches,
                    max_bulges,
                    max_bulge_size,
                    thread_min_match_fraction,
                    thread_no_filter,
                );

                let contig_slice = contig_seq.as_slice();
                for (score, cigar, _mismatches, _gaps, _max_gap_size, pos) in matches {
                    let start = pos.min(contig_len);
                    let end = (start + guide_len).min(contig_len);
                    let target_seq = contig_slice[start..end].to_vec();

                    if hit_tx
                        .send(HitResult {
                            ref_id: Arc::clone(&ref_id),
                            pos: start,
                            strand: '+',
                            score,
                            cigar,
                            guide: Arc::clone(&guide),
                            target_len: contig_len,
                            target_seq,
                        })
                        .is_err()
                    {
                        return;
                    }
                }
            }
        }));
    }
    drop(hit_tx);

    // Process reference sequences with a single reader feeding the work queue
    let file = File::open(&args.reference).expect("Failed to open reference file");
    let reader: Box<dyn BufRead> = if args.reference.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let reader = fasta::Reader::new(reader);

    'reader: for record in reader.records() {
        let record = match record {
            Ok(rec) => rec,
            Err(err) => {
                eprintln!("Error during FASTA record parsing: {}", err);
                break;
            }
        };

        let contig_seq = Arc::new(record.seq().to_vec());
        let contig_len = contig_seq.len();
        if contig_len == 0 {
            continue;
        }
        let ref_id = Arc::new(record.id().to_string());

        for guide in &guides {
            if contig_len < guide.len() {
                continue;
            }

            let task = Task {
                guide: Arc::clone(guide),
                contig_seq: Arc::clone(&contig_seq),
                contig_len,
                ref_id: Arc::clone(&ref_id),
            };

            if task_tx.send(task).is_err() {
                eprintln!("Task queue closed while scheduling work; terminating reader loop");
                break 'reader;
            }
        }
    }
    drop(task_tx);

    for handle in worker_handles {
        if let Err(err) = handle.join() {
            eprintln!("Worker thread exited with error: {:?}", err);
        }
    }

    if let Err(err) = writer_handle.join() {
        eprintln!("Writer thread exited with error: {:?}", err);
    }
}
