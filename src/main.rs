use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader};
use flate2::read::MultiGzDecoder;
use std::sync::Arc;
use std::collections::HashMap;
use clap::Parser;
use bio::io::fasta;
use sassy::profiles::Dna;
use sassy::search::Searcher;
// Remove the broken imports for now - we'll add correct ones later
// use sassy::{search, Alphabet, SearchConfig};
use std::fmt::Write;
use rayon::prelude::*;

mod cfd_score;

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        _ => b'N',  // Convert any unexpected bases to N
    }).collect()
}

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
    cfd_score: Option<f64>,  // Add CFD score field
    target_seq: Vec<u8>,     // Add target sequence for CFD calculation
}

impl Hit {
    fn quality_score(&self) -> i32 {
        // Count matches and other stats
        let matches = self.cigar.chars()
            .filter(|&c| c == 'M' || c == '=')
            .count();
            
        let mismatches = self.cigar.chars()
            .filter(|&c| c == 'X')
            .count();
            
        let gaps = self.cigar.chars()
            .filter(|&c| c == 'I' || c == 'D')
            .count();
            
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
        self.strand == other.strand && 
        self.ref_id == other.ref_id &&
        self.pos < other.end_pos() && 
        other.pos < self.end_pos()
    }
}

// Replace your entire report_hit function with this corrected version:

fn report_hit(ref_id: &str, pos: usize, _len: usize, strand: char, 
              _score: i32, cigar: &str, guide: &[u8], target_len: usize,
              _max_mismatches: u32, _max_bulges: u32, _max_bulge_size: u32,
              target_seq: &[u8], pam: &str) {
    
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
                        },
                        'X' => {
                            mismatches += count;
                        },
                        'I' | 'D' => {
                            gaps += 1;
                            max_gap_size = max_gap_size.max(count);
                        },
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
                },
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
        None => "\tcf:f:0.0000".to_string()
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
use rand::{SeedableRng, RngCore, rngs::SmallRng};

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
    fn test_perfect_match_sassy() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGATCGAT";
        
        let result = scan_window_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_some());
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "10=");
    }

    #[test]
    fn test_with_mismatches_sassy() {
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGAT";  // Single mismatch at position 5
        
        let result = scan_window_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single mismatch");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "4=1X5=");
    }

    #[test]
    fn test_with_bulge_sassy() {
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGAATCGAT";  // Single base insertion after position 4
        
        let result = scan_window_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single base bulge");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences_sassy() {
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        
        let result = scan_window_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_none());
    }

    #[test]
    fn test_perfect_match_with_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let target = create_flanked_sequence(&mut rng, guide, 500);
        
        let result = scan_window_sassy(guide, &target[500..510], 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should match perfectly even with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "10=");
    }

    #[test]
    fn test_with_mismatches_and_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGAT";  // Single mismatch at position 5
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window_sassy(guide, &target[500..510], 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single mismatch with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "4=1X5=");
    }
    
    #[test]
    fn test_with_bulge_and_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGAATCGAT";  // Single base insertion after position 4
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window_sassy(guide, &target[500..511], 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single base bulge with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences_with_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window_sassy(guide, &target[500..510], 1, 1, 1, 0.75, false);
        assert!(result.is_none(), "Should reject sequence with too many mismatches even with flanks");
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
            cigar: "MMMMMMMMMM".to_string(),  // 10 perfect matches
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
            pos: 105,  // Overlaps with perfect_hit
            strand: '+',
            score: 3,  // Higher score (worse)
            cigar: "MMMMXMMMMM".to_string(),  // 9 matches, 1 mismatch
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
            pos: 110,  // Doesn't overlap with others
            strand: '+',
            score: 6,  // Even higher score (worse)
            cigar: "MMMDMMMMM".to_string(),  // Gap
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
        };
        assert_eq!(perfect_hit.end_pos(), 110, "End position should be pos + matches");
        assert_eq!(mismatch_hit.end_pos(), 115, "End position includes mismatches");
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
    #[arg(short, long)]
    guide: String,

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

    /// Size of sequence window to scan (bp, default: 4x guide length)
    #[arg(short = 'w', long)]
    window_size: Option<usize>,

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
    
    cigar.to_string()  // Just return the CIGAR as-is
}

fn scan_window_sassy(
    guide: &[u8], 
    window: &[u8], 
    max_mismatches: u32, 
    max_bulges: u32, 
    max_bulge_size: u32,
    min_match_fraction: f32, 
    no_filter: bool
) -> Option<(i32, String, u32, u32, u32, usize)> {
    
    // Calculate maximum allowed errors
    let max_errors = (max_mismatches + max_bulges) as usize;
    
    // Create SASSY searcher with DNA profile
    let mut searcher: Searcher<Dna> = Searcher::new(false, None);
    
    // Convert window to a Vec so it implements SearchAble
    let window_vec = window.to_vec();
    
    // Search for matches using real SASSY
    let matches = searcher.search(guide, &window_vec, max_errors);

    if matches.is_empty() {
        return None;
    }
    
    // Take the best match (lowest cost)
    let best_match = matches.into_iter().min_by_key(|m| m.cost)?;
    
    let score = best_match.cost as i32;
    
    // Simple CIGAR generation based on alignment cost
    let cigar_str = if best_match.cost == 0 {
        format!("{}=", guide.len())
    } else {
        let matches = guide.len().saturating_sub(best_match.cost as usize);
        let mismatches = best_match.cost as usize;
        
        if matches > 0 && mismatches > 0 {
            format!("{}={}X", matches, mismatches)
        } else if mismatches > 0 {
            format!("{}X", mismatches)
        } else {
            format!("{}=", guide.len())
        }
    };

    // Calculate statistics from CIGAR
    let (matches_count, mismatches, gaps, max_gap_size) = parse_cigar_stats(&cigar_str);

    // Apply filtering
    let non_n_positions = guide.iter().filter(|&&b| b != b'N').count();
    let match_percentage = if non_n_positions > 0 {
        (matches_count as f32 / non_n_positions as f32) * 100.0
    } else {
        0.0
    };

    if no_filter || (
        matches_count >= 1 && 
        match_percentage >= min_match_fraction * 100.0 && 
        mismatches <= max_mismatches && 
        gaps <= max_bulges && 
        max_gap_size <= max_bulge_size
    ) {
        // DEBUG: Show what SASSY found
        let actual_match_pos = best_match.start.1 as usize;
    

    
        // Show the actual sequences being compared
        if actual_match_pos + guide.len() <= window.len() {
            let found_seq = &window[actual_match_pos..actual_match_pos + guide.len()];
        } else {
            println!("ERROR: Position {} + {} > window length {}", actual_match_pos, guide.len(), window.len());
        }
    
        Some((score, cigar_str, mismatches, gaps, max_gap_size, actual_match_pos))
    } else {
        None
    }
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
                        },
                        'X' => {
                            mismatches += count;
                            current_gap_size = 0;
                        },
                        'I' | 'D' => {
                            if current_gap_size == 0 {
                                gaps += 1;
                            }
                            current_gap_size += count;
                            max_gap_size = max_gap_size.max(current_gap_size);
                        },
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
                },
                'X' => {
                    mismatches += 1;
                    current_gap_size = 0;
                },
                'I' | 'D' => {
                    if current_gap_size == 0 {
                        gaps += 1;
                    }
                    current_gap_size += 1;
                    max_gap_size = max_gap_size.max(current_gap_size);
                },
                _ => {}
            }
        }
    }
    
    (matches_count, mismatches, gaps, max_gap_size)
}

fn main() {
    let args = Args::parse();
    
    // **FIXED: Better CFD initialization with more informative error handling**
    match cfd_score::init_score_matrices(
        args.mismatch_scores.to_str().unwrap_or("mismatch_scores.txt"),
        args.pam_scores.to_str().unwrap_or("pam_scores.txt")
    ) {
        Ok(()) => {
            eprintln!("CFD scoring initialized successfully");
        }
        Err(e) => {
            eprintln!("Warning: CFD scoring disabled - {}", e);
            eprintln!("Expected files: {} and {}", 
                      args.mismatch_scores.display(), 
                      args.pam_scores.display());
        }
    }
    
    // Prepare guide sequences (forward and reverse complement)
    let guide_fwd = Arc::new(args.guide.as_bytes().to_vec());
    let _guide_rc = Arc::new(reverse_complement(&guide_fwd));
    let guide_len = guide_fwd.len();

    // Set thread pool size if specified
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to initialize thread pool");
    }

    // Process reference sequences
    let file = File::open(&args.reference).expect("Failed to open reference file");
    let reader: Box<dyn BufRead> = if args.reference.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let reader = fasta::Reader::new(reader);
    
    for result in reader.records() {
        let record = result.expect("Error during FASTA record parsing");
        let seq = record.seq().to_vec();
        let seq_len = seq.len();
        let record_id = record.id().to_string();
        
        // Calculate window size as 4x guide length if not specified
        let window_size = args.window_size.unwrap_or(guide_len * 4);
        let step_size = window_size / 2;
        let windows: Vec<_> = (0..seq.len())
            .step_by(step_size)
            .map(|i| (i, (i + window_size).min(seq.len())))
            .collect();

        // Process windows in parallel and collect all hits
        let hits: Vec<Hit> = windows.into_par_iter()
            .map_init(
                || (),
                |_unit, (window_start, end)| {
                    let window = &seq[window_start..end];
                    if window.len() < guide_len { return None; }
    
                    // Try forward orientation
                    if let Some((score, cigar, _mismatches, _gaps, _max_gap_size, match_offset_in_window)) = 
                        scan_window_sassy(&guide_fwd, window,
                                args.max_mismatches, args.max_bulges, args.max_bulge_size,
                                args.min_match_fraction, args.no_filter) {
                
                        // DEBUG: Show the coordinate calculation

                        let actual_pos = window_start + match_offset_in_window;
                
                        return Some(Hit {
                            ref_id: record_id.clone(),
                            pos: actual_pos,
                            strand: '+',
                            score,
                            cigar: cigar.clone(),
                            guide: Arc::clone(&guide_fwd),
                            target_len: seq_len,
                            max_mismatches: args.max_mismatches,
                            max_bulges: args.max_bulges,
                            max_bulge_size: args.max_bulge_size,
                            cfd_score: None,
                            target_seq: {
                                // Extract target sequence from the correct position
                                let start = actual_pos;
                                let end = (actual_pos + guide_len).min(seq_len);
                                seq[start..end].to_vec()
                            },
                        });
                    }
            
                    None
                })
            .filter_map(|x| x)
            .collect();

        // Report hits directly (simplified for now)
        for hit in hits {
            report_hit(
                &hit.ref_id, 
                hit.pos, 
                hit.guide.len(), 
                hit.strand, 
                hit.score, 
                &hit.cigar, 
                &hit.guide, 
                hit.target_len,
                hit.max_mismatches, 
                hit.max_bulges, 
                hit.max_bulge_size,
                &hit.target_seq,
                &args.pam
            );
        }
    }
}
