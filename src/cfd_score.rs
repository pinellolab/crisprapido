//! # Cutting frequency determination (CFD) score calculator
//! Module for calculating CFD scores for CRISPR guide RNA off-target sites
//! Adapted from the Python implementation by Linda Lin 3/23/2025

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;
use std::sync::Once;
use std::sync::Mutex;
use lazy_static::lazy_static;

// Static matrices for CFD scoring
lazy_static! {
    static ref MISMATCH_SCORES: Mutex<Option<HashMap<String, f64>>> = Mutex::new(None);
    static ref PAM_SCORES: Mutex<Option<HashMap<String, f64>>> = Mutex::new(None);
    static ref INIT: Once = Once::new();
}

/// Initialize the scoring matrices from the provided file paths
pub fn init_score_matrices(mismatch_path: &str, pam_path: &str) -> Result<(), String> {
    INIT.call_once(|| {
        let mm_matrix = parse_scoring_matrix(mismatch_path)
            .map_err(|e| format!("Failed to load mismatch scores: {}", e));
        
        let pam_matrix = parse_scoring_matrix(pam_path)
            .map_err(|e| format!("Failed to load PAM scores: {}", e));
        
        if let (Ok(mm), Ok(pam)) = (mm_matrix, pam_matrix) {
            *MISMATCH_SCORES.lock().unwrap() = Some(mm);
            *PAM_SCORES.lock().unwrap() = Some(pam);
        }
    });
    
    // Check if matrices were successfully loaded
    let mm_loaded = MISMATCH_SCORES.lock().unwrap().is_some();
    let pam_loaded = PAM_SCORES.lock().unwrap().is_some();
    
    if mm_loaded && pam_loaded {
        Ok(())
    } else {
        Err("Failed to initialize scoring matrices".to_string())
    }
}

/// Calculate CFD score for aligned sequences
pub fn calculate_cfd(spacer: &str, protospacer: &str, pam: &str) -> Result<f64, String> {
    
    // Handle different guide lengths by taking first 20bp
    let spacer_20bp = if spacer.len() >= 20 {
        if spacer.contains('-') {
            // Handle gaps - for now, truncate to 20 characters including gaps
            let truncated = if spacer.len() > 20 { &spacer[0..20] } else { spacer };
            truncated.to_string()
        } else {
            spacer[0..20].to_string() // Take first 20bp
        }
    } else {
        return Err(format!("Spacer too short: {} bp, expected at least 20 bp", spacer.len()));
    };

    let protospacer_20bp = if protospacer.len() >= 20 {
        if protospacer.contains('-') {
            let truncated = if protospacer.len() > 20 { &protospacer[0..20] } else { protospacer };
            truncated.to_string()
        } else {
            protospacer[0..20].to_string() // Take first 20bp
        }
    } else {
        return Err(format!("Protospacer too short: {} bp, expected at least 20 bp", protospacer.len()));
    };


    // Validate PAM length (allow NGG-style 3-mers by dropping the first base)
    let pam_clean = pam.to_uppercase();
    let pam_processed = match pam_clean.len() {
        2 => pam_clean,
        3 => pam_clean[1..].to_string(),
        _ => return Err(format!("PAM must be 2 or 3 nucleotides, got {} bp", pam.len())),
    };

    // Get locked references to scoring matrices

    let mm_scores_lock = MISMATCH_SCORES.lock().unwrap();
    let pam_scores_lock = PAM_SCORES.lock().unwrap();

    // Verify matrices are initialized
    let mm_scores = mm_scores_lock.as_ref()
        .ok_or_else(|| "Mismatch scores not initialized".to_string())?;
    let pam_scores = pam_scores_lock.as_ref()
        .ok_or_else(|| "PAM scores not initialized".to_string())?;

    // Pre-process sequences (convert T to U for RNA)
    let spacer_list: Vec<char> = spacer_20bp.to_uppercase().replace("T", "U").chars().collect();
    let protospacer_list: Vec<char> = protospacer_20bp.to_uppercase().replace("T", "U").chars().collect();

    // Ensure both sequences are exactly 20bp after processing
    if spacer_list.len() != 20 || protospacer_list.len() != 20 {
        return Err(format!("Processed sequences must be 20bp: spacer={}, protospacer={}", 
                          spacer_list.len(), protospacer_list.len()));
    }

    // Calculate CFD score
    let mut score = 1.0;
    
    for (i, (&spacer_nt, &proto_nt)) in spacer_list.iter().zip(protospacer_list.iter()).enumerate() {
        if spacer_nt == proto_nt {
            continue;
        }

        if is_ambiguous_nt(spacer_nt) || is_ambiguous_nt(proto_nt) {
            continue;
        }

        if i == 0 && (spacer_nt == '-' || proto_nt == '-') {
            continue;
        }

        let key = format!("r{}:d{},{}", spacer_nt, reverse_complement_nt(proto_nt), i + 1);
        if let Some(penalty) = mm_scores.get(&key) {
            score *= penalty;
        }
    }

    // Apply PAM penalty
    match pam_scores.get(&pam_processed) {
        Some(pam_penalty) => {
            score *= pam_penalty;
        },
        None => {
            return Ok(0.0); // Unknown PAM gets score 0
        }
    }

    Ok(score)
}

/// Get CFD score using CIGAR-based alignment
/// 
/// # Arguments
/// * `guide` - Guide RNA sequence as byte array
/// * `target` - Target DNA sequence as byte array
/// * `cigar` - CIGAR string representing the alignment
/// * `pam` - 2nt PAM sequence
/// 
/// # Returns
/// * `Option<f64>` - CFD score if calculation succeeds

/// Prepare aligned spacer and protospacer sequences for CFD calculation
fn prepare_aligned_sequences(guide: &[u8], target: &[u8], cigar: &str) -> (String, String) {
    let mut spacer = String::new();
    let mut protospacer = String::new();
    
    // Handle empty CIGAR by assuming perfect match
    if cigar.is_empty() {
        let guide_str = String::from_utf8_lossy(guide);
        let target_str = String::from_utf8_lossy(target);
        
        // Take first 20bp of each sequence
        let spacer_20 = if guide_str.len() >= 20 { &guide_str[0..20] } else { &guide_str };
        let target_20 = if target_str.len() >= 20 { &target_str[0..20] } else { &target_str };
        
        return (spacer_20.to_string(), target_20.to_string());
    }
    
    let mut guide_pos = 0;
    let mut target_pos = 0;
    
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
                if let Ok(count) = num_str.parse::<usize>() {
                    match op {
                        'M' | '=' => {
                            // Match operations
                            for _ in 0..count {
                                if guide_pos < guide.len() && target_pos < target.len() {
                                    spacer.push(guide[guide_pos] as char);
                                    protospacer.push(target[target_pos] as char);
                                    guide_pos += 1;
                                    target_pos += 1;
                                } else {
                                    break;
                                }
                            }
                        },
                        'X' => {
                            // Mismatch operations
                            for _ in 0..count {
                                if guide_pos < guide.len() && target_pos < target.len() {
                                    spacer.push(guide[guide_pos] as char);
                                    protospacer.push(target[target_pos] as char);
                                    guide_pos += 1;
                                    target_pos += 1;
                                } else {
                                    break;
                                }
                            }
                        },
                        'I' => {
                            // Insertion in query (gap in target)
                            for _ in 0..count {
                                if guide_pos < guide.len() {
                                    spacer.push(guide[guide_pos] as char);
                                    protospacer.push('-');
                                    guide_pos += 1;
                                } else {
                                    break;
                                }
                            }
                        },
                        'D' => {
                            // Deletion in query (gap in query) 
                            for _ in 0..count {
                                if target_pos < target.len() {
                                    spacer.push('-');
                                    protospacer.push(target[target_pos] as char);
                                    target_pos += 1;
                                } else {
                                    break;
                                }
                            }
                        },
                        _ => {
                            eprintln!("Warning: Unknown CIGAR operation: {}", op);
                        }
                    }
                }
            }
        } else {
            // Handle single character operations (legacy format)
            let op = chars.next().unwrap();
            match op {
                'M' | '=' | 'X' => {
                    if guide_pos < guide.len() && target_pos < target.len() {
                        spacer.push(guide[guide_pos] as char);
                        protospacer.push(target[target_pos] as char);
                        guide_pos += 1;
                        target_pos += 1;
                    }
                },
                'I' => {
                    if guide_pos < guide.len() {
                        spacer.push(guide[guide_pos] as char);
                        protospacer.push('-');
                        guide_pos += 1;
                    }
                },
                'D' => {
                    if target_pos < target.len() {
                        spacer.push('-');
                        protospacer.push(target[target_pos] as char);
                        target_pos += 1;
                    }
                },
                _ => {}
            }
        }
    }
    
    // Ensure we have exactly 20 characters by padding or truncating
    while spacer.len() < 20 {
        spacer.push('-');
    }
    while protospacer.len() < 20 {
        protospacer.push('-');
    }
    
    // Truncate to exactly 20bp
    let spacer_final = spacer.chars().take(20).collect();
    let protospacer_final = protospacer.chars().take(20).collect();
    
    (spacer_final, protospacer_final)
}


/// Get reverse complement of a single nucleotide (supports bulges)
fn reverse_complement_nt(nucleotide: char) -> char {
    match nucleotide {
        'A' => 'T',
        'C' => 'G',
        'T' | 'U' => 'A',
        'G' => 'C',
        '-' => '-',
        _ => nucleotide,
    }
}

fn is_ambiguous_nt(nt: char) -> bool {
    matches!(nt, 'N' | 'n' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V')
}

/// Parse scoring matrix from space-delimited file
fn parse_scoring_matrix(file_path: &str) -> Result<HashMap<String, f64>, String> {
    // Open file
    let file = File::open(file_path)
        .map_err(|e| format!("Cannot open {}: {}", file_path, e))?;
    
    // Read file
    let reader = BufReader::new(file);
    let mut matrix = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading line: {}", e))?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let score = parts[1].parse::<f64>()
                .map_err(|e| format!("Invalid score format: {}", e))?;
            matrix.insert(parts[0].to_string(), score);
        }
    }
    Ok(matrix)
}

#[cfg(test)]
mod cfd_unit_tests {
    use super::*;

    /// Test perfect match with canonical PAM gives score 1.0
    #[test]
    fn test_perfect_match() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let spacer = "GAAACAGTCGATTTTATCAC";
        let score = calculate_cfd(spacer, spacer, "GG").unwrap();
        assert!((score - 1.0).abs() < 1e-9, "Perfect match should be 1.0, got {}", score);
    }

    /// Test that non-canonical PAMs reduce score
    #[test]
    fn test_pam_effects() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let spacer = "GAAACAGTCGATTTTATCAC";

        // GG should be highest (1.0 for perfect match)
        let gg_score = calculate_cfd(spacer, spacer, "GG").unwrap();
        let ag_score = calculate_cfd(spacer, spacer, "AG").unwrap();
        let cg_score = calculate_cfd(spacer, spacer, "CG").unwrap();

        assert!(gg_score > ag_score, "GG should score higher than AG");
        assert!(ag_score > cg_score, "AG should score higher than CG");
    }

    /// Test gap at position 1 (PAM-distal) has no penalty
    #[test]
    fn test_position_1_gap_no_penalty() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        // Gap in spacer at position 1
        let score = calculate_cfd("-AAACAGTCGATTTTATCAC", "GAAACAGTCGATTTTATCAC", "GG").unwrap();
        assert!((score - 1.0).abs() < 1e-9, "Gap at position 1 should have no penalty, got {}", score);

        // Gap in protospacer at position 1
        let score = calculate_cfd("GAAACAGTCGATTTTATCAC", "-AAACAGTCGATTTTATCAC", "GG").unwrap();
        assert!((score - 1.0).abs() < 1e-9, "Gap at position 1 should have no penalty, got {}", score);
    }

    /// Test gaps at other positions have penalties
    #[test]
    fn test_gaps_at_other_positions_have_penalties() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        // Gap at position 7 should have a penalty (score 0.0 per scoring matrix)
        let score = calculate_cfd("GAAACA-TCGATTTTATCAC", "GAAACAGTCGATTTTATCAC", "GG").unwrap();
        assert!((score - 0.0).abs() < 1e-9, "Gap at position 7 should be 0.0, got {}", score);

        // Gap at position 2 should reduce score
        let score = calculate_cfd("G-AACAGTCGATTTTATCAC", "GAAACAGTCGATTTTATCAC", "GG").unwrap();
        assert!(score < 1.0, "Gap at position 2 should reduce score, got {}", score);
    }

    /// Test case insensitivity
    #[test]
    fn test_case_insensitivity() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let upper = "GAAACAGTCGATTTTATCAC";
        let lower = "gaaacagtcgattttatcac";

        let score_upper = calculate_cfd(upper, upper, "GG").unwrap();
        let score_mixed = calculate_cfd(lower, upper, "gg").unwrap();

        assert!((score_upper - score_mixed).abs() < 1e-9,
            "Case should not matter: upper={}, mixed={}", score_upper, score_mixed);
    }

    /// Test T/U equivalence
    #[test]
    fn test_t_u_equivalence() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let with_t = "GAAACAGTCGATTTTATCAC";
        let with_u = "GAAACAGUCGAUUUUAUCAC";

        let score_t = calculate_cfd(with_t, with_t, "GG").unwrap();
        let score_u = calculate_cfd(with_u, with_t, "GG").unwrap();

        assert!((score_t - score_u).abs() < 1e-9,
            "T and U should be equivalent: T={}, U={}", score_t, score_u);
    }

    /// Test that scores are in valid range [0, 1]
    #[test]
    fn test_scores_in_valid_range() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let base = "GAAACAGTCGATTTTATCAC";
        let bases = ['A', 'C', 'G', 'T'];

        // Test a variety of mismatches
        for pos in [1, 5, 10, 15, 20] {
            for &new_base in &bases {
                let mut proto: Vec<char> = base.chars().collect();
                proto[pos - 1] = new_base;
                let proto_str: String = proto.into_iter().collect();

                let score = calculate_cfd(base, &proto_str, "GG").unwrap();
                assert!(score >= 0.0 && score <= 1.0,
                    "Score out of range [0,1]: {} for pos {} -> {}", score, pos, new_base);
            }
        }
    }

    /// Test CIGAR-based alignment preparation
    #[test]
    fn test_prepare_aligned_sequences() {
        let guide = b"ATCGATCGATCGATCGATCG";
        let target = b"ATCGATCGATCGATCGATCG";

        // Perfect match
        let (spacer, proto) = prepare_aligned_sequences(guide, target, "20=");
        assert_eq!(spacer.len(), 20);
        assert_eq!(proto.len(), 20);
        assert_eq!(spacer, proto);

        // Insertion (gap in target)
        let (spacer, proto) = prepare_aligned_sequences(guide, &target[..19], "10=1I9=");
        assert!(proto.contains('-'), "Insertion should create gap in protospacer");

        // Deletion (gap in query)
        let guide_short = b"ATCGATCGACGATCGATCG"; // 19bp
        let (spacer, proto) = prepare_aligned_sequences(guide_short, target, "10=1D9=");
        assert!(spacer.contains('-'), "Deletion should create gap in spacer");
    }
}

/// Get CFD score using CIGAR-based alignment
///
/// This function takes raw guide and target sequences along with a CIGAR string
/// and produces properly aligned sequences with gap characters for CFD scoring.
pub fn get_cfd_score(guide: &[u8], target_seq: &[u8], cigar: &str, pam: &str) -> Option<f64> {
    // Use CIGAR to create properly aligned sequences with gap characters
    let (spacer, protospacer) = prepare_aligned_sequences(guide, target_seq, cigar);

    match calculate_cfd(&spacer, &protospacer, pam) {
        Ok(score) => Some(score),
        Err(_) => None,
    }
}
