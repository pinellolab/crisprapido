extern crate crisprapido;

use crisprapido::cfd_score;
use std::fs;
use std::path::Path;

/// Ensure score files exist for testing
fn ensure_score_files() {
    // We'll skip the file creation since you already have the files
    // Make sure mismatch_scores.txt and pam_scores.txt are in the root directory
}

/// Test function to validate our CFD score calculation against Python implementation
#[test]
#[ignore]
fn test_cfd_score_against_python() {
    // Initialize scoring matrices
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");

    // Test case 1: Perfect match
    let spacer1 = "ATCGATCGATCGATCGATCG";
    let protospacer1 = "ATCGATCGATCGATCGATCG"; 
    let pam1 = "GG";
    
    println!("Test case 1: Perfect match");
    println!("Spacer:      {}", spacer1);
    println!("Protospacer: {}", protospacer1);
    println!("PAM:         {}", pam1);
    
    let score1 = cfd_score::calculate_cfd(spacer1, protospacer1, pam1)
        .expect("Perfect match calculation failed");
    println!("Score: {:.6} (expected 1.000000)", score1);
    
    assert!((score1 - 1.0).abs() < 0.0001, "Perfect match should be 1.0, got {}", score1);

    // Test case 2: Single mismatch at position 1 (A->T)
    let spacer2 = "ATCGATCGATCGATCGATCG";
    let protospacer2 = "TTCGATCGATCGATCGATCG";
    let pam2 = "GG";
    
    println!("\nTest case 2: Single mismatch at position 1");
    println!("Spacer:      {}", spacer2);
    println!("Protospacer: {}", protospacer2);
    println!("PAM:         {}", pam2);
    
    let score2 = cfd_score::calculate_cfd(spacer2, protospacer2, pam2)
        .expect("Mismatch calculation failed");
    println!("Score: {:.6} (expected < 1.0)", score2);
    
    // The score should be less than 1.0 for a mismatch
    assert!(score2 < 1.0, "Mismatch should result in score < 1.0, got {}", score2);
    
    // Test CIGAR-based calculation
    let guide2 = b"ATCGATCGATCGATCGATCG";
    let target2 = b"TTCGATCGATCGATCGATCG";
    let cigar2 = "1X19=";
    
    println!("\nTesting CIGAR-based calculation:");
    let cigar_score2 = cfd_score::get_cfd_score(guide2, target2, cigar2, pam2)
        .expect("CIGAR calculation failed");
    println!("CIGAR score: {:.6}", cigar_score2);
    
    // Both methods should give similar results
    assert!((score2 - cigar_score2).abs() < 0.001, 
            "Direct and CIGAR calculations should match: direct={}, cigar={}", 
            score2, cigar_score2);
}

/// Test with varied number of mismatches and their positions
#[test]
fn test_positional_effects() {
    // Initialize score matrices
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");
    
    // Base perfect-match sequence
    let base_spacer = "ATCGATCGATCGATCGATCG";
    
    // Test mismatches at each position (1-based)
    for pos in 1..=20 {
        // Create a sequence with a single mismatch at position pos
        let mut protospacer = base_spacer.to_string();
        let base_at_pos = protospacer.chars().nth(pos-1).unwrap();
        let mismatch_base = match base_at_pos {
            'A' => 'T',
            'T' => 'G',
            'G' => 'C',
            'C' => 'A',
            _ => panic!("Unexpected base: {}", base_at_pos),
        };
        
        // Replace the base at position pos-1 (0-indexed)
        let mut chars: Vec<char> = protospacer.chars().collect();
        chars[pos-1] = mismatch_base;
        protospacer = chars.into_iter().collect();
        
        // Calculate CIGAR string
        let mut cigar = "M".repeat(20);
        let mut cigar_chars: Vec<char> = cigar.chars().collect();
        cigar_chars[pos-1] = 'X';
        cigar = cigar_chars.into_iter().collect();
        
        // Calculate CFD score
        if let Some(score) = cfd_score::get_cfd_score(
            base_spacer.as_bytes(), 
            protospacer.as_bytes(), 
            &cigar, 
            "GG"
        ) {
            println!("Position {} mismatch: CFD score = {:.6}", pos, score);
        } else {
            println!("CFD score calculation failed for position {}", pos);
        }
    }
}

/// Test with different PAM sequences
#[test]
fn test_pam_effects() {
    // Initialize score matrices
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");
    
    // Base perfect-match sequence
    let spacer = "ATCGATCGATCGATCGATCG";
    let protospacer = "ATCGATCGATCGATCGATCG";
    let cigar = "M".repeat(20);
    
    // Test different PAM sequences
    let pams = vec![
        ("GG", 1.0),               // Canonical PAM - highest score
        ("AG", 0.25925925899999996), // Non-canonical but functional
        ("TG", 0.038961038999999996), // Non-canonical but somewhat functional
        ("CG", 0.107142857),       // Non-canonical
        ("AT", 0.0),               // Non-functional
    ];
    
    for (pam, expected_score) in pams {
        if let Some(score) = cfd_score::get_cfd_score(
            spacer.as_bytes(), 
            protospacer.as_bytes(), 
            &cigar, 
            pam
        ) {
            println!("PAM {}: CFD score = {:.6} (expected {:.6})", 
                     pam, score, expected_score);
            
            // Allow some floating point tolerance
            let tolerance = 0.0001;
            assert!((score - expected_score).abs() < tolerance,
                    "PAM {} test failed: got {:.6} but expected {:.6}",
                    pam, score, expected_score);
        } else {
            println!("CFD score calculation failed for PAM {}", pam);
        }
    }
}

/// Test with bulges (insertions/deletions)
#[test]
fn test_bulge_effects() {
    // Initialize score matrices
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");
    
    // Base perfect-match sequence
    let base_spacer = "ATCGATCGATCGATCGATCG";
    
    // Test bulges at different positions
    for pos in 1..=20 {
        // Create a sequence with a deletion at position pos
        let mut del_proto = base_spacer.to_string();
        let mut del_chars: Vec<char> = del_proto.chars().collect();
        del_chars.remove(pos-1);
        del_chars.push('-'); // Add placeholder to keep length
        del_proto = del_chars.into_iter().collect();
        
        // Create CIGAR string for deletion
        let mut del_cigar = String::with_capacity(20);
        for i in 0..20 {
            if i == pos-1 {
                del_cigar.push('D'); // Deletion
            } else {
                del_cigar.push('M'); // Match
            }
        }
        
        // Create a sequence with an insertion at position pos
        let mut ins_proto = base_spacer.to_string();
        let mut ins_chars: Vec<char> = ins_proto.chars().collect();
        ins_chars.insert(pos-1, '-');
        ins_chars.pop(); // Remove last char to keep length
        ins_proto = ins_chars.into_iter().collect();
        
        // Create CIGAR string for insertion
        let mut ins_cigar = String::with_capacity(20);
        for i in 0..20 {
            if i == pos-1 {
                ins_cigar.push('I'); // Insertion
            } else {
                ins_cigar.push('M'); // Match
            }
        }
        
        // Calculate CFD score for deletion
        let del_score = cfd_score::get_cfd_score(
            base_spacer.as_bytes(), 
            del_proto.as_bytes(), 
            &del_cigar, 
            "GG"
        );
        
        // Calculate CFD score for insertion
        let ins_score = cfd_score::get_cfd_score(
            base_spacer.as_bytes(), 
            ins_proto.as_bytes(), 
            &ins_cigar, 
            "GG"
        );
        
        // Print results
        println!("Position {} bulge:", pos);
        println!("  Deletion: CFD score = {:?}", del_score);
        println!("  Insertion: CFD score = {:?}", ins_score);
    }
}
