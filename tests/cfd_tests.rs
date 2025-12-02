//! Integration tests for CFD score calculation
//!
//! These tests verify basic CFD functionality.
//! For comprehensive Python reference comparison, see tests/cfd_python_comparison.rs

extern crate crisprapido;

use crisprapido::cfd_score;

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
        let base_at_pos = protospacer.chars().nth(pos - 1).unwrap();
        let mismatch_base = match base_at_pos {
            'A' => 'T',
            'T' => 'G',
            'G' => 'C',
            'C' => 'A',
            _ => panic!("Unexpected base: {}", base_at_pos),
        };

        // Replace the base at position pos-1 (0-indexed)
        let mut chars: Vec<char> = protospacer.chars().collect();
        chars[pos - 1] = mismatch_base;
        protospacer = chars.into_iter().collect();

        // Calculate CFD score using direct method
        let score = cfd_score::calculate_cfd(base_spacer, &protospacer, "GG")
            .expect("CFD calculation failed");

        // Score should be in valid range
        assert!(
            score >= 0.0 && score <= 1.0,
            "Position {} mismatch: score {} out of range [0,1]",
            pos,
            score
        );
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

    // Test different PAM sequences (expected values from Python reference)
    let pams = vec![
        ("GG", 1.0),                   // Canonical PAM - highest score
        ("AG", 0.25925925899999996),   // Non-canonical but functional
        ("TG", 0.038961038999999996),  // Non-canonical
        ("CG", 0.107142857),           // Non-canonical
        ("AT", 0.0),                   // Non-functional
        ("NGG", 1.0),                  // 3nt PAM should collapse to GG
        ("NAG", 0.25925925899999996),  // 3nt variant matching AG penalty
    ];

    for (pam, expected_score) in pams {
        let score = cfd_score::calculate_cfd(spacer, spacer, pam)
            .expect("CFD calculation failed");

        // Allow some floating point tolerance
        let tolerance = 1e-6;
        assert!(
            (score - expected_score).abs() < tolerance,
            "PAM {} test failed: got {:.9} but expected {:.9}",
            pam,
            score,
            expected_score
        );
    }
}

#[test]
fn test_terminal_ns_are_neutral() {
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");

    let spacer = "CACACAGGTCCAAAGGCGGG";
    let target = "NNNNNNNNNNNNNNNNNNGG";

    let score = cfd_score::calculate_cfd(spacer, target, "GG").unwrap();
    assert_eq!(score, 1.0, "Ambiguous Ns at distal end should not penalize score");
}

/// Test CIGAR-based calculation matches direct calculation
#[test]
fn test_cigar_matches_direct() {
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");

    let spacer = "ATCGATCGATCGATCGATCG";
    let protospacer = "ATCGATCGATCGATCGATCG";

    // Direct calculation
    let direct_score = cfd_score::calculate_cfd(spacer, protospacer, "GG").unwrap();

    // CIGAR-based calculation (20 matches)
    let cigar_score =
        cfd_score::get_cfd_score(spacer.as_bytes(), protospacer.as_bytes(), "20=", "GG").unwrap();

    assert!(
        (direct_score - cigar_score).abs() < 1e-9,
        "Direct ({}) and CIGAR ({}) calculations should match",
        direct_score,
        cigar_score
    );
}

/// Test that gap handling works correctly via CIGAR
#[test]
fn test_cigar_gap_handling() {
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");

    let guide = b"ATCGATCGATCGATCGATCG";
    let target = b"ATCGATCGATCGATCGATCG";

    // Test with insertion at position 10 (gap in target)
    // Guide:  ATCGATCGATCGATCGATCG
    // Target: ATCGATCGAT-GATCGATCG (with gap)
    let score = cfd_score::get_cfd_score(guide, &target[..19], "10=1I9=", "GG");
    assert!(score.is_some(), "CIGAR with insertion should work");

    // Position 11 insertion should have a penalty
    if let Some(s) = score {
        // The score should be reduced due to the gap
        assert!(s <= 1.0, "Gap should not increase score above 1.0");
    }
}
