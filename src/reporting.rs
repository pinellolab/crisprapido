use std::collections::HashMap;
use std::fmt::Write;

use crate::cfd_score;
use crate::hit::Hit;
use crate::reverse_complement;
use crate::verification::cigar_operations;

fn convert_to_minimap2_cigar(cigar: &str) -> String {
    let mut result = String::new();
    let mut count = 0;
    let mut current_op = None;

    for c in cigar.chars() {
        let op = match c {
            'M' => '=',
            'X' | 'I' | 'D' => c,
            _ => continue,
        };

        if Some(op) == current_op {
            count += 1;
        } else {
            if count > 0 {
                write!(result, "{}{}", count, current_op.unwrap()).unwrap();
            }
            current_op = Some(op);
            count = 1;
        }
    }

    if count > 0 && current_op.is_some() {
        write!(result, "{}{}", count, current_op.unwrap()).unwrap();
    }

    result
}

pub(crate) fn report_filtered_hits(hits: Vec<Hit>, pam: &str) {
    print!("{}", format_filtered_hits(hits, pam));
}

fn pam_matches_requested(hit: &Hit, requested_pam: &str) -> bool {
    let Some(pam) = hit.pam_seq.as_deref() else {
        return false;
    };
    pam.eq_ignore_ascii_case(requested_pam)
        || (hit.strand == '-'
            && pam.eq_ignore_ascii_case(&String::from_utf8_lossy(&reverse_complement(
                requested_pam.as_bytes(),
            ))))
}

fn preferred_overlap_hit(candidate: &Hit, current: &Hit, requested_pam: &str) -> bool {
    let candidate_requested_pam = pam_matches_requested(candidate, requested_pam);
    let current_requested_pam = pam_matches_requested(current, requested_pam);

    if candidate_requested_pam != current_requested_pam {
        return candidate_requested_pam;
    }

    if candidate_requested_pam
        && candidate.pos == current.pos
        && candidate.end_pos() != current.end_pos()
    {
        return candidate.end_pos() < current.end_pos();
    }

    candidate.quality_score() > current.quality_score()
}

fn format_filtered_hits(hits: Vec<Hit>, pam: &str) -> String {
    let mut hits_by_group: HashMap<(String, char), Vec<Hit>> = HashMap::new();
    for hit in hits {
        hits_by_group
            .entry((hit.ref_id.clone(), hit.strand))
            .or_insert_with(Vec::new)
            .push(hit);
    }

    let mut groups: Vec<_> = hits_by_group.into_iter().collect();
    groups.sort_by(
        |((left_ref, left_strand), _), ((right_ref, right_strand), _)| {
            left_ref
                .cmp(right_ref)
                .then_with(|| left_strand.cmp(right_strand))
        },
    );

    let mut output = String::new();
    for (_, mut group_hits) in groups {
        group_hits.sort_by_key(|hit| hit.pos);

        let _filtered_hits: Vec<Hit> = Vec::new();
        let mut i = 0;
        while i < group_hits.len() {
            let mut best_idx = i;
            let mut j = i + 1;

            while j < group_hits.len() && group_hits[j].pos < group_hits[i].end_pos() {
                if group_hits[j].overlaps_with(&group_hits[i]) {
                    if preferred_overlap_hit(&group_hits[j], &group_hits[best_idx], pam) {
                        best_idx = j;
                    }
                }
                j += 1;
            }

            let best_hit = &group_hits[best_idx];
            output.push_str(&format_hit(
                &best_hit.ref_id,
                best_hit.pos,
                best_hit.guide.len(),
                best_hit.strand,
                best_hit.score,
                &best_hit.cigar,
                &best_hit.guide,
                best_hit.target_len,
                best_hit.max_mismatches,
                best_hit.max_bulges,
                best_hit.max_bulge_size,
                &best_hit.target_seq,
                best_hit.pam_seq.as_deref().unwrap_or(""),
            ));
            output.push('\n');

            i = j;
        }
    }

    output
}

pub(crate) fn format_hit(
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
) -> String {
    let mut ref_pos = pos;
    let mut ref_consumed = 0;
    let mut query_start = 0;
    let mut query_consumed = 0;

    let cigar_ops = cigar_operations(cigar).unwrap_or_default();

    let leading_dels: usize = cigar_ops
        .iter()
        .take_while(|(_, op)| *op == 'D')
        .map(|(length, _)| *length)
        .sum();
    ref_pos += leading_dels;

    let mut mismatches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut pos = 0;
    for (length, op) in &cigar_ops {
        match op {
            'X' => {
                for _ in 0..*length {
                    if pos < guide.len() && guide[pos] != b'N' {
                        mismatches += 1;
                    }
                    pos += 1;
                }
                ref_consumed += length;
                query_consumed += length;
            }
            'I' => {
                current_gap_size += length;
                if current_gap_size == *length {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                query_consumed += length;
            }
            'D' => {
                current_gap_size += length;
                if current_gap_size == *length {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                ref_consumed += length;
                query_start += length;
            }
            'M' | '=' => {
                current_gap_size = 0;
                ref_consumed += length;
                query_consumed += length;
                pos += length;
            }
            _ => (),
        }
    }

    let mut adjusted_score = 0;
    let mut in_gap = false;
    let mut pos = 0;
    for (length, op) in &cigar_ops {
        match op {
            'X' => {
                for _ in 0..*length {
                    if pos < guide.len() && guide[pos] != b'N' {
                        adjusted_score += 3;
                    }
                    pos += 1;
                }
            }
            'I' | 'D' => {
                if !in_gap {
                    adjusted_score += 5;
                    in_gap = true;
                }
                adjusted_score += length;
                if *op == 'I' {
                    pos += length;
                }
            }
            'M' | '=' => {
                in_gap = false;
                pos += length;
            }
            _ => (),
        }
    }

    let matches: usize = cigar_ops
        .iter()
        .filter_map(|(length, op)| matches!(op, 'M' | '=').then_some(*length))
        .sum();

    let block_len: usize = cigar_ops.iter().map(|(length, _)| *length).sum();
    let guide_len = guide.len();

    macro_rules! debug {
        ($($arg:tt)*) => {
            #[cfg(feature = "debug")]
            eprintln!($($arg)*);
        };
    }

    debug!("Window scan debug:");
    debug!("  CIGAR: {}", cigar);
    debug!("  N-adjusted mismatches: {} (max: 4)", mismatches);
    debug!("  Gaps: {} (max: 1)", gaps);
    debug!("  Max gap size: {} (max: 2)", max_gap_size);
    debug!("  Guide sequence: {}", String::from_utf8_lossy(guide));
    debug!("  Passes filters: true");
    debug!("");

    let cfd_score = if !target_seq.is_empty() && pam.len() == 2 {
        cfd_score::get_cfd_score(guide, target_seq, cigar, pam)
    } else {
        None
    };

    let cfd_tag = if let Some(score) = cfd_score {
        format!("\tcf:f:{:.4}", score)
    } else {
        String::new()
    };

    format!("Guide\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tas:i:{}\tnm:i:{}\tng:i:{}\tbs:i:{}\tcg:Z:{}{}",
        guide_len,
        query_start,
        query_start + query_consumed,
        strand,
        ref_id,
        target_len,
        ref_pos,
        ref_pos + ref_consumed,
        matches,
        block_len,
        adjusted_score,
        mismatches,
        gaps,
        max_gap_size,
        convert_to_minimap2_cigar(cigar),
        cfd_tag
    )
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;

    #[test]
    fn test_formatted_output_contains_cigar_and_cfd_tags() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let line = format_hit(
            "chr1",
            0,
            guide.len(),
            '+',
            0,
            "MMMMMMMMMMMMMMMMMMMM",
            guide,
            guide.len(),
            4,
            1,
            2,
            guide,
            "GG",
        );

        assert!(line.contains("\tcg:Z:20="));
        assert!(line.contains("\tcf:f:"));
    }

    fn test_hit(ref_id: &str, strand: char, pos: usize) -> Hit {
        let guide = Arc::new(b"GAGTCCGAGCAGAAGAAGAA".to_vec());
        Hit {
            ref_id: ref_id.to_string(),
            pos,
            strand,
            score: 0,
            cigar: "MMMMMMMMMMMMMMMMMMMM".to_string(),
            guide: Arc::clone(&guide),
            target_len: 100,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: guide.as_ref().clone(),
            pam_seq: Some("GG".to_string()),
        }
    }

    #[test]
    fn test_filtered_output_deterministic_for_different_insertion_orders() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");
        let hits_a = vec![
            test_hit("chr2", '+', 20),
            test_hit("chr1", '-', 10),
            test_hit("chr1", '+', 30),
            test_hit("chr1", '+', 5),
        ];
        let hits_b = vec![
            test_hit("chr1", '+', 5),
            test_hit("chr2", '+', 20),
            test_hit("chr1", '+', 30),
            test_hit("chr1", '-', 10),
        ];

        assert_eq!(
            format_filtered_hits(hits_a, "GG"),
            format_filtered_hits(hits_b, "GG")
        );
    }

    #[test]
    fn test_requested_pam_hit_wins_overlap_over_nonrequested_pam() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");
        let guide = Arc::new(b"GAGTCCGAGCAGAAGAAGAA".to_vec());
        let nonrequested = Hit {
            ref_id: "chr1".to_string(),
            pos: 10,
            strand: '+',
            score: 0,
            cigar: "MMMMMMMMMMMMMMMMMMMM".to_string(),
            guide: Arc::clone(&guide),
            target_len: 100,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: guide.as_ref().clone(),
            pam_seq: Some("TG".to_string()),
        };
        let requested = Hit {
            ref_id: "chr1".to_string(),
            pos: 12,
            strand: '+',
            score: 7,
            cigar: "IIMMMMMMMMMMMMMMMMMM".to_string(),
            guide: Arc::clone(&guide),
            target_len: 100,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: b"--GTCCGAGCAGAAGAAGAA".to_vec(),
            pam_seq: Some("GG".to_string()),
        };

        let output = format_filtered_hits(vec![nonrequested, requested], "GG");

        assert_eq!(output.lines().count(), 1);
        assert!(output.contains("	12	30	"));
        assert!(output.contains("	cg:Z:2I18="));
    }

    #[test]
    fn test_same_start_requested_pam_prefers_shorter_reference_span() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");
        let guide = Arc::new(b"GAGTCCGAGCAGAAGAAGAA".to_vec());
        let longer = Hit {
            ref_id: "chr1".to_string(),
            pos: 10,
            strand: '+',
            score: 6,
            cigar: "MMMMDMMMMMMMMMMMMMMMM".to_string(),
            guide: Arc::clone(&guide),
            target_len: 100,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: guide.as_ref().clone(),
            pam_seq: Some("GG".to_string()),
        };
        let shorter = Hit {
            ref_id: "chr1".to_string(),
            pos: 10,
            strand: '+',
            score: 7,
            cigar: "IIMMMMMMMMMMMMMMMMMM".to_string(),
            guide: Arc::clone(&guide),
            target_len: 100,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: b"--GTCCGAGCAGAAGAAGAA".to_vec(),
            pam_seq: Some("GG".to_string()),
        };

        let output = format_filtered_hits(vec![longer, shorter], "GG");

        assert_eq!(output.lines().count(), 1);
        assert!(output.contains("	10	28	"));
        assert!(output.contains("	cg:Z:2I18="));
    }
}
