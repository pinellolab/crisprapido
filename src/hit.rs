use std::sync::Arc;

#[derive(Clone)]
pub(crate) struct Hit {
    pub(crate) ref_id: String,
    pub(crate) pos: usize,
    pub(crate) strand: char,
    pub(crate) score: i32,
    pub(crate) cigar: String,
    pub(crate) guide: Arc<Vec<u8>>,
    pub(crate) target_len: usize,
    pub(crate) max_mismatches: u32,
    pub(crate) max_bulges: u32,
    pub(crate) max_bulge_size: u32,
    pub(crate) cfd_score: Option<f64>,
    pub(crate) target_seq: Vec<u8>,
    pub(crate) pam_seq: Option<String>,
}

impl Hit {
    pub(crate) fn quality_score(&self) -> i32 {
        let matches = self.cigar.chars().filter(|&c| c == 'M' || c == '=').count();

        let mismatches = self.cigar.chars().filter(|&c| c == 'X').count();

        let gaps = self.cigar.chars().filter(|&c| c == 'I' || c == 'D').count();

        matches as i32 - (mismatches as i32) - (gaps as i32 * 2) - self.score
    }

    pub(crate) fn end_pos(&self) -> usize {
        let mut ref_consumed = 0;
        for c in self.cigar.chars() {
            match c {
                'M' | '=' | 'X' | 'D' => ref_consumed += 1,
                _ => {}
            }
        }
        self.pos + ref_consumed
    }

    pub(crate) fn overlaps_with(&self, other: &Hit) -> bool {
        self.strand == other.strand
            && self.ref_id == other.ref_id
            && self.pos < other.end_pos()
            && other.pos < self.end_pos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hit_quality_scoring_and_filtering() {
        let guide_seq = Arc::new(b"ATCGATCGAT".to_vec());

        let perfect_hit = Hit {
            ref_id: "chr1".to_string(),
            pos: 100,
            strand: '+',
            score: 0,
            cigar: "MMMMMMMMMM".to_string(),
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
            pam_seq: None,
        };

        let mismatch_hit = Hit {
            ref_id: "chr1".to_string(),
            pos: 105,
            strand: '+',
            score: 3,
            cigar: "MMMMXMMMMM".to_string(),
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
            pam_seq: None,
        };

        let bulge_hit = Hit {
            ref_id: "chr1".to_string(),
            pos: 110,
            strand: '+',
            score: 6,
            cigar: "MMMDMMMMM".to_string(),
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
            pam_seq: None,
        };

        assert!(
            perfect_hit.overlaps_with(&mismatch_hit),
            "Hits should overlap"
        );
        assert!(
            mismatch_hit.overlaps_with(&perfect_hit),
            "Overlap should be symmetric"
        );
        assert!(
            !perfect_hit.overlaps_with(&bulge_hit),
            "These hits shouldn't overlap"
        );

        assert!(
            perfect_hit.quality_score() > mismatch_hit.quality_score(),
            "Perfect match should have higher quality than mismatch"
        );
        assert!(
            mismatch_hit.quality_score() > bulge_hit.quality_score(),
            "Mismatch should have higher quality than bulge"
        );

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
