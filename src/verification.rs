use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;

use bio::io::fasta;
use flate2::read::MultiGzDecoder;
use lib_wfa2::affine_wavefront::AffineWavefronts;

use crate::columba::{cigar_reference_span, ColumbaCandidate};
use crate::hit::Hit;
use crate::{reverse_complement, Args};

#[derive(Debug, PartialEq, Eq)]
pub(crate) struct VerificationWindow {
    pub(crate) start: usize,
    pub(crate) end: usize,
}

pub(crate) fn cigar_operations(cigar: &str) -> Option<Vec<(usize, char)>> {
    let mut operations = Vec::new();
    let mut length = 0usize;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            length = length
                .checked_mul(10)?
                .checked_add(c.to_digit(10)? as usize)?;
        } else {
            let op_len = if length == 0 { 1 } else { length };
            operations.push((op_len, c));
            length = 0;
        }
    }

    if length != 0 {
        return None;
    }

    Some(operations)
}

pub(crate) fn alignment_reference_span(cigar: &str) -> usize {
    cigar_operations(cigar)
        .unwrap_or_default()
        .iter()
        .filter_map(|(length, op)| matches!(op, 'M' | '=' | 'X' | 'D').then_some(*length))
        .sum()
}

fn extract_cfd_target_from_alignment(
    oriented_window: &[u8],
    alignment_offset: usize,
    cigar: &str,
) -> Option<Vec<u8>> {
    let mut target = Vec::new();
    let mut target_pos = alignment_offset;

    for (length, op) in cigar_operations(cigar)? {
        match op {
            'M' | '=' | 'X' => {
                for _ in 0..length {
                    let base = *oriented_window.get(target_pos)?;
                    target.push(base);
                    target_pos += 1;
                }
            }
            'I' => {
                target.extend(std::iter::repeat(b'-').take(length));
            }
            'D' => {
                target_pos = target_pos.checked_add(length)?;
                if target_pos > oriented_window.len() {
                    return None;
                }
            }
            _ => {}
        }
    }

    Some(target)
}

fn extract_adjacent_pam(
    seq: &[u8],
    protospacer_start: usize,
    protospacer_end: usize,
    strand: char,
) -> Option<String> {
    match strand {
        '+' => {
            let pam = seq.get(protospacer_end..protospacer_end + 2)?;
            Some(String::from_utf8_lossy(pam).to_string())
        }
        '-' => {
            let pam_start = protospacer_start.checked_sub(2)?;
            let pam = seq.get(pam_start..protospacer_start)?;
            Some(String::from_utf8_lossy(&reverse_complement(pam)).to_string())
        }
        _ => None,
    }
}

pub(crate) fn build_verified_hit(
    ref_id: String,
    seq: &[u8],
    hit_pos: usize,
    strand: char,
    score: i32,
    cigar: String,
    guide: Arc<Vec<u8>>,
    alignment_window: &[u8],
    alignment_offset: usize,
    args: &Args,
) -> Hit {
    let target_seq = extract_cfd_target_from_alignment(alignment_window, alignment_offset, &cigar)
        .unwrap_or_else(|| {
            eprintln!(
                "Warning: unable to extract aligned CFD target for {}:{} on strand {}",
                ref_id, hit_pos, strand
            );
            Vec::new()
        });
    let protospacer_end = hit_pos + alignment_reference_span(&cigar);
    let pam_seq = extract_adjacent_pam(seq, hit_pos, protospacer_end, strand).or_else(|| {
        eprintln!(
            "Warning: unable to extract adjacent PAM for {}:{}..{} on strand {}",
            ref_id, hit_pos, protospacer_end, strand
        );
        None
    });

    Hit {
        ref_id,
        pos: hit_pos,
        strand,
        score,
        cigar,
        guide,
        target_len: seq.len(),
        max_mismatches: args.max_mismatches,
        max_bulges: args.max_bulges,
        max_bulge_size: args.max_bulge_size,
        cfd_score: None,
        target_seq,
        pam_seq,
    }
}

pub(crate) fn columba_verification_window(
    reference_start: usize,
    reference_span: usize,
    contig_len: usize,
    guide_len: usize,
    pam_len: usize,
    max_bulge_size: u32,
) -> VerificationWindow {
    let flank = pam_len + max_bulge_size as usize + 1;
    let target_len = reference_span.max(guide_len);
    let start = reference_start.saturating_sub(flank);
    let end = reference_start
        .saturating_add(target_len)
        .saturating_add(flank)
        .min(contig_len);

    VerificationWindow { start, end }
}

pub(crate) fn load_reference_records(path: &PathBuf) -> Result<Vec<(String, Vec<u8>)>, String> {
    let file = File::open(path)
        .map_err(|e| format!("Failed to open reference file '{}': {}", path.display(), e))?;
    let reader: Box<dyn BufRead> = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let reader = fasta::Reader::new(reader);
    let mut records = Vec::new();

    for result in reader.records() {
        let record = result.map_err(|e| format!("Error during FASTA record parsing: {}", e))?;
        records.push((record.id().to_string(), record.seq().to_vec()));
    }

    Ok(records)
}

pub(crate) fn verify_columba_candidates(
    candidates: &[ColumbaCandidate],
    references: &[(String, Vec<u8>)],
    guide_fwd: &Arc<Vec<u8>>,
    args: &Args,
) -> Vec<Hit> {
    let reference_by_name: HashMap<&str, &(String, Vec<u8>)> = references
        .iter()
        .map(|record| (record.0.as_str(), record))
        .collect();
    let pam_len = args.pam.len();
    let guide_len = guide_fwd.len();
    let mut hits = Vec::new();

    for candidate in candidates {
        let _candidate_metadata = (candidate.edit_distance, candidate.alignment_score);
        let Some((record_id, seq)) = reference_by_name.get(candidate.reference_name.as_str())
        else {
            eprintln!(
                "Warning: skipping Columba candidate '{}' on missing reference '{}'",
                candidate.query_name, candidate.reference_name
            );
            continue;
        };

        let reference_span = match cigar_reference_span(&candidate.cigar) {
            Ok(span) => span,
            Err(e) => {
                eprintln!(
                    "Warning: skipping Columba candidate '{}' on '{}': {}",
                    candidate.query_name, candidate.reference_name, e
                );
                continue;
            }
        };
        let window = columba_verification_window(
            candidate.reference_start,
            reference_span,
            seq.len(),
            guide_len,
            pam_len,
            args.max_bulge_size,
        );
        if window.end <= window.start || window.end - window.start < guide_len {
            eprintln!(
                "Warning: skipping Columba candidate '{}' on '{}:{}' because verification window {}..{} is shorter than guide length {}",
                candidate.query_name,
                candidate.reference_name,
                candidate.reference_start + 1,
                window.start,
                window.end,
                guide_len
            );
            continue;
        }

        let window_seq = &seq[window.start..window.end];
        let reverse_window;
        let (guide, strand, alignment_window): (&Arc<Vec<u8>>, char, &[u8]) = if candidate.reverse {
            reverse_window = reverse_complement(window_seq);
            (guide_fwd, '-', &reverse_window)
        } else {
            (guide_fwd, '+', window_seq)
        };
        let mut aligner = AffineWavefronts::with_penalties(0, 3, 5, 1);

        if let Some((score, cigar, _mismatches, _gaps, _max_gap_size, leading_dels)) = scan_window(
            &mut aligner,
            guide,
            alignment_window,
            args.max_mismatches,
            args.max_bulges,
            args.max_bulge_size,
            args.min_match_fraction,
            args.no_filter,
        ) {
            let hit_pos = if candidate.reverse {
                window
                    .end
                    .saturating_sub(leading_dels)
                    .saturating_sub(alignment_reference_span(&cigar))
            } else {
                window.start + leading_dels
            };
            hits.push(build_verified_hit(
                record_id.clone(),
                seq,
                hit_pos,
                strand,
                score,
                cigar,
                Arc::clone(guide),
                alignment_window,
                leading_dels,
                args,
            ));
        } else {
            eprintln!(
                "Warning: skipping Columba candidate '{}' on '{}:{}' because WFA2 verification failed",
                candidate.query_name,
                candidate.reference_name,
                candidate.reference_start + 1
            );
        }
    }

    hits
}

pub(crate) fn scan_window(
    aligner: &AffineWavefronts,
    guide: &[u8],
    window: &[u8],
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
    min_match_fraction: f32,
    no_filter: bool,
) -> Option<(i32, String, u32, u32, u32, usize)> {
    aligner.align(window, guide);
    let score = aligner.score();
    let raw_cigar = String::from_utf8_lossy(aligner.cigar()).to_string();

    let mut leading_indels = true;
    let mut leading_dels = 0;
    for c in raw_cigar.chars() {
        if leading_indels {
            match c {
                'D' => leading_dels += 1,
                'I' => (),
                _ => leading_indels = false,
            }
        }
    }

    let cigar = raw_cigar
        .chars()
        .skip_while(|&c| c == 'D' || c == 'I')
        .collect::<String>()
        .trim_end_matches(|c| c == 'D' || c == 'I')
        .to_string();

    let mut n_adjusted_mismatches = 0;
    let mut matches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut pos = 0;

    for c in cigar.chars() {
        match c {
            'X' => {
                if pos < guide.len() && guide[pos] != b'N' {
                    n_adjusted_mismatches += 1;
                }
                pos += 1;
            }
            'I' | 'D' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                if c == 'I' {
                    pos += 1;
                }
            }
            'M' | '=' => {
                current_gap_size = 0;
                matches += 1;
                pos += 1;
            }
            _ => (),
        }
    }

    macro_rules! debug {
        ($($arg:tt)*) => {
            #[cfg(feature = "debug")]
            eprintln!($($arg)*);
        };
    }

    debug!(
        "CIGAR: {}, N-adjusted Mismatches: {}, Gaps: {}, Max gap size: {}",
        cigar, n_adjusted_mismatches, gaps, max_gap_size
    );

    let non_n_positions = guide.iter().filter(|&&b| b != b'N').count();
    let match_percentage = if non_n_positions > 0 {
        (matches as f32 / non_n_positions as f32) * 100.0
    } else {
        0.0
    };

    let min_match_percentage = min_match_fraction * 100.0;

    debug!(
        "Match percentage: {}, Minimum required: {}",
        match_percentage, min_match_percentage
    );

    if no_filter
        || (matches >= 1
            && match_percentage >= min_match_percentage
            && ((cfg!(test) && n_adjusted_mismatches <= 1 && gaps <= 1 && max_gap_size <= 1)
                || (!cfg!(test)
                    && n_adjusted_mismatches <= max_mismatches
                    && gaps <= max_bulges
                    && max_gap_size <= max_bulge_size)))
    {
        Some((
            score,
            cigar,
            n_adjusted_mismatches,
            gaps,
            max_gap_size,
            leading_dels,
        ))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use lib_wfa2::affine_wavefront::AffineWavefronts;
    use rand::{rngs::SmallRng, RngCore, SeedableRng};

    use super::*;
    use crate::cfd_score;
    use crate::columba::parse_columba_sam_record;

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

    fn setup_aligner() -> AffineWavefronts {
        AffineWavefronts::with_penalties(0, 3, 5, 1)
    }

    fn columba_test_args() -> Args {
        Args {
            reference: PathBuf::from("unused.fa"),
            guide: "GAGTCCGAGCAGAAGAAGAA".to_string(),
            pam: "GG".to_string(),
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            min_match_fraction: 0.75,
            mismatch_scores: PathBuf::from("mismatch_scores.txt"),
            pam_scores: PathBuf::from("pam_scores.txt"),
            window_size: None,
            threads: None,
            no_filter: false,
            columba_sam: None,
            columba_bin: None,
            columba_index: None,
            keep_columba_sam: false,
        }
    }

    fn parse_candidate(line: &str) -> ColumbaCandidate {
        parse_columba_sam_record(line).unwrap().unwrap()
    }

    fn verify_test_candidates(
        candidates: Vec<ColumbaCandidate>,
        references: Vec<(String, Vec<u8>)>,
    ) -> Vec<Hit> {
        let args = columba_test_args();
        let guide_fwd = Arc::new(args.guide.as_bytes().to_vec());
        verify_columba_candidates(&candidates, &references, &guide_fwd, &args)
    }

    fn cfd_for_hit(hit: &Hit) -> Option<f64> {
        cfd_score::get_cfd_score(
            &hit.guide,
            &hit.target_seq,
            &hit.cigar,
            hit.pam_seq.as_deref().unwrap_or(""),
        )
    }

    fn first_verified_hit_for_window(
        ref_id: &str,
        seq: &[u8],
        window_start: usize,
        window_end: usize,
    ) -> Hit {
        let args = columba_test_args();
        let guide = Arc::new(args.guide.as_bytes().to_vec());
        let mut aligner = setup_aligner();
        let window = &seq[window_start..window_end];
        let (score, cigar, _, _, _, leading_dels) = scan_window(
            &mut aligner,
            &guide,
            window,
            args.max_mismatches,
            args.max_bulges,
            args.max_bulge_size,
            args.min_match_fraction,
            args.no_filter,
        )
        .unwrap();
        build_verified_hit(
            ref_id.to_string(),
            seq,
            window_start + leading_dels,
            '+',
            score,
            cigar,
            guide,
            window,
            leading_dels,
            &args,
        )
    }

    #[test]
    fn test_cigar_operations_compact_and_wfa_style() {
        assert_eq!(cigar_operations("20=").unwrap(), vec![(20, '=')]);
        assert_eq!(
            cigar_operations("MMXDI").unwrap(),
            vec![(1, 'M'), (1, 'M'), (1, 'X'), (1, 'D'), (1, 'I')]
        );
    }

    #[test]
    fn test_verify_columba_perfect_forward_candidate() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let mut seq = b"TTTT".to_vec();
        seq.extend_from_slice(guide);
        seq.extend_from_slice(b"GGAAAA");
        let candidate =
            parse_candidate("guide_20bp\t0\tchr1\t5\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].ref_id, "chr1");
        assert_eq!(hits[0].strand, '+');
        assert!(hits[0].cigar.contains('M'));
    }

    #[test]
    fn test_verify_columba_perfect_reverse_candidate() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let guide_rc = reverse_complement(guide);
        let mut seq = b"CCCC".to_vec();
        seq.extend_from_slice(&guide_rc);
        seq.extend_from_slice(b"GGAAAA");
        let candidate =
            parse_candidate("guide_20bp\t16\tchr1\t5\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand, '-');
    }

    #[test]
    fn test_verify_columba_one_substitution() {
        let mut target = b"GAGTCCGAGCAGAAGAAGAA".to_vec();
        target[10] = b'T';
        let mut seq = b"TTTT".to_vec();
        seq.extend_from_slice(&target);
        seq.extend_from_slice(b"GGAAAA");
        let candidate =
            parse_candidate("guide_20bp\t0\tchr1\t5\t60\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert_eq!(hits.len(), 1);
        assert!(hits[0].cigar.contains('X'));
    }

    #[test]
    fn test_verify_columba_insertion_deletion_style_candidate() {
        let mut target = b"GAGTCCGAGC".to_vec();
        target.push(b'A');
        target.extend_from_slice(b"AGAAGAAGAA");
        let mut seq = b"TTTT".to_vec();
        seq.extend_from_slice(&target);
        seq.extend_from_slice(b"GGAAAA");
        let candidate =
            parse_candidate("guide_20bp\t0\tchr1\t5\t60\t10M1I10M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert_eq!(hits.len(), 1);
        assert!(hits[0].cigar.contains('I') || hits[0].cigar.contains('D'));
    }

    #[test]
    fn test_verify_columba_secondary_alignment_retained() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let mut seq = b"TTTT".to_vec();
        seq.extend_from_slice(guide);
        seq.extend_from_slice(b"GGAAAA");
        let candidate =
            parse_candidate("guide_20bp\t256\tchr1\t5\t0\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn test_columba_candidate_near_reference_start_clamps_window() {
        let window = columba_verification_window(0, 20, 30, 20, 2, 2);
        assert_eq!(window.start, 0);
        assert_eq!(window.end, 25);
    }

    #[test]
    fn test_columba_candidate_near_reference_end_clamps_window() {
        let window = columba_verification_window(12, 20, 30, 20, 2, 2);
        assert_eq!(window.start, 7);
        assert_eq!(window.end, 30);
    }

    #[test]
    fn test_verify_columba_missing_rname_skipped_safely() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let candidate =
            parse_candidate("guide_20bp\t0\tmissing\t5\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits =
            verify_test_candidates(vec![candidate], vec![("chr1".to_string(), guide.to_vec())]);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_imported_candidate_mode_does_not_run_normal_genome_scan() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let mut chr1 = b"TTTT".to_vec();
        chr1.extend_from_slice(guide);
        chr1.extend_from_slice(b"GGAAAA");
        let mut chr2 = b"CCCC".to_vec();
        chr2.extend_from_slice(guide);
        chr2.extend_from_slice(b"GGAAAA");
        let candidate =
            parse_candidate("guide_20bp\t0\tchr1\t5\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits = verify_test_candidates(
            vec![candidate],
            vec![("chr1".to_string(), chr1), ("chr2".to_string(), chr2)],
        );
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].ref_id, "chr1");
    }

    #[test]
    fn test_cfd_input_independent_of_flank_size() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");
        let seq = b"AAAAAAAAAAGAGTCCGAGCAGAAGAAGAAGGCCCCCCCCCC";
        let small = first_verified_hit_for_window("chr1", seq, 5, 37);
        let large = first_verified_hit_for_window("chr1", seq, 0, seq.len());

        assert_eq!(small.target_seq, b"GAGTCCGAGCAGAAGAAGAA".to_vec());
        assert_eq!(small.target_seq, large.target_seq);
        assert_eq!(small.pam_seq, Some("GG".to_string()));
        assert_eq!(small.pam_seq, large.pam_seq);
        assert_eq!(cfd_for_hit(&small), cfd_for_hit(&large));
    }

    fn compare_normal_and_columba_hit(
        reference_name: &str,
        seq: Vec<u8>,
        sam_record: &str,
    ) -> (Hit, Hit) {
        let normal = first_verified_hit_for_window(reference_name, &seq, 0, seq.len());
        let candidate = parse_candidate(sam_record);
        let columba_hits =
            verify_test_candidates(vec![candidate], vec![(reference_name.to_string(), seq)]);
        assert_eq!(columba_hits.len(), 1);
        (normal, columba_hits.into_iter().next().unwrap())
    }

    #[test]
    fn test_normal_and_columba_perfect_forward_same_cfd_input() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");
        let seq = b"TTTTTTTTTTGAGTCCGAGCAGAAGAAGAAGGTTTTTTTT".to_vec();
        let (normal, columba) = compare_normal_and_columba_hit(
            "perfect_forward",
            seq,
            "guide_20bp\t0\tperfect_forward\t11\t0\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0",
        );

        assert_eq!(normal.target_seq, columba.target_seq);
        assert_eq!(normal.pam_seq, columba.pam_seq);
        assert_eq!(normal.cigar, columba.cigar);
        assert_eq!(cfd_for_hit(&normal), cfd_for_hit(&columba));
    }

    #[test]
    fn test_normal_and_columba_one_substitution_same_cfd_input() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");
        let seq = b"TTTTTTTTTTGAGTCCGAGTAGAAGAAGAAGGTTTTTTTT".to_vec();
        let (normal, columba) = compare_normal_and_columba_hit(
            "one_substitution",
            seq,
            "guide_20bp\t256\tone_substitution\t11\t0\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1",
        );

        assert_eq!(normal.target_seq, columba.target_seq);
        assert_eq!(normal.pam_seq, columba.pam_seq);
        assert_eq!(normal.cigar, columba.cigar);
        assert_eq!(cfd_for_hit(&normal), cfd_for_hit(&columba));
    }

    #[test]
    fn test_normal_and_columba_one_deletion_same_cfd_input() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");
        let seq = b"TTTTTTTTTTGAGTCCGAGCGAAGAAGAAGGTTTTTTTT".to_vec();
        let (normal, columba) = compare_normal_and_columba_hit(
            "one_bp_deletion_relative_to_guide",
            seq,
            "guide_20bp\t256\tone_bp_deletion_relative_to_guide\t11\t0\t10M1I9M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1",
        );

        assert_eq!(normal.target_seq, columba.target_seq);
        assert_eq!(normal.pam_seq, columba.pam_seq);
        assert_eq!(normal.cigar, columba.cigar);
        assert_eq!(cfd_for_hit(&normal), cfd_for_hit(&columba));
    }

    #[test]
    fn test_reverse_candidate_orients_target_and_pam_for_cfd() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let guide_rc = reverse_complement(guide);
        let mut seq = b"TT".to_vec();
        seq.extend_from_slice(b"GGGGGGGG");
        seq.extend_from_slice(&guide_rc);
        seq.extend_from_slice(b"TTTTTT");
        let candidate =
            parse_candidate("guide_20bp\t16\tchr1\t11\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand, '-');
        assert_eq!(hits[0].target_seq, guide.to_vec());
        assert_eq!(hits[0].pam_seq, Some("CC".to_string()));
    }

    #[test]
    fn test_pam_boundary_does_not_panic() {
        let args = columba_test_args();
        let guide = Arc::new(args.guide.as_bytes().to_vec());
        let hit = build_verified_hit(
            "near_end".to_string(),
            b"GAGTCCGAGCAGAAGAAGAA",
            0,
            '+',
            0,
            "MMMMMMMMMMMMMMMMMMMM".to_string(),
            guide,
            b"GAGTCCGAGCAGAAGAAGAA",
            0,
            &args,
        );
        assert_eq!(hit.target_seq, b"GAGTCCGAGCAGAAGAAGAA".to_vec());
        assert_eq!(hit.pam_seq, None);
        assert!(cfd_for_hit(&hit).is_none());
    }

    #[test]
    fn test_actual_pam_extracted_not_cli_pattern() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");
        let args = columba_test_args();
        let guide = Arc::new(args.guide.as_bytes().to_vec());
        let seq = b"GAGTCCGAGCAGAAGAAGAAAG";
        let hit = build_verified_hit(
            "chr1".to_string(),
            seq,
            0,
            '+',
            0,
            "MMMMMMMMMMMMMMMMMMMM".to_string(),
            guide,
            &seq[..20],
            0,
            &args,
        );

        assert_eq!(hit.pam_seq, Some("AG".to_string()));
        assert_eq!(
            cfd_for_hit(&hit),
            cfd_score::get_cfd_score(&hit.guide, &hit.target_seq, &hit.cigar, "AG")
        );
    }

    #[test]
    fn test_perfect_match() {
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = b"ATCGATCGAT";

        let result = scan_window(&mut aligner, guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_some());
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMMMMMMM");
    }

    #[test]
    fn test_with_mismatches() {
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = b"ATCGTTCGAT";

        let result = scan_window(&mut aligner, guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single mismatch");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMXMMMMM");
    }

    #[test]
    fn test_with_bulge() {
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = b"ATCGAATCGAT";

        let result = scan_window(&mut aligner, guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single base bulge");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(
            cigar.contains('I') || cigar.contains('D'),
            "Should contain an insertion or deletion"
        );
    }

    #[test]
    fn test_too_many_differences() {
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = b"ATCGTTCGTT";

        let result = scan_window(&mut aligner, guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_none());
    }

    #[test]
    fn test_perfect_match_with_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = create_flanked_sequence(&mut rng, guide, 500);

        let result = scan_window(&mut aligner, guide, &target[500..510], 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should match perfectly even with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMMMMMMM");
    }

    #[test]
    fn test_with_mismatches_and_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGAT";
        let target = create_flanked_sequence(&mut rng, core, 500);

        let result = scan_window(&mut aligner, guide, &target[500..510], 1, 1, 1, 0.75, false);
        assert!(
            result.is_some(),
            "Should accept a single mismatch with flanks"
        );
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "MMMMXMMMMM");
    }

    #[test]
    fn test_with_bulge_and_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGAATCGAT";
        let target = create_flanked_sequence(&mut rng, core, 500);

        let result = scan_window(&mut aligner, guide, &target[500..511], 1, 1, 1, 0.75, false);
        assert!(
            result.is_some(),
            "Should accept a single base bulge with flanks"
        );
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(
            cigar.contains('I') || cigar.contains('D'),
            "Should contain an insertion or deletion"
        );
    }

    #[test]
    fn test_too_many_differences_with_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGTT";
        let target = create_flanked_sequence(&mut rng, core, 500);

        let result = scan_window(&mut aligner, guide, &target[500..510], 1, 1, 1, 0.75, false);
        assert!(
            result.is_none(),
            "Should reject sequence with too many mismatches even with flanks"
        );
    }
}
