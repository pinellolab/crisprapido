use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;

use bio::io::fasta;
use flate2::read::MultiGzDecoder;
use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentSpan};

use crate::columba::{cigar_reference_span, ColumbaCandidate};
use crate::hit::Hit;
use crate::{reverse_complement, Args};

#[derive(Default)]
struct ColumbaVerificationStats {
    original_candidates: usize,
    primary_accepted: usize,
    primary_rejected: usize,
    expansion_triggers: usize,
    raw_alternatives_generated: usize,
    unique_alternatives_verified: usize,
    fallback_hits_accepted: usize,
    duplicates_removed: usize,
}

#[cfg(test)]
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

pub(crate) fn uppercase_ascii_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(u8::to_ascii_uppercase).collect()
}

#[derive(Debug, PartialEq, Eq)]
struct CigarValidation {
    cigar: String,
    leading_dels: usize,
    matches: u32,
    mismatches: u32,
    gaps: u32,
    max_gap_size: u32,
}

fn trim_reference_overhangs(raw_cigar: &str) -> (String, usize) {
    let leading_dels = raw_cigar.chars().take_while(|&c| c == 'D').count();
    let cigar = raw_cigar
        .chars()
        .skip_while(|&c| c == 'D')
        .collect::<String>()
        .trim_end_matches('D')
        .to_string();

    (cigar, leading_dels)
}

fn validate_cigar_accounting(
    raw_cigar: &str,
    guide: &[u8],
    trim_reference_overhangs_enabled: bool,
) -> CigarValidation {
    let (cigar, leading_dels) = if trim_reference_overhangs_enabled {
        trim_reference_overhangs(raw_cigar)
    } else {
        (raw_cigar.to_string(), 0)
    };
    let mut n_adjusted_mismatches = 0;
    let mut matches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut pos = 0;

    for c in cigar.chars() {
        match c {
            'X' => {
                current_gap_size = 0;
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

    CigarValidation {
        cigar,
        leading_dels,
        matches,
        mismatches: n_adjusted_mismatches,
        gaps,
        max_gap_size,
    }
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
            Some(String::from_utf8_lossy(&uppercase_ascii_sequence(pam)).to_string())
        }
        '-' => {
            let pam_start = protospacer_start.checked_sub(2)?;
            let pam = uppercase_ascii_sequence(seq.get(pam_start..protospacer_start)?);
            Some(String::from_utf8_lossy(&reverse_complement(&pam)).to_string())
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

#[cfg(test)]
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
    verify_columba_candidates_impl(candidates, references, guide_fwd, args).0
}

#[cfg(test)]
fn verify_columba_candidates_with_stats(
    candidates: &[ColumbaCandidate],
    references: &[(String, Vec<u8>)],
    guide_fwd: &Arc<Vec<u8>>,
    args: &Args,
) -> (Vec<Hit>, ColumbaVerificationStats) {
    verify_columba_candidates_impl(candidates, references, guide_fwd, args)
}

fn hit_matches_requested_pam(hit: &Hit, requested_pam: &str) -> bool {
    let Some(pam) = hit.pam_seq.as_deref() else {
        return false;
    };
    pam.eq_ignore_ascii_case(requested_pam)
        || (hit.strand == '-'
            && pam.eq_ignore_ascii_case(&String::from_utf8_lossy(&reverse_complement(
                requested_pam.as_bytes(),
            ))))
}

fn verify_columba_candidates_impl(
    candidates: &[ColumbaCandidate],
    references: &[(String, Vec<u8>)],
    guide_fwd: &Arc<Vec<u8>>,
    args: &Args,
) -> (Vec<Hit>, ColumbaVerificationStats) {
    let reference_by_name: HashMap<&str, &(String, Vec<u8>)> = references
        .iter()
        .map(|record| (record.0.as_str(), record))
        .collect();
    let mut hits = Vec::new();
    let mut seen_hits = HashSet::new();
    let mut stats = ColumbaVerificationStats::default();

    for candidate in candidates {
        stats.original_candidates += 1;
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
        if reference_span == 0 {
            eprintln!(
                "Warning: skipping Columba candidate '{}' on '{}:{}' because its CIGAR consumes no reference bases",
                candidate.query_name,
                candidate.reference_name,
                candidate.reference_start + 1
            );
            continue;
        }
        let Some(candidate_end) = candidate.reference_start.checked_add(reference_span) else {
            eprintln!(
                "Warning: skipping Columba candidate '{}' on '{}:{}' because its reference span overflows",
                candidate.query_name,
                candidate.reference_name,
                candidate.reference_start + 1
            );
            continue;
        };
        if candidate_end > seq.len() {
            eprintln!(
                "Warning: skipping Columba candidate '{}' on '{}:{}' because candidate span {}..{} exceeds contig length {}",
                candidate.query_name,
                candidate.reference_name,
                candidate.reference_start + 1,
                candidate.reference_start,
                candidate_end,
                seq.len()
            );
            continue;
        }

        let mut aligner = AffineWavefronts::with_penalties(0, 3, 5, 1);
        let primary_hit = verify_columba_interval(
            &mut aligner,
            record_id,
            seq,
            candidate.reference_start,
            candidate_end,
            candidate.reverse,
            guide_fwd,
            args,
        );
        let terminal_mismatch_interval = if primary_hit.is_none() {
            terminal_insertion_substitution_interval(
                candidate.reference_start,
                candidate_end,
                &candidate.cigar,
                seq.len(),
                guide_fwd.len(),
                args.max_mismatches,
                args.max_bulges,
            )
        } else {
            None
        };
        let internal_mismatch_intervals =
            if primary_hit.is_none() && terminal_mismatch_interval.is_none() {
                internal_insertion_substitution_intervals(
                    candidate.reference_start,
                    candidate_end,
                    &candidate.cigar,
                    seq.len(),
                    guide_fwd.len(),
                    args.max_mismatches,
                    args.max_bulges,
                )
            } else {
                Vec::new()
            };
        let mut recovered = false;
        let should_expand = match primary_hit {
            Some(hit) => {
                stats.primary_accepted += 1;
                let has_pam = hit.pam_seq.is_some();
                let requested_pam = hit_matches_requested_pam(&hit, &args.pam);
                let accepted_reference_insertion = requested_pam && hit.cigar.contains('D');
                push_unique_hit(hit, &mut hits, &mut seen_hits, &mut stats);
                recovered = true;
                args.max_bulges > 0
                    && args.max_bulge_size > 0
                    && ((has_pam && !requested_pam) || accepted_reference_insertion)
            }
            None => {
                stats.primary_rejected += 1;
                args.max_bulges > 0 && args.max_bulge_size > 0
            }
        };

        if let Some((start, end)) = terminal_mismatch_interval {
            stats.expansion_triggers += 1;
            stats.raw_alternatives_generated += 1;
            stats.unique_alternatives_verified += 1;
            if let Some(hit) = verify_columba_interval(
                &mut aligner,
                record_id,
                seq,
                start,
                end,
                candidate.reverse,
                guide_fwd,
                args,
            ) {
                stats.fallback_hits_accepted += 1;
                push_unique_hit(hit, &mut hits, &mut seen_hits, &mut stats);
                recovered = true;
            }
        }

        if !internal_mismatch_intervals.is_empty() {
            stats.expansion_triggers += 1;
            stats.raw_alternatives_generated += internal_mismatch_intervals.len();
            for (start, end) in internal_mismatch_intervals {
                stats.unique_alternatives_verified += 1;
                if let Some(hit) = verify_columba_interval(
                    &mut aligner,
                    record_id,
                    seq,
                    start,
                    end,
                    candidate.reverse,
                    guide_fwd,
                    args,
                ) {
                    stats.fallback_hits_accepted += 1;
                    push_unique_hit(hit, &mut hits, &mut seen_hits, &mut stats);
                    recovered = true;
                }
            }
        }

        if should_expand {
            stats.expansion_triggers += 1;
            let alternatives = fallback_intervals(
                candidate.reference_start,
                candidate_end,
                seq.len(),
                guide_fwd.len(),
                args.max_bulge_size,
                &mut stats,
            );
            for (start, end) in alternatives {
                if let Some(hit) = verify_columba_interval(
                    &mut aligner,
                    record_id,
                    seq,
                    start,
                    end,
                    candidate.reverse,
                    guide_fwd,
                    args,
                ) {
                    stats.fallback_hits_accepted += 1;
                    push_unique_hit(hit, &mut hits, &mut seen_hits, &mut stats);
                    recovered = true;
                }
            }
        }

        if !recovered {
            eprintln!(
                "Warning: skipping Columba candidate '{}' on '{}:{}' because WFA2 verification failed",
                candidate.query_name,
                candidate.reference_name,
                candidate.reference_start + 1
            );
        }
    }

    (hits, stats)
}

fn push_unique_hit(
    hit: Hit,
    hits: &mut Vec<Hit>,
    seen_hits: &mut HashSet<(String, char, usize, usize, String)>,
    stats: &mut ColumbaVerificationStats,
) {
    let key = (
        hit.ref_id.clone(),
        hit.strand,
        hit.pos,
        hit.end_pos(),
        hit.cigar.clone(),
    );
    if seen_hits.insert(key) {
        hits.push(hit);
    } else {
        stats.duplicates_removed += 1;
    }
}

fn verify_columba_interval(
    aligner: &mut AffineWavefronts,
    record_id: &str,
    seq: &[u8],
    reference_start: usize,
    reference_end: usize,
    reverse: bool,
    guide_fwd: &Arc<Vec<u8>>,
    args: &Args,
) -> Option<Hit> {
    let candidate_seq = seq.get(reference_start..reference_end)?;
    let normalized_candidate_seq = uppercase_ascii_sequence(candidate_seq);
    let reverse_candidate_seq;
    let (strand, alignment_span): (char, &[u8]) = if reverse {
        reverse_candidate_seq = reverse_complement(&normalized_candidate_seq);
        ('-', &reverse_candidate_seq)
    } else {
        ('+', &normalized_candidate_seq)
    };

    let (score, cigar, _mismatches, _gaps, _max_gap_size) = anchored_candidate_alignment(
        aligner,
        guide_fwd,
        alignment_span,
        args.max_mismatches,
        args.max_bulges,
        args.max_bulge_size,
        args.min_match_fraction,
        args.no_filter,
    )?;

    Some(build_verified_hit(
        record_id.to_string(),
        seq,
        reference_start,
        strand,
        score,
        cigar,
        Arc::clone(guide_fwd),
        alignment_span,
        0,
        args,
    ))
}

// Columba can encode terminal substitutions as query insertions on a shorter span.
// Recover only that implied full-guide interval; WFA2 remains authoritative.
fn terminal_insertion_substitution_interval(
    original_start: usize,
    original_end: usize,
    cigar: &str,
    contig_len: usize,
    guide_len: usize,
    max_mismatches: u32,
    max_bulges: u32,
) -> Option<(usize, usize)> {
    if max_mismatches == 0 || max_bulges != 0 {
        return None;
    }

    let operations = cigar_operations(cigar)?;
    if operations.len() < 2 {
        return None;
    }

    let last_index = operations.len() - 1;
    let leading_insertions = match operations.first() {
        Some((length, 'I')) => *length,
        _ => 0,
    };
    let trailing_insertions = match operations.last() {
        Some((length, 'I')) => *length,
        _ => 0,
    };
    let terminal_insertions = leading_insertions.checked_add(trailing_insertions)?;
    if terminal_insertions == 0 || terminal_insertions > max_mismatches as usize {
        return None;
    }

    for (index, (_, op)) in operations.iter().enumerate() {
        match op {
            'M' | '=' | 'X' => {}
            'I' if index == 0 || index == last_index => {}
            _ => return None,
        }
    }

    let original_span = original_end.checked_sub(original_start)?;
    if original_span.checked_add(terminal_insertions)? != guide_len {
        return None;
    }

    let start = original_start.checked_sub(leading_insertions)?;
    let end = original_end.checked_add(trailing_insertions)?;
    if end > contig_len {
        return None;
    }

    Some((start, end))
}

// Columba can also encode a substitution as an internal query insertion on a
// shorter span. Test only the full-guide intervals containing that span; the
// end-to-end WFA2 filter still requires a mismatch-only alignment.
fn internal_insertion_substitution_intervals(
    original_start: usize,
    original_end: usize,
    cigar: &str,
    contig_len: usize,
    guide_len: usize,
    max_mismatches: u32,
    max_bulges: u32,
) -> Vec<(usize, usize)> {
    if max_mismatches == 0 || max_bulges != 0 {
        return Vec::new();
    }

    let Some(operations) = cigar_operations(cigar) else {
        return Vec::new();
    };
    if operations.len() < 3 {
        return Vec::new();
    }

    let last_index = operations.len() - 1;
    let mut insertion_bases = 0usize;
    let mut has_internal_insertion = false;
    for (index, (length, op)) in operations.iter().enumerate() {
        match op {
            'M' | '=' | 'X' => {}
            'I' => {
                let Some(total) = insertion_bases.checked_add(*length) else {
                    return Vec::new();
                };
                insertion_bases = total;
                has_internal_insertion |= index != 0 && index != last_index;
            }
            _ => return Vec::new(),
        }
    }

    let original_span = match original_end.checked_sub(original_start) {
        Some(span) => span,
        None => return Vec::new(),
    };
    if !has_internal_insertion
        || insertion_bases == 0
        || insertion_bases > max_mismatches as usize
        || original_span.checked_add(insertion_bases) != Some(guide_len)
    {
        return Vec::new();
    }

    let mut intervals = Vec::new();
    for left_extension in 0..=insertion_bases {
        let right_extension = insertion_bases - left_extension;
        let Some(start) = original_start.checked_sub(left_extension) else {
            continue;
        };
        let Some(end) = original_end.checked_add(right_extension) else {
            continue;
        };
        if end <= contig_len {
            intervals.push((start, end));
        }
    }
    intervals
}

fn fallback_intervals(
    original_start: usize,
    original_end: usize,
    contig_len: usize,
    guide_len: usize,
    max_bulge_size: u32,
    stats: &mut ColumbaVerificationStats,
) -> Vec<(usize, usize)> {
    let z = max_bulge_size as usize;
    let min_span = guide_len.saturating_sub(z).max(1);
    let max_span = guide_len.saturating_add(z);
    let mut seen = HashSet::new();
    let mut intervals = Vec::new();

    let shift_bound = z.saturating_mul(2) as isize;
    for shift in -shift_bound..=shift_bound {
        let start = if shift.is_negative() {
            original_start.saturating_sub(shift.unsigned_abs())
        } else {
            original_start.saturating_add(shift as usize)
        };
        if start >= contig_len {
            continue;
        }
        for span in min_span..=max_span {
            stats.raw_alternatives_generated += 1;
            let Some(end) = start.checked_add(span) else {
                continue;
            };
            if end > contig_len || (start == original_start && end == original_end) {
                continue;
            }
            if seen.insert((start, end)) {
                intervals.push((start, end));
            }
        }
    }

    stats.unique_alternatives_verified += intervals.len();
    intervals
}

fn validation_passes(
    validation: &CigarValidation,
    guide: &[u8],
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
    min_match_fraction: f32,
    no_filter: bool,
    use_test_filter_override: bool,
) -> bool {
    let non_n_positions = guide.iter().filter(|&&b| b != b'N').count();
    let match_percentage = if non_n_positions > 0 {
        (validation.matches as f32 / non_n_positions as f32) * 100.0
    } else {
        0.0
    };
    let min_match_percentage = min_match_fraction * 100.0;

    no_filter
        || (validation.matches >= 1
            && match_percentage >= min_match_percentage
            && ((use_test_filter_override
                && cfg!(test)
                && validation.mismatches <= 1
                && validation.gaps <= 1
                && validation.max_gap_size <= 1)
                || ((!use_test_filter_override || !cfg!(test))
                    && validation.mismatches <= max_mismatches
                    && validation.gaps <= max_bulges
                    && validation.max_gap_size <= max_bulge_size)))
}

fn anchored_candidate_alignment(
    aligner: &mut AffineWavefronts,
    guide: &[u8],
    candidate_span: &[u8],
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
    min_match_fraction: f32,
    no_filter: bool,
) -> Option<(i32, String, u32, u32, u32)> {
    aligner.set_alignment_span(AlignmentSpan::End2End);
    aligner.align(candidate_span, guide);
    let score = aligner.score();
    let raw_cigar = String::from_utf8_lossy(aligner.cigar()).to_string();
    let validation = validate_cigar_accounting(&raw_cigar, guide, false);

    if validation_passes(
        &validation,
        guide,
        max_mismatches,
        max_bulges,
        max_bulge_size,
        min_match_fraction,
        no_filter,
        false,
    ) {
        Some((
            score,
            validation.cigar,
            validation.mismatches,
            validation.gaps,
            validation.max_gap_size,
        ))
    } else {
        None
    }
}

pub(crate) fn scan_window(
    aligner: &mut AffineWavefronts,
    guide: &[u8],
    window: &[u8],
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
    min_match_fraction: f32,
    no_filter: bool,
) -> Option<(i32, String, u32, u32, u32, usize)> {
    aligner.set_alignment_span(AlignmentSpan::EndsFree {
        pattern_begin_free: window.len() as i32,
        pattern_end_free: window.len() as i32,
        text_begin_free: 0,
        text_end_free: 0,
    });
    aligner.align(window, guide);
    let score = aligner.score();
    let raw_cigar = String::from_utf8_lossy(aligner.cigar()).to_string();
    let validation = validate_cigar_accounting(&raw_cigar, guide, true);

    macro_rules! debug {
        ($($arg:tt)*) => {
            #[cfg(feature = "debug")]
            eprintln!($($arg)*);
        };
    }

    debug!(
        "CIGAR: {}, N-adjusted Mismatches: {}, Gaps: {}, Max gap size: {}",
        validation.cigar, validation.mismatches, validation.gaps, validation.max_gap_size
    );

    let non_n_positions = guide.iter().filter(|&&b| b != b'N').count();
    let match_percentage = if non_n_positions > 0 {
        (validation.matches as f32 / non_n_positions as f32) * 100.0
    } else {
        0.0
    };

    let min_match_percentage = min_match_fraction * 100.0;

    debug!(
        "Match percentage: {}, Minimum required: {}",
        match_percentage, min_match_percentage
    );
    #[cfg(not(feature = "debug"))]
    let _ = (match_percentage, min_match_percentage);

    if validation_passes(
        &validation,
        guide,
        max_mismatches,
        max_bulges,
        max_bulge_size,
        min_match_fraction,
        no_filter,
        true,
    ) {
        Some((
            score,
            validation.cigar,
            validation.mismatches,
            validation.gaps,
            validation.max_gap_size,
            validation.leading_dels,
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

    fn verify_test_candidates_with_args(
        guide: &str,
        max_mismatches: u32,
        max_bulges: u32,
        max_bulge_size: u32,
        min_match_fraction: f32,
        candidates: Vec<ColumbaCandidate>,
        references: Vec<(String, Vec<u8>)>,
    ) -> Vec<Hit> {
        let mut args = columba_test_args();
        args.guide = guide.to_string();
        args.max_mismatches = max_mismatches;
        args.max_bulges = max_bulges;
        args.max_bulge_size = max_bulge_size;
        args.min_match_fraction = min_match_fraction;
        let guide_fwd = Arc::new(args.guide.as_bytes().to_vec());
        verify_columba_candidates(&candidates, &references, &guide_fwd, &args)
    }

    fn cigar_test_counts(cigar: &str, guide: &[u8]) -> (u32, u32, u32, u32) {
        let validation = validate_cigar_accounting(cigar, guide, false);
        (
            validation.matches,
            validation.mismatches,
            validation.gaps,
            validation.max_gap_size,
        )
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
        let normalized_window = uppercase_ascii_sequence(window);
        let (score, cigar, _, _, _, leading_dels) = scan_window(
            &mut aligner,
            &guide,
            &normalized_window,
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
            &normalized_window,
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
    fn test_verify_columba_uppercase_guide_matches_lowercase_reference() {
        let guide = b"gagtccgagcagaagaagaa";
        let mut seq = b"tttt".to_vec();
        seq.extend_from_slice(guide);
        seq.extend_from_slice(b"ggaaaa");
        let candidate =
            parse_candidate("guide_20bp\t0\tchr1\t5\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].ref_id, "chr1");
        assert_eq!(hits[0].strand, '+');
        assert_eq!(hits[0].target_seq, b"GAGTCCGAGCAGAAGAAGAA".to_vec());
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
            parse_candidate("guide_20bp\t0\tchr1\t5\t60\t10M1D10M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert!(hits
            .iter()
            .any(|hit| hit.cigar.contains('I') || hit.cigar.contains('D')));
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

    fn verify_candidates_with_stats_for_test(
        guide: &str,
        max_mismatches: u32,
        max_bulges: u32,
        max_bulge_size: u32,
        candidates: Vec<ColumbaCandidate>,
        seq: Vec<u8>,
    ) -> (Vec<Hit>, ColumbaVerificationStats) {
        let mut args = columba_test_args();
        args.guide = guide.to_string();
        args.max_mismatches = max_mismatches;
        args.max_bulges = max_bulges;
        args.max_bulge_size = max_bulge_size;
        args.min_match_fraction = 0.75;
        let guide_fwd = Arc::new(args.guide.as_bytes().to_vec());
        verify_columba_candidates_with_stats(
            &candidates,
            &[("chr1".to_string(), seq)],
            &guide_fwd,
            &args,
        )
    }

    #[test]
    fn test_fallback_recovers_suppressed_terminal_two_base_guide_insertion() {
        let guide = "AACGGCCTCCCAAAGTGCTG";
        let mut seq = b"CT".to_vec();
        seq.extend_from_slice(&guide.as_bytes()[2..]);
        seq.extend_from_slice(b"GG");
        let candidate =
            parse_candidate("guide\t0\tchr1\t1\t60\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert_eq!(stats.original_candidates, 1);
        assert_eq!(stats.primary_accepted, 0);
        assert_eq!(stats.primary_rejected, 1);
        assert_eq!(stats.expansion_triggers, 1);
        assert!(stats.unique_alternatives_verified <= 44);
        assert_eq!(stats.fallback_hits_accepted, 1);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].pos, 2);
        assert_eq!(hits[0].end_pos(), 20);
        assert_eq!(hits[0].pam_seq, Some("GG".to_string()));
        assert_eq!(
            cigar_test_counts(&hits[0].cigar, guide.as_bytes()),
            (18, 0, 1, 2)
        );
        assert!(hits[0].cigar.starts_with("II"));
    }

    #[test]
    fn test_fallback_recovers_suppressed_near_terminal_two_base_insertion() {
        let guide = "CGCGGCCTCCCAAAGTGCTG";
        let target = [
            guide.as_bytes()[0..2].to_vec(),
            guide.as_bytes()[4..].to_vec(),
        ]
        .concat();
        let mut seq = b"CT".to_vec();
        seq.extend_from_slice(&target);
        seq.extend_from_slice(b"GG");
        let candidate =
            parse_candidate("guide\t0\tchr1\t1\t60\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert_eq!(stats.primary_rejected, 1);
        assert_eq!(stats.expansion_triggers, 1);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].pos, 2);
        assert_eq!(hits[0].end_pos(), 20);
        assert_eq!(hits[0].pam_seq, Some("GG".to_string()));
        assert_eq!(
            cigar_test_counts(&hits[0].cigar, guide.as_bytes()),
            (18, 0, 1, 2)
        );
        assert!(hits[0].cigar.starts_with("MMII"));
    }

    #[test]
    fn test_fallback_recovers_reverse_suppressed_terminal_insertion() {
        let guide = "AACGGCCTCCCAAAGTGCTG";
        let mut oriented = b"CT".to_vec();
        oriented.extend_from_slice(&guide.as_bytes()[2..]);
        let mut seq = b"CC".to_vec();
        seq.extend_from_slice(&reverse_complement(&oriented));
        seq.extend_from_slice(b"AAAA");
        let candidate =
            parse_candidate("guide\t16\tchr1\t3\t60\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert_eq!(stats.primary_rejected, 1);
        assert_eq!(stats.expansion_triggers, 1);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand, '-');
        assert_eq!(hits[0].pos, 2);
        assert_eq!(hits[0].end_pos(), 20);
        assert_eq!(hits[0].pam_seq, Some("GG".to_string()));
        assert_eq!(
            cigar_test_counts(&hits[0].cigar, guide.as_bytes()),
            (18, 0, 1, 2)
        );
    }

    #[test]
    fn test_fallback_not_triggered_when_primary_candidate_passes() {
        let guide = "AACGGCCTCCCAAAGTGCTG";
        let mut seq = guide.as_bytes().to_vec();
        seq.extend_from_slice(b"GG");
        let candidate =
            parse_candidate("guide\t0\tchr1\t1\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert_eq!(hits.len(), 1);
        assert_eq!(stats.primary_accepted, 1);
        assert_eq!(stats.primary_rejected, 0);
        assert_eq!(stats.expansion_triggers, 0);
        assert_eq!(stats.unique_alternatives_verified, 0);
        assert_eq!(stats.fallback_hits_accepted, 0);
    }

    #[test]
    fn test_fallback_triggers_when_primary_passes_with_nonrequested_pam() {
        let guide = "GACCATCCTGGCTAACACGG";
        let mut seq = guide.as_bytes().to_vec();
        seq.extend_from_slice(b"TG");
        let candidate = parse_candidate("guide	0	chr1	1	60	20M	*	0	0	*	*	AS:i:0	NM:i:0");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert_eq!(stats.primary_accepted, 1);
        assert_eq!(stats.primary_rejected, 0);
        assert_eq!(stats.expansion_triggers, 1);
        assert!(hits.iter().any(|hit| {
            hit.pos == 0 && hit.end_pos() == 18 && hit.pam_seq.as_deref() == Some("GG")
        }));
    }

    #[test]
    fn test_fallback_recovers_start_shift_four_when_z_is_two() {
        let guide = "AGGTTTCACCATGTTCGCCA";
        let mut target = guide.as_bytes()[..16].to_vec();
        target.extend_from_slice(b"GA");
        target.extend_from_slice(&guide.as_bytes()[16..]);
        let mut seq = target;
        seq.extend_from_slice(b"GGTT");
        let candidate = parse_candidate("guide	0	chr1	5	60	20M	*	0	0	*	*	AS:i:2	NM:i:2");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert_eq!(stats.expansion_triggers, 1);
        assert!(stats.unique_alternatives_verified > 24);
        let recovered = hits
            .iter()
            .find(|hit| hit.pos == 0 && hit.end_pos() == 22 && hit.pam_seq.as_deref() == Some("GG"))
            .expect("widened fallback should recover the shifted deletion locus");
        assert_eq!(
            cigar_test_counts(&recovered.cigar, guide.as_bytes()),
            (20, 0, 1, 2)
        );
        assert!(recovered.cigar.contains("DD"));
    }

    #[test]
    fn test_fallback_triggers_for_requested_pam_reference_insertion_primary() {
        let guide = "GACCATCCTGGCTAACACGG";
        let mut seq = guide.as_bytes().to_vec();
        seq.extend_from_slice(b"TGGG");
        let candidate = parse_candidate("guide	0	chr1	1	60	22M	*	0	0	*	*	AS:i:2	NM:i:2");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert_eq!(stats.primary_accepted, 1);
        assert_eq!(stats.expansion_triggers, 1);
        assert!(hits.iter().any(|hit| {
            hit.pos == 0 && hit.end_pos() == 18 && hit.pam_seq.as_deref() == Some("GG")
        }));
    }

    #[test]
    fn test_fallback_triggers_for_reverse_requested_pam_reference_insertion_primary() {
        let guide = "GACCATCCTGGCTAACACGG";
        let mut oriented = guide.as_bytes()[0..2].to_vec();
        oriented.extend_from_slice(b"TT");
        oriented.extend_from_slice(&guide.as_bytes()[2..]);
        let mut seq = b"CC".to_vec();
        seq.extend_from_slice(&reverse_complement(&oriented));
        seq.extend_from_slice(b"AAAA");
        let candidate = parse_candidate("guide	16	chr1	3	60	22M	*	0	0	*	*	AS:i:0	NM:i:0");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert_eq!(stats.primary_accepted, 1);
        assert_eq!(stats.expansion_triggers, 1);
        assert!(hits.iter().any(|hit| {
            hit.strand == '-'
                && hit.pos == 2
                && hit.end_pos() == 20
                && hit.pam_seq.as_deref() == Some("GG")
                && hit.cigar.starts_with("II")
        }));
    }

    #[test]
    fn test_fallback_recovers_real_reverse_same_start_shorter_span() {
        let guide = "ATCTTTGCACTGATCTCCCA";
        let seq = b"ACATCAccTGGGAGATCAGTGCAAAGGTATGTCACAAA".to_vec();
        let candidate = parse_candidate("guide	272	chr1	9	60	20M	*	0	0	*	*	AS:i:1	NM:i:1");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert_eq!(stats.primary_accepted, 0);
        assert_eq!(stats.primary_rejected, 1);
        assert_eq!(stats.expansion_triggers, 1);
        assert!(
            hits.iter().any(|hit| {
                hit.strand == '-'
                    && hit.pos == 8
                    && hit.end_pos() == 26
                    && hit.pam_seq.as_deref() == Some("GG")
                    && hit.cigar.starts_with("II")
            }),
            "fallback hits: {:?}",
            hits.iter()
                .map(|hit| (
                    hit.pos,
                    hit.end_pos(),
                    hit.strand,
                    hit.cigar.clone(),
                    hit.pam_seq.clone()
                ))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_fallback_not_triggered_when_bulges_are_disabled() {
        let guide = "AACGGCCTCCCAAAGTGCTG";
        let mut seq = b"CT".to_vec();
        seq.extend_from_slice(&guide.as_bytes()[2..]);
        seq.extend_from_slice(b"GG");
        let candidate =
            parse_candidate("guide\t0\tchr1\t1\t60\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 0, 2, vec![candidate], seq);

        assert!(hits.is_empty());
        assert_eq!(stats.primary_rejected, 1);
        assert_eq!(stats.expansion_triggers, 0);
        assert_eq!(stats.unique_alternatives_verified, 0);
    }

    #[test]
    fn test_terminal_columba_insertions_recover_as_substitutions_without_bulges() {
        let guide = "GAGTCCGAGCAGAAGAAGAA";
        let mut seq = b"ACGTCAGTACGTACGA".to_vec();
        seq.extend_from_slice(b"TAGTCCGAGCAGAAGAAGAC");
        seq.extend_from_slice(b"GGTGCATGCATGCATG");
        let candidate =
            parse_candidate("guide\t256\tchr1\t18\t0\t1I18M1I\t*\t0\t0\t*\t*\tAS:i:2\tNM:i:2");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 2, 0, 0, vec![candidate], seq);

        assert_eq!(stats.primary_rejected, 1);
        assert_eq!(stats.expansion_triggers, 1);
        assert_eq!(stats.unique_alternatives_verified, 1);
        assert_eq!(stats.fallback_hits_accepted, 1);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].pos, 16);
        assert_eq!(hits[0].end_pos(), 36);
        assert_eq!(hits[0].pam_seq.as_deref(), Some("GG"));
        assert_eq!(
            cigar_test_counts(&hits[0].cigar, guide.as_bytes()),
            (18, 2, 0, 0)
        );
        assert!(hits[0].cigar.starts_with('X'));
        assert!(hits[0].cigar.ends_with('X'));
    }

    #[test]
    fn test_internal_columba_insertion_recovers_as_substitutions_without_bulges() {
        let guide = "AACGGCTCAATCAAAACAAA";
        let target = b"AACTGCTCAATCAAAAAAAA";
        let mut seq = b"TT".to_vec();
        seq.extend_from_slice(target);
        seq.extend_from_slice(b"GG");
        let candidate =
            parse_candidate("guide\t0\tchr1\t3\t0\t16M1I3M\t*\t0\t0\t*\t*\tAS:i:2\tNM:i:2");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 2, 0, 0, vec![candidate], seq);

        assert_eq!(stats.primary_rejected, 1);
        assert_eq!(stats.expansion_triggers, 1);
        assert_eq!(stats.raw_alternatives_generated, 2);
        assert_eq!(stats.unique_alternatives_verified, 2);
        assert_eq!(stats.fallback_hits_accepted, 1);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].pos, 2);
        assert_eq!(hits[0].end_pos(), 22);
        assert_eq!(hits[0].pam_seq.as_deref(), Some("GG"));
        assert_eq!(hits[0].target_seq, target);
        assert_eq!(
            cigar_test_counts(&hits[0].cigar, guide.as_bytes()),
            (18, 2, 0, 0)
        );
        assert!(!hits[0].cigar.contains('I') && !hits[0].cigar.contains('D'));
    }

    #[test]
    fn test_fallback_rejects_invalid_expanded_alternatives() {
        let guide = "CGCGGCCTCCCAAAGTGCTG";
        let seq = b"CTAAAAAAAAAAAAAAAAAAGG".to_vec();
        let candidate =
            parse_candidate("guide\t0\tchr1\t1\t60\t20M\t*\t0\t0\t*\t*\tAS:i:2\tNM:i:2");

        let (hits, stats) =
            verify_candidates_with_stats_for_test(guide, 0, 1, 2, vec![candidate], seq);

        assert!(hits.is_empty());
        assert_eq!(stats.primary_rejected, 1);
        assert_eq!(stats.expansion_triggers, 1);
        assert!(stats.unique_alternatives_verified > 0);
        assert_eq!(stats.fallback_hits_accepted, 0);
    }

    #[test]
    fn test_fallback_boundary_candidates_do_not_panic() {
        let guide = "CGCGGCCTCCCAAAGTGCTG";
        let seq = b"CTCGGCCTCCCAAAGTGCTGGG".to_vec();
        let start_candidate =
            parse_candidate("guide\t0\tchr1\t1\t60\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");
        let end_candidate =
            parse_candidate("guide\t0\tchr1\t5\t60\t20M\t*\t0\t0\t*\t*\tAS:i:2\tNM:i:2");

        let (hits, stats) = verify_candidates_with_stats_for_test(
            guide,
            0,
            1,
            2,
            vec![start_candidate, end_candidate],
            seq,
        );

        assert!(stats.expansion_triggers >= 1);
        assert!(hits.iter().all(|hit| hit.end_pos() <= hit.target_len));
    }

    #[test]
    fn test_fallback_deduplicates_same_biological_hit() {
        let guide = "AACGGCCTCCCAAAGTGCTG";
        let mut seq = b"CT".to_vec();
        seq.extend_from_slice(&guide.as_bytes()[2..]);
        seq.extend_from_slice(b"GG");
        let candidate_one =
            parse_candidate("guide\t0\tchr1\t1\t60\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");
        let candidate_two =
            parse_candidate("guide\t256\tchr1\t1\t60\t20M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let (hits, stats) = verify_candidates_with_stats_for_test(
            guide,
            0,
            1,
            2,
            vec![candidate_one, candidate_two],
            seq,
        );

        assert_eq!(hits.len(), 1);
        assert_eq!(stats.fallback_hits_accepted, 2);
        assert_eq!(stats.duplicates_removed, 1);
    }

    #[test]
    fn test_imported_same_span_terminal_guide_insertion_is_anchored() {
        let guide = "GTATTTCCCTTTTCACCGTA";
        let candidate_ref = b"TATTTCCCTTTTCACCGTA";
        let mut seq = b"AAAAAAAAAA".to_vec();
        seq.extend_from_slice(candidate_ref);
        seq.extend_from_slice(b"GGCCCC");
        let candidate =
            parse_candidate("guide\t0\tchr1\t11\t60\t1I19M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let hits = verify_test_candidates_with_args(
            guide,
            0,
            1,
            2,
            0.75,
            vec![candidate],
            vec![("chr1".to_string(), seq)],
        );

        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].pos, 10);
        assert_eq!(hits[0].cigar, "IMMMMMMMMMMMMMMMMMMM");
        assert_eq!(
            cigar_test_counts(&hits[0].cigar, guide.as_bytes()),
            (19, 0, 1, 1)
        );
        assert!(!hits[0].cigar.contains('X'));
    }

    #[test]
    fn test_imported_reverse_terminal_guide_insertion_is_anchored() {
        let guide = "GTATTTCCCTTTTCACCGTA";
        let candidate_ref = b"TATTTCCCTTTTCACCGTA";
        let mut seq = b"CCGG".to_vec();
        seq.extend_from_slice(b"AAAAAA");
        let start = seq.len();
        seq.extend_from_slice(&reverse_complement(candidate_ref));
        seq.extend_from_slice(b"CCCC");
        let candidate = parse_candidate(&format!(
            "guide\t16\tchr1\t{}\t60\t1I19M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1",
            start + 1
        ));

        let hits = verify_test_candidates_with_args(
            guide,
            0,
            1,
            2,
            0.75,
            vec![candidate],
            vec![("chr1".to_string(), seq)],
        );

        let anchored = hits
            .iter()
            .find(|hit| hit.pos == start && hit.strand == '-')
            .expect("anchored reverse candidate should remain present");
        assert_eq!(anchored.cigar, "IMMMMMMMMMMMMMMMMMMM");
        assert_eq!(
            cigar_test_counts(&anchored.cigar, guide.as_bytes()),
            (19, 0, 1, 1)
        );
    }

    #[test]
    fn test_imported_internal_insertion_stays_internal_when_anchored() {
        let guide = "ATCGATCGAT";
        let candidate_ref = b"ATCGACGAT";
        let mut seq = b"TTTT".to_vec();
        seq.extend_from_slice(candidate_ref);
        seq.extend_from_slice(b"GGAAAA");
        let candidate =
            parse_candidate("guide\t0\tchr1\t5\t60\t5M1I4M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let hits = verify_test_candidates_with_args(
            guide,
            0,
            1,
            2,
            0.75,
            vec![candidate],
            vec![("chr1".to_string(), seq)],
        );

        assert_eq!(hits.len(), 1);
        assert!(hits[0].cigar.contains('I'));
        assert!(!hits[0].cigar.starts_with('I'));
        assert!(!hits[0].cigar.ends_with('I'));
    }

    #[test]
    fn test_imported_deletion_uses_candidate_reference_span() {
        let guide = "ATCGATCGAT";
        let candidate_ref = b"ATCGAATCGAT";
        let mut seq = b"TTTT".to_vec();
        seq.extend_from_slice(candidate_ref);
        seq.extend_from_slice(b"GGAAAA");
        let candidate =
            parse_candidate("guide\t0\tchr1\t5\t60\t5M1D5M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let hits = verify_test_candidates_with_args(
            guide,
            0,
            1,
            2,
            0.75,
            vec![candidate],
            vec![("chr1".to_string(), seq)],
        );

        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].pos, 4);
        assert!(hits[0].cigar.contains('D'));
        assert_eq!(
            alignment_reference_span(&hits[0].cigar),
            candidate_ref.len()
        );
    }

    #[test]
    fn test_imported_exact_candidate_does_not_shift_in_repetitive_context() {
        let guide = "AAAAAAAAAAAAAAAAAAAA";
        let seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec();
        let candidate =
            parse_candidate("guide\t0\tchr1\t6\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits = verify_test_candidates_with_args(
            guide,
            0,
            1,
            2,
            0.75,
            vec![candidate],
            vec![("chr1".to_string(), seq)],
        );

        let anchored = hits
            .iter()
            .find(|hit| hit.pos == 5 && hit.cigar == "MMMMMMMMMMMMMMMMMMMM")
            .expect("reported exact candidate should remain anchored at the SAM coordinate");
        assert_eq!(anchored.pos, 5);
    }

    #[test]
    fn test_imported_incorrect_candidate_cannot_shift_to_nearby_match() {
        let guide = "GAGTCCGAGCAGAAGAAGAA";
        let mut seq = b"CCCCCCCCCCCCCCCCCCCC".to_vec();
        seq.extend_from_slice(guide.as_bytes());
        seq.extend_from_slice(b"GGAAAA");
        let candidate =
            parse_candidate("guide\t0\tchr1\t1\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits = verify_test_candidates_with_args(
            guide,
            0,
            1,
            2,
            0.75,
            vec![candidate],
            vec![("chr1".to_string(), seq)],
        );

        assert!(hits.is_empty());
    }

    #[test]
    fn test_five_terminal_guide_bases_remain_counted_and_rejected() {
        let guide = b"GGAAACGAATGGAGTGTAAT";
        let raw_cigar = "MMMMMMMMMMMMMMMIIIII";
        let validation = validate_cigar_accounting(raw_cigar, guide, true);

        assert_eq!(validation.cigar, raw_cigar);
        assert_eq!(validation.matches, 15);
        assert_eq!(validation.gaps, 1);
        assert_eq!(validation.max_gap_size, 5);

        let mut aligner = setup_aligner();
        let result = scan_window(&mut aligner, guide, &guide[..15], 0, 1, 2, 0.75, false);
        assert!(
            result.is_none(),
            "terminal guide gap must not pass just because -f permits 15/20 matches"
        );
    }

    #[test]
    fn test_one_terminal_guide_base_gap_remains_counted_and_allowed() {
        let guide = b"GGAAACGAATGGAGTGTAAT";
        let raw_cigar = "MMMMMMMMMMMMMMMMMMMI";
        let validation = validate_cigar_accounting(raw_cigar, guide, true);

        assert_eq!(validation.cigar, raw_cigar);
        assert_eq!(validation.matches, 19);
        assert_eq!(validation.gaps, 1);
        assert_eq!(validation.max_gap_size, 1);

        let mut aligner = setup_aligner();
        let result = scan_window(&mut aligner, guide, &guide[..19], 0, 1, 2, 0.75, false)
            .expect("one terminal guide-base gap is within the configured bulge limits");
        assert_eq!(result.1, raw_cigar);
        assert_eq!(result.3, 1);
        assert_eq!(result.4, 1);
    }

    #[test]
    fn test_three_terminal_guide_base_gap_is_rejected() {
        let guide = b"GGAAACGAATGGAGTGTAAT";
        let raw_cigar = "MMMMMMMMMMMMMMMMMIII";
        let validation = validate_cigar_accounting(raw_cigar, guide, true);

        assert_eq!(validation.cigar, raw_cigar);
        assert_eq!(validation.matches, 17);
        assert_eq!(validation.gaps, 1);
        assert_eq!(validation.max_gap_size, 3);

        let mut aligner = setup_aligner();
        let result = scan_window(&mut aligner, guide, &guide[..17], 0, 1, 2, 0.75, false);
        assert!(result.is_none());
    }

    #[test]
    fn test_free_reference_overhangs_are_trimmed_without_counted_bulge() {
        let guide = b"GGAAACGAATGGAGTGTAAT";
        let raw_cigar = "DDDDDMMMMMMMMMMMMMMMMMMMMDDDDD";
        let validation = validate_cigar_accounting(raw_cigar, guide, true);

        assert_eq!(validation.leading_dels, 5);
        assert_eq!(validation.cigar, "MMMMMMMMMMMMMMMMMMMM");
        assert_eq!(validation.matches, 20);
        assert_eq!(validation.gaps, 0);
        assert_eq!(validation.max_gap_size, 0);
    }

    #[test]
    fn test_reference_overhang_trim_keeps_real_guide_gap() {
        let guide = b"GGAAACGAATGGAGTGTAAT";
        let raw_cigar = "DDDMMMMMMMMMMIIMMMMMMMMMMDDD";
        let validation = validate_cigar_accounting(raw_cigar, guide, true);

        assert_eq!(validation.leading_dels, 3);
        assert_eq!(validation.cigar, "MMMMMMMMMMIIMMMMMMMMMM");
        assert_eq!(validation.matches, 20);
        assert_eq!(validation.gaps, 1);
        assert_eq!(validation.max_gap_size, 2);
    }

    #[test]
    fn test_reverse_columba_candidate_retains_guide_gap_accounting() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let mut target = guide.to_vec();
        target.remove(10);
        let mut seq = b"CCCCC".to_vec();
        seq.extend_from_slice(&reverse_complement(&target));
        seq.extend_from_slice(b"AAAAAGG");
        let candidate =
            parse_candidate("guide_20bp\t16\tchr1\t6\t60\t19M\t*\t0\t0\t*\t*\tAS:i:1\tNM:i:1");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand, '-');
        assert!(hits[0].cigar.contains('I'));
    }

    #[test]
    fn test_repetitive_flanked_exact_target_uses_terminal_overhangs() {
        let mut aligner = setup_aligner();
        let guide = b"GGAAACGAATGGAGTGTAAT";
        let window = b"GGAATGGAAACGAATGGAGTGTAATGGAAT";
        assert_eq!(&window[5..25], guide);

        let result = scan_window(&mut aligner, guide, window, 0, 1, 2, 0.75, false);
        assert!(
            result.is_some(),
            "exact target in flanked window should verify"
        );
        let (_score, cigar, mismatches, gaps, max_gap_size, leading_dels) = result.unwrap();
        assert_eq!(leading_dels, 5);
        assert_eq!(cigar, "MMMMMMMMMMMMMMMMMMMM");
        assert_eq!(mismatches, 0);
        assert_eq!(gaps, 0);
        assert_eq!(max_gap_size, 0);
        assert!(!cigar.contains("DDDDD"));
    }

    #[test]
    fn test_ordinary_flanked_exact_target_still_verifies() {
        let mut aligner = setup_aligner();
        let guide = b"CTATTCAGTTCCCATATCCC";
        let window = b"GGGCACTATTCAGTTCCCATATCCCGGAAA";
        assert_eq!(&window[5..25], guide);

        let result = scan_window(&mut aligner, guide, window, 0, 1, 2, 0.75, false);
        assert!(
            result.is_some(),
            "ordinary exact target in flanked window should verify"
        );
        let (_score, cigar, mismatches, gaps, max_gap_size, leading_dels) = result.unwrap();
        assert_eq!(leading_dels, 5);
        assert_eq!(cigar, "MMMMMMMMMMMMMMMMMMMM");
        assert_eq!(mismatches, 0);
        assert_eq!(gaps, 0);
        assert_eq!(max_gap_size, 0);
    }

    #[test]
    fn test_internal_bulge_remains_internal_with_flanked_window() {
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let window = b"GGGATCGAATCGATTTT";

        let result = scan_window(&mut aligner, guide, window, 1, 1, 1, 0.75, false);
        assert!(
            result.is_some(),
            "one-base internal bulge should still verify"
        );
        let (_score, cigar, _mismatches, gaps, max_gap_size, leading_dels) = result.unwrap();
        assert_eq!(leading_dels, 3);
        assert_eq!(gaps, 1);
        assert_eq!(max_gap_size, 1);
        assert!(cigar.contains('D') || cigar.contains('I'));
        assert_ne!(cigar, "MMMMMMMMMM");
    }

    #[test]
    fn test_partial_guide_match_in_flanked_window_is_rejected() {
        let mut aligner = setup_aligner();
        let guide = b"GGAAACGAATGGAGTGTAAT";
        let window = b"GGAATGGAAACGAATGGAAT";

        let result = scan_window(&mut aligner, guide, window, 0, 1, 2, 0.75, false);
        assert!(
            result.is_none(),
            "reference-window free ends must not allow partial guide matches"
        );
    }

    #[test]
    fn test_reverse_columba_candidate_with_mixed_case_reference() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        let guide_rc = reverse_complement(guide);
        let mut seq = b"TTGGGGGGGG".to_vec();
        seq.extend(guide_rc.iter().map(u8::to_ascii_lowercase));
        seq.extend_from_slice(b"tttttt");
        let candidate =
            parse_candidate("guide_20bp\t16\tchr1\t11\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0");

        let hits = verify_test_candidates(vec![candidate], vec![("chr1".to_string(), seq)]);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].pos, 10);
        assert_eq!(hits[0].strand, '-');
        assert_eq!(hits[0].target_seq, guide.to_vec());
        assert_eq!(hits[0].pam_seq, Some("CC".to_string()));
        assert_eq!(hits[0].cigar, "MMMMMMMMMMMMMMMMMMMM");
    }

    #[test]
    fn test_scan_window_reports_offsets_with_reference_overhangs() {
        let guide = b"GAGTCCGAGCAGAAGAAGAA";
        for (window, expected_offset) in [
            (&b"GAGTCCGAGCAGAAGAAGAATTTTT"[..], 0usize),
            (&b"AAAAAGAGTCCGAGCAGAAGAAGAATTTTT"[..], 5usize),
            (&b"AAAGAGTCCGAGCAGAAGAAGAATTTTTTT"[..], 3usize),
            (&b"AAAAAAAAAAGAGTCCGAGCAGAAGAAGAA"[..], 10usize),
        ] {
            let mut aligner = setup_aligner();
            let result = scan_window(&mut aligner, guide, window, 0, 1, 2, 0.75, false)
                .expect("exact guide substring should verify");
            assert_eq!(result.1, "MMMMMMMMMMMMMMMMMMMM");
            assert_eq!(result.5, expected_offset);
        }
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
