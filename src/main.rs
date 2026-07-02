use bio::io::fasta;
use clap::Parser;
use crossbeam_channel::{bounded, unbounded, Sender};
use flate2::read::MultiGzDecoder;
use num_cpus;
use sassy::profiles::Iupac;
use sassy::Searcher;
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc,
};
use std::thread;

mod cfd_score;

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            _ => b'N', // Convert any unexpected bases to N
        })
        .collect()
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
    cfd_score: Option<f64>,
    target_seq: Vec<u8>,
}

#[derive(Clone)]
struct OutputHit {
    ref_id: String,
    pos: usize,
    guide_len: usize,
    strand: char,
    score: i32,
    cigar: String,
    guide: Vec<u8>,
    target_len: usize,
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
    target_seq: Vec<u8>,
    pam: String,
}

#[derive(Clone)]
struct GuideJob {
    forward: Vec<u8>,
    reverse: Vec<u8>,
}

enum GuideSource {
    Inline(Vec<Vec<u8>>),
    File(PathBuf),
}

impl OutputHit {
    fn write_output(&self, writer: &mut impl Write) {
        report_hit(
            writer,
            &self.ref_id,
            self.pos,
            self.guide_len,
            self.strand,
            self.score,
            &self.cigar,
            &self.guide,
            self.target_len,
            self.max_mismatches,
            self.max_bulges,
            self.max_bulge_size,
            &self.target_seq,
            &self.pam,
        );
    }
}

// Structure to hold contig data
struct ContigData {
    id: String,
    seq: Vec<u8>,
}

impl ContigData {
    fn from_record(record: fasta::Record) -> Self {
        ContigData {
            id: record.id().to_string(),
            seq: record.seq().to_vec(),
        }
    }

    fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn id(&self) -> &str {
        &self.id
    }

    fn len(&self) -> usize {
        self.seq.len()
    }
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
    writer: &mut impl Write,
    ref_id: &str,
    pos: usize,
    guide_len: usize,
    strand: char,
    score: i32,
    cigar: &str,
    guide: &[u8],
    target_len: usize,
    _max_mismatches: u32,
    _max_bulges: u32,
    _max_bulge_size: u32,
    target_seq: &[u8],
    pam: &str,
) {
    let (matches, mismatches, gaps, max_gap_size) = parse_cigar_stats(cigar);

    let effective_cigar = if cigar.is_empty() {
        format!("{}=", guide_len)
    } else {
        cigar.to_string()
    };

    let guide_window_len = guide_len.min(guide.len());
    let target_window_len = guide_len.min(target_seq.len());

    let guide_window = &guide[..guide_window_len];
    let target_window = &target_seq[..target_window_len];

    let spacer_for_cfd = if strand == '-' {
        reverse_complement(guide_window)
    } else {
        guide_window.to_vec()
    };
    let protospacer_for_cfd = target_window.to_vec();

    let cfd_score =
        cfd_score::get_cfd_score(&spacer_for_cfd, &protospacer_for_cfd, &effective_cigar, pam)
            .unwrap_or(0.0);

    let guide_str = String::from_utf8_lossy(guide_window);
    let target_str = String::from_utf8_lossy(target_window);

    let _ = writeln!(
        writer,
        "Guide\t{}\t0\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tas:i:{}\tnm:i:{}\tng:i:{}\tbs:i:{}\tcg:Z:{}\tcf:f:{:.4}\tqs:Z:{}\tts:Z:{}",
        guide_len,
        guide_window_len,
        strand,
        ref_id,
        target_len,
        pos,
        pos + guide_window_len,
        matches,
        score,
        mismatches,
        gaps,
        max_gap_size,
        effective_cigar,
        cfd_score,
        guide_str,
        target_str
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

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty());
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = &results[0];
        assert_eq!(cigar, "10=");
    }

    #[test]
    fn test_with_mismatches_sassy() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGTTCGAT"; // Single mismatch at position 5

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty(), "Should accept a single mismatch");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = &results[0];
        assert_eq!(cigar, "4=1X5=");
    }

    #[test]
    fn test_with_bulge_sassy() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGAATCGAT"; // Single base insertion after position 4

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false, true);
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

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false, true);
        assert!(results.is_empty());
    }

    #[test]
    fn test_perfect_match_with_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let target = create_flanked_sequence(&mut rng, guide, 500);

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false, true);
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

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false, true);
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

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false, true);
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

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false, true);

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

    #[test]
    fn test_handles_ns_in_target() {
        // Test that scan_contig_sassy handles N's in the target sequence (via Iupac profile)
        let guide = b"ATCGATCGAT";
        let target_with_n = b"ATCGATCNATCGATCGAT"; // Has an N in the middle

        // Should not panic - Iupac profile handles N's natively
        let results = scan_contig_sassy(guide, target_with_n, 1, 1, 1, 0.75, false, true);
        assert!(
            !results.is_empty(),
            "Should handle N's in target and find match"
        );
        assert!(results
            .iter()
            .any(|(_, _, _, _, _, pos)| alignment_overlaps_ambiguous(
                target_with_n,
                *pos,
                guide.len()
            )));

        // By default, ambiguous hits should be dropped
        let filtered = scan_contig_sassy(guide, target_with_n, 1, 1, 1, 0.75, false, false);
        assert!(filtered
            .iter()
            .all(|(_, _, _, _, _, pos)| !alignment_overlaps_ambiguous(
                target_with_n,
                *pos,
                guide.len()
            )));
    }

    #[test]
    fn test_handles_ns_in_guide() {
        // Test that scan_contig_sassy handles N's in the guide sequence (via Iupac profile)
        let guide_with_n = b"ATCNATCGAT"; // Has an N at position 3
        let target = b"ATCGATCGAT";

        // Should not panic - Iupac profile handles N's natively (N matches any base)
        let results = scan_contig_sassy(guide_with_n, target, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty(), "Should handle N's in guide");

        let filtered = scan_contig_sassy(guide_with_n, target, 1, 1, 1, 0.75, false, false);
        assert!(
            filtered.is_empty(),
            "Ambiguous guide bases should be ignored unless requested"
        );
    }

    #[test]
    fn test_cfd_scores_are_strand_invariant() {
        cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let guide = b"GAAACAGTCGATTTTATCAC".to_vec();
        let pam = "GG";

        let plus_line = capture_report_line('+', guide.as_slice(), guide.as_slice(), pam);
        let plus_cf = extract_cf_score(&plus_line);

        let guide_rc = reverse_complement(guide.as_slice());
        let minus_line = capture_report_line('-', guide_rc.as_slice(), guide.as_slice(), pam);
        let minus_cf = extract_cf_score(&minus_line);

        assert!(
            (plus_cf - minus_cf).abs() < 1e-9,
            "CFD scores should match on both strands: plus={}, minus={}",
            plus_cf,
            minus_cf
        );
        assert!(
            (plus_cf - 1.0).abs() < 1e-9,
            "Perfect match should stay at 1.0, got {}",
            plus_cf
        );
    }

    fn capture_report_line(strand: char, guide: &[u8], target: &[u8], pam: &str) -> String {
        let mut buffer = Vec::new();
        report_hit(
            &mut buffer,
            "chr_demo",
            0,
            guide.len(),
            strand,
            0,
            "20=",
            guide,
            target.len(),
            0,
            0,
            0,
            target,
            pam,
        );
        String::from_utf8(buffer).expect("report_hit output should be valid UTF-8")
    }

    fn extract_cf_score(line: &str) -> f64 {
        line.split('\t')
            .find_map(|field| field.strip_prefix("cf:f:"))
            .and_then(|value| value.parse::<f64>().ok())
            .expect("cf:f field missing or unparsable")
    }
}

#[derive(Parser, Debug, Clone)]
#[command(author, version, about = "CRISPR guide RNA off-target scanner")]
struct Args {
    /// Input reference FASTA file (-r)
    #[arg(short, long)]
    reference: PathBuf,

    /// Guide RNA sequence (without PAM) (-g)
    #[arg(short, long, conflicts_with = "guides_file")]
    guide: Option<String>,

    /// File containing guide RNA sequences, one per line (without PAM)
    #[arg(long = "guides-file", conflicts_with = "guide")]
    guides_file: Option<PathBuf>,

    /// PAM sequence (to use for CFD scoring)
    #[arg(short = 'p', long, default_value = "GG")]
    pam: String,

    /// Maximum number of mismatches allowed (0 = perfect match only, up to 3)
    #[arg(short, long, default_value = "3", value_parser = clap::value_parser!(u32).range(0..=3))]
    max_mismatches: u32,

    /// Maximum number of bulges (indel events) allowed — only used by sassy backend
    #[arg(short = 'b', long, default_value = "1")]
    max_bulges: u32,

    /// Maximum size of each bulge in bp — only used by sassy backend
    #[arg(short = 'z', long, default_value = "2")]
    max_bulge_size: u32,

    /// Allow insertions and deletions (indels/bulges) in matches.
    /// When absent, only SNP-style mismatches are reported (bulges filtered out).
    #[arg(long)]
    allow_indels: bool,

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

    /// Include hits that overlap ambiguous bases (N/R/Y etc.)
    #[arg(long)]
    include_ambiguous: bool,

    /// Use ropebwt3 as the search backend instead of sassy.
    /// Builds an FM-index (.fmd) from the reference once at startup,
    /// then queries it for each guide RNA. Faster for large genomes.
    #[arg(long)]
    use_ropebwt3: bool,

    /// Path to the ropebwt3 binary
    #[arg(long, default_value = "./ropebwt3/ropebwt3")]
    ropebwt3_bin: PathBuf,
}

fn convert_to_minimap2_cigar(cigar: &str) -> String {
    // Remove the debug print line
    if cigar.is_empty() {
        return "".to_string();
    }

    cigar.to_string() // Just return the CIGAR as-is
}

// Thread-local SASSY searcher to avoid repeated allocation
// Using Iupac profile for native N handling (faster than Dna + N->A conversion)
thread_local! {
    static SEARCHER: RefCell<Searcher<Iupac>> = RefCell::new(Searcher::new_fwd());
}

fn process_single_guide(
    args: &Arc<Args>,
    contigs: &Arc<Vec<ContigData>>,
    tx_hits: &Sender<OutputHit>,
    job: GuideJob,
    fmd_path: Option<Arc<PathBuf>>,
) {
    if args.use_ropebwt3 {
        // ── ropebwt3 backend ──────────────────────────────────────────────
        // ropebwt3 searches both strands internally (BWT includes reverse complement),
        // so we only query the forward guide and let ropebwt3 report strand via PAF.
        let fmd = fmd_path.as_ref().expect("fmd_path must be set when use_ropebwt3=true");
        let threads = args.threads.unwrap_or_else(num_cpus::get).max(1);

        let hits = query_ropebwt3(
            &args.ropebwt3_bin,
            fmd,
            &job.forward,
            args.max_mismatches,
            args.allow_indels,
            threads,
        );

        // Build a contig lookup map for fast seq retrieval
        let contig_map: std::collections::HashMap<&str, &ContigData> =
            contigs.iter().map(|c| (c.id(), c)).collect();

        for (ref_id, pos, cigar, score) in hits {
            // Filter indels if not allowed
            let (_, _, gaps, _) = parse_cigar_stats(&cigar);
            if !args.allow_indels && gaps > 0 {
                continue;
            }

            // Count mismatches and enforce limit
            let (_, mismatches, _, _) = parse_cigar_stats(&cigar);
            if mismatches > args.max_mismatches {
                continue;
            }

            let contig = match contig_map.get(ref_id.as_str()) {
                Some(c) => c,
                None => continue,
            };
            let seq = contig.seq();
            let guide_len = job.forward.len();
            let end = (pos + guide_len).min(seq.len());
            if pos >= seq.len() {
                continue;
            }
            let target_seq = seq[pos..end].to_vec();

            if !args.include_ambiguous && target_seq.iter().any(|&b| is_ambiguous_base(b)) {
                continue;
            }

            // ropebwt3 reports strand in PAF but we unify to '+' here (index covers both strands)
            let strand = '+';

            let output_hit = OutputHit {
                ref_id: ref_id.clone(),
                pos,
                guide_len,
                strand,
                score,
                cigar,
                guide: job.forward.clone(),
                target_len: seq.len(),
                max_mismatches: args.max_mismatches,
                max_bulges: args.max_bulges,
                max_bulge_size: args.max_bulge_size,
                target_seq,
                pam: args.pam.clone(),
            };

            tx_hits.send(output_hit).expect("Failed to send hit");
        }
    } else {
        // ── sassy backend (original) ──────────────────────────────────────
        for contig in contigs.iter() {
            let seq = contig.seq();
            if seq.len() < job.forward.len() {
                continue;
            }

            for (guide, strand) in [(&job.forward, '+'), (&job.reverse, '-')] {
                let effective_bulges = if args.allow_indels { args.max_bulges } else { 0 };
                let effective_bulge_size = if args.allow_indels { args.max_bulge_size } else { 0 };

                let matches = scan_contig_sassy(
                    guide,
                    seq,
                    args.max_mismatches,
                    effective_bulges,
                    effective_bulge_size,
                    args.min_match_fraction,
                    args.no_filter,
                    args.include_ambiguous,
                );

                for (score, cigar, _mm, gaps, _max_gap, pos) in matches {
                    // Enforce --allow-indels at the filter level
                    if !args.allow_indels && gaps > 0 {
                        continue;
                    }

                    let guide_len = guide.len();
                    let start = pos;
                    let end = (pos + guide_len).min(seq.len());
                    let mut target_seq = seq[start..end].to_vec();

                    if strand == '-' {
                        target_seq = reverse_complement(&target_seq);
                    }

                    if !args.include_ambiguous && target_seq.iter().any(|&b| is_ambiguous_base(b)) {
                        continue;
                    }

                    let output_hit = OutputHit {
                        ref_id: contig.id().to_string(),
                        pos,
                        guide_len,
                        strand,
                        score,
                        cigar: if strand == '-' {
                            reverse_cigar(&cigar)
                        } else {
                            cigar.clone()
                        },
                        guide: guide.clone(),
                        target_len: seq.len(),
                        max_mismatches: args.max_mismatches,
                        max_bulges: args.max_bulges,
                        max_bulge_size: args.max_bulge_size,
                        target_seq,
                        pam: args.pam.clone(),
                    };

                    tx_hits.send(output_hit).expect("Failed to send hit");
                }
            }
        }
    }
}

fn scan_contig_sassy(
    guide: &[u8],
    contig: &[u8],
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
    min_match_fraction: f32,
    no_filter: bool,
    include_ambiguous: bool,
) -> Vec<(i32, String, u32, u32, u32, usize)> {
    if !include_ambiguous {
        if guide.iter().any(|&b| is_ambiguous_base(b)) {
            return Vec::new();
        }
    }

    // Calculate maximum allowed errors
    let max_errors = (max_mismatches + max_bulges) as usize;

    // Use thread-local SASSY searcher (reused across calls in same thread)
    // ZERO-COPY: SASSY accepts &[u8] via RcSearchAble trait
    let matches =
        SEARCHER.with(|searcher| searcher.borrow_mut().search(guide, &contig, max_errors));

    if matches.is_empty() {
        return Vec::new();
    }

    // Process ALL matches, not just the best one
    matches
        .into_iter()
        .filter_map(|sassy_match| {
            let score = sassy_match.cost as i32;
            let pos = sassy_match.text_start as usize;

            // Use SASSY's CIGAR and normalize it to always include counts
            let cigar_str = normalize_cigar(&sassy_match.cigar.to_string());

            // Skip alignments that touch ambiguous bases unless explicitly requested
            if !include_ambiguous && alignment_overlaps_ambiguous(contig, pos, guide.len()) {
                return None;
            }

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
    while let Some(&_ch) = chars.peek() {
        let count = if _ch.is_ascii_digit() {
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

fn reverse_cigar(cigar: &str) -> String {
    if cigar.is_empty() {
        return String::new();
    }

    let mut ops: Vec<(String, char)> = Vec::new();
    let mut chars = cigar.chars().peekable();
    while let Some(&_ch) = chars.peek() {
        let mut num_str = String::new();
        while let Some(&digit_ch) = chars.peek() {
            if digit_ch.is_ascii_digit() {
                num_str.push(chars.next().unwrap());
            } else {
                break;
            }
        }

        if num_str.is_empty() {
            break;
        }

        if let Some(op) = chars.next() {
            ops.push((num_str, op));
        }
    }

    ops.reverse();
    ops.into_iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect::<String>()
}

fn alignment_overlaps_ambiguous(contig: &[u8], start: usize, guide_len: usize) -> bool {
    let end = start.saturating_add(guide_len).min(contig.len());
    contig[start..end].iter().any(|&b| is_ambiguous_base(b))
}

fn is_ambiguous_base(b: u8) -> bool {
    matches!(
        b,
        b'N' | b'n'
            | b'R'
            | b'r'
            | b'Y'
            | b'y'
            | b'S'
            | b's'
            | b'W'
            | b'w'
            | b'K'
            | b'k'
            | b'M'
            | b'm'
            | b'B'
            | b'b'
            | b'D'
            | b'd'
            | b'H'
            | b'h'
            | b'V'
            | b'v'
    )
}

fn parse_cigar_stats(cigar: &str) -> (usize, u32, u32, u32) {
    let mut matches_count = 0;
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut max_gap_size = 0;
    let mut current_gap_size = 0;

    let mut chars = cigar.chars().peekable();
    while let Some(&ch) = chars.peek() {
        if ch.is_ascii_digit() {
            let mut num_str = String::new();
            while let Some(&digit_ch) = chars.peek() {
                if digit_ch.is_ascii_digit() {
                    num_str.push(chars.next().unwrap());
                } else {
                    break;
                }
            }

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


// ─────────────────────────────────────────────────────────────────────────────
// ropebwt3 backend
// ─────────────────────────────────────────────────────────────────────────────

/// Build a ropebwt3 FM-index (.fmd) from the reference FASTA.
/// Returns the path to the .fmd file.
/// If the index already exists, skips building (fast re-runs).
fn build_ropebwt3_index(bin: &PathBuf, reference: &PathBuf) -> PathBuf {
    let fmd_path = reference.with_extension("fmd");

    if fmd_path.exists() {
        eprintln!(
            "ropebwt3: index already exists at {}, skipping build",
            fmd_path.display()
        );
        return fmd_path;
    }

    eprintln!(
        "ropebwt3: building FM-index from {} → {}",
        reference.display(),
        fmd_path.display()
    );

    let status = Command::new(bin)
        .args(["build", "-o"])
        .arg(&fmd_path)
        .arg(reference)
        .status()
        .expect("Failed to launch ropebwt3 build — check --ropebwt3-bin path");

    if !status.success() {
        eprintln!("ropebwt3 build failed with exit code: {:?}", status.code());
        std::process::exit(1);
    }

    eprintln!("ropebwt3: index built successfully");
    fmd_path
}

/// Query the ropebwt3 FM-index for one guide RNA sequence.
/// Uses `ropebwt3 sw` for inexact matching (mismatches + optional indels).
/// Returns hits as (ref_id, pos, cigar, score).
fn query_ropebwt3(
    bin: &PathBuf,
    fmd_path: &PathBuf,
    guide: &[u8],
    max_mismatches: u32,
    allow_indels: bool,
    threads: usize,
) -> Vec<(String, usize, String, i32)> {
    // Write the guide to a temp FASTA for ropebwt3 input
    let tmp_dir = std::env::temp_dir();
    let tmp_fa = tmp_dir.join(format!(
        "crisprapido_guide_{}.fa",
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .subsec_nanos()
    ));

    {
        let mut f = File::create(&tmp_fa).expect("Failed to create temp guide FASTA");
        writeln!(f, ">guide").unwrap();
        writeln!(f, "{}", String::from_utf8_lossy(guide)).unwrap();
    }

    // ropebwt3 sw: local alignment with affine-gap penalty
    // -e     : end-to-end mode (full guide must align)
    // -N<k>  : maximum edit distance (mismatches + gaps)
    // -t<n>  : threads
    let max_errors = if allow_indels {
        max_mismatches
    } else {
        max_mismatches // gaps will be filtered post-hoc when allow_indels=false
    };

    let output = Command::new(bin)
        .arg("sw")
        .arg("-e")
        .arg(format!("-N{}", max_errors))
        .arg(format!("-t{}", threads))
        .arg(fmd_path)
        .arg(&tmp_fa)
        .output()
        .expect("Failed to launch ropebwt3 sw");

    let _ = std::fs::remove_file(&tmp_fa);

    if !output.status.success() {
        eprintln!(
            "ropebwt3 sw warning: exit code {:?}",
            output.status.code()
        );
        return Vec::new();
    }

    // Parse PAF output: qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, matches, block, mapq, [tags]
    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut hits = Vec::new();

    for line in stdout.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let strand_char = fields[4];
        let ref_id = fields[5].to_string();
        let tstart: usize = fields[7].parse().unwrap_or(0);
        let matches: i32 = fields[9].parse().unwrap_or(0);
        let block_len: i32 = fields[10].parse().unwrap_or(0);
        let score = block_len - matches; // edit distance approximation

        // Extract CIGAR from cg:Z: tag if present
        let cigar = fields
            .iter()
            .skip(12)
            .find_map(|f| f.strip_prefix("cg:Z:"))
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("{}=", guide.len()));

        let cigar = normalize_cigar(&cigar);

        // Flip strand: ropebwt3 indexes both strands; strand '+' in PAF means forward
        let _strand = if strand_char == "+" { '+' } else { '-' };

        hits.push((ref_id, tstart, cigar, score));
    }

    hits
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

    // Load all reference sequences into memory (shared, immutable)
    let file = File::open(&args.reference).expect("Failed to open reference file");
    let reader: Box<dyn BufRead> = if args.reference.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let fasta_reader = fasta::Reader::new(reader);

    // Collect all contigs as Arc-wrapped ContigData for shared access
    // N's are replaced with A's and tracked in bit vectors
    let contigs: Arc<Vec<ContigData>> = Arc::new(
        fasta_reader
            .records()
            .map(|r| r.expect("Error during FASTA record parsing"))
            .map(ContigData::from_record)
            .collect(),
    );

    eprintln!("Loaded {} contig(s)", contigs.len());

    let guide_queue_size = num_cpus::get().max(2);
    let (guide_tx, guide_rx) = bounded::<GuideJob>(guide_queue_size);

    let total_guides = Arc::new(AtomicUsize::new(0));
    let guide_counter = total_guides.clone();

    let guide_source = if let Some(ref guide_str) = args.guide {
        GuideSource::Inline(vec![guide_str.as_bytes().to_vec()])
    } else if let Some(ref guides_file) = args.guides_file {
        GuideSource::File(guides_file.clone())
    } else {
        eprintln!("Error: Must provide either --guide (-g) or --guides-file");
        std::process::exit(1);
    };

    // Build ropebwt3 index once at startup if requested
    let fmd_path: Option<Arc<PathBuf>> = if args.use_ropebwt3 {
        let fmd = build_ropebwt3_index(&args.ropebwt3_bin, &args.reference);
        Some(Arc::new(fmd))
    } else {
        None
    };

    let consumer_contigs = contigs.clone();
    let consumer_args = Arc::new(args.clone());
    let (tx, rx) = unbounded::<OutputHit>();

    // Spawn output thread first
    let output_thread = thread::spawn(move || {
        let stdout = std::io::stdout();
        let mut writer = BufWriter::with_capacity(256 * 1024, stdout.lock());
        for hit in rx {
            hit.write_output(&mut writer);
        }
        writer.flush().expect("Failed to flush output");
    });

    // Producer: feed guides into bounded queue
    let producer_handle = thread::spawn(move || {
        let enqueue = |guide: Vec<u8>| {
            let reverse = reverse_complement(&guide);
            guide_counter.fetch_add(1, Ordering::Relaxed);
            guide_tx
                .send(GuideJob {
                    forward: guide,
                    reverse,
                })
                .ok();
        };

        match guide_source {
            GuideSource::Inline(list) => {
                for guide in list {
                    enqueue(guide);
                }
            }
            GuideSource::File(path) => {
                let file = File::open(&path).expect("Failed to open guides file");
                let reader: Box<dyn BufRead> = if path.extension().map_or(false, |ext| ext == "gz")
                {
                    Box::new(BufReader::new(MultiGzDecoder::new(file)))
                } else {
                    Box::new(BufReader::new(file))
                };

                for line in reader.lines().filter_map(Result::ok) {
                    let trimmed = line.trim();
                    if trimmed.is_empty() || trimmed.starts_with('#') {
                        continue;
                    }
                    enqueue(trimmed.as_bytes().to_vec());
                }
            }
        }
        // close the channel to signal completion
        drop(guide_tx);
    });

    let worker_count = args.threads.unwrap_or_else(num_cpus::get);
    let mut workers = Vec::with_capacity(worker_count);

    for _ in 0..worker_count {
        let rx_guides = guide_rx.clone();
        let contigs = consumer_contigs.clone();
        let tx_hits = tx.clone();
        let args = Arc::clone(&consumer_args);
        let fmd = fmd_path.clone();

        workers.push(thread::spawn(move || {
            for job in rx_guides {
                process_single_guide(&args, &contigs, &tx_hits, job, fmd.clone());
            }
        }));
    }

    drop(tx);

    for worker in workers {
        let _ = worker.join();
    }

    let _ = producer_handle.join();
    output_thread.join().expect("Output thread panicked");
}
