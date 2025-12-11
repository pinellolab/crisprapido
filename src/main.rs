use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::{Arc, atomic::{AtomicBool, AtomicU64, Ordering}};
use std::time::{Duration, Instant};
use flate2::read::MultiGzDecoder;
use clap::Parser;
use bio::io::fasta;
use sassy::profiles::Iupac;
use sassy::Searcher;
use rayon::prelude::*;
use crossbeam_channel::{bounded, Sender};
use std::thread;
use std::cell::RefCell;

mod cfd_score;

#[derive(Default)]
struct PerfCounters {
    hits_sent: AtomicU64,
    hits_written: AtomicU64,
    bytes_written: AtomicU64,
    send_block_ns: AtomicU64,
    send_block_count: AtomicU64,
    task_total_ns: AtomicU64,
    task_count: AtomicU64,
    queue_max_depth: AtomicU64,
}

impl PerfCounters {
    fn record_hit_sent(&self, blocked: Duration) {
        if blocked.as_nanos() > 0 {
            self.send_block_ns.fetch_add(blocked.as_nanos() as u64, Ordering::Relaxed);
            self.send_block_count.fetch_add(1, Ordering::Relaxed);
        }
        let sent = self.hits_sent.fetch_add(1, Ordering::Relaxed) + 1;
        let written = self.hits_written.load(Ordering::Relaxed);
        let depth = sent.saturating_sub(written);
        self.update_max_depth(depth);
    }

    fn record_hit_written(&self) {
        self.hits_written.fetch_add(1, Ordering::Relaxed);
    }

    fn record_bytes(&self, bytes: usize) {
        self.bytes_written.fetch_add(bytes as u64, Ordering::Relaxed);
    }

    fn record_task_duration(&self, duration: Duration) {
        self.task_total_ns.fetch_add(duration.as_nanos() as u64, Ordering::Relaxed);
        self.task_count.fetch_add(1, Ordering::Relaxed);
    }

    fn update_max_depth(&self, depth: u64) {
        loop {
            let current = self.queue_max_depth.load(Ordering::Relaxed);
            if depth <= current {
                break;
            }
            if self
                .queue_max_depth
                .compare_exchange(current, depth, Ordering::Relaxed, Ordering::Relaxed)
                .is_ok()
            {
                break;
            }
        }
    }

    fn print_summary(&self) {
        let hits = self.hits_written.load(Ordering::Relaxed);
        let bytes = self.bytes_written.load(Ordering::Relaxed);
        let tasks = self.task_count.load(Ordering::Relaxed);
        let avg_task_ms = if tasks > 0 {
            let total_ns = self.task_total_ns.load(Ordering::Relaxed);
            (total_ns as f64 / tasks as f64) / 1_000_000.0
        } else {
            0.0
        };
        let avg_block_ms = {
            let blocks = self.send_block_count.load(Ordering::Relaxed);
            if blocks > 0 {
                let total_ns = self.send_block_ns.load(Ordering::Relaxed);
                (total_ns as f64 / blocks as f64) / 1_000_000.0
            } else {
                0.0
            }
        };
        eprintln!(
            "[perf] total_hits={} total_bytes={:.2}MB avg_task_ms={:.2} avg_send_wait_ms={:.3} max_queue_depth={}",
            hits,
            bytes as f64 / (1024.0 * 1024.0),
            avg_task_ms,
            avg_block_ms,
            self.queue_max_depth.load(Ordering::Relaxed)
        );
    }
}

struct CountingWriter<W: Write> {
    inner: W,
    perf: Option<Arc<PerfCounters>>,
}

impl<W: Write> CountingWriter<W> {
    fn new(inner: W, perf: Option<Arc<PerfCounters>>) -> Self {
        Self { inner, perf }
    }
}

impl<W: Write> Write for CountingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buf)?;
        if let Some(perf) = &self.perf {
            perf.record_bytes(written);
        }
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

fn spawn_perf_monitor(
    perf: Arc<PerfCounters>,
    stop_flag: Arc<AtomicBool>,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        let interval = Duration::from_secs(5);
        let mut last_hits = 0;
        let mut last_bytes = 0;
        let mut last_tasks = 0;
        let mut last_blocks = 0;
        while !stop_flag.load(Ordering::Relaxed) {
            thread::sleep(interval);
            if stop_flag.load(Ordering::Relaxed) {
                break;
            }
            let hits = perf.hits_written.load(Ordering::Relaxed);
            let bytes = perf.bytes_written.load(Ordering::Relaxed);
            let tasks = perf.task_count.load(Ordering::Relaxed);
            let blocks = perf.send_block_count.load(Ordering::Relaxed);
            let interval_secs = interval.as_secs_f64();
            let hits_rate = (hits - last_hits) as f64 / interval_secs;
            let bytes_rate = (bytes - last_bytes) as f64 / (1024.0 * 1024.0) / interval_secs;
            let task_rate = (tasks - last_tasks) as f64 / interval_secs;
            let block_rate = (blocks - last_blocks) as f64 / interval_secs;
            let avg_task_ms = if tasks > 0 {
                let total_ns = perf.task_total_ns.load(Ordering::Relaxed);
                (total_ns as f64 / tasks as f64) / 1_000_000.0
            } else {
                0.0
            };
            eprintln!(
                "[perf] hits/s={:.1} MB/s={:.2} tasks/s={:.1} queue_depth={} max_queue={} block/s={:.1} avg_task_ms={:.2}",
                hits_rate,
                bytes_rate,
                task_rate,
                perf.hits_sent.load(Ordering::Relaxed)
                    .saturating_sub(perf.hits_written.load(Ordering::Relaxed)),
                perf.queue_max_depth.load(Ordering::Relaxed),
                block_rate,
                avg_task_ms,
            );
            last_hits = hits;
            last_bytes = bytes;
            last_tasks = tasks;
            last_blocks = blocks;
        }
    })
}


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

// Structure for passing output data through the channel
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

fn report_hit(writer: &mut impl Write, ref_id: &str, pos: usize, _len: usize, strand: char,
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
    let _ = writeln!(writer, "Guide\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tas:i:{}\tnm:i:{}\tng:i:{}\tbs:i:{}\tcg:Z:{}{}{}",
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
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGAT";  // Single mismatch at position 5

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty(), "Should accept a single mismatch");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = &results[0];
        assert_eq!(cigar, "4=1X5=");
    }

    #[test]
    fn test_with_bulge_sassy() {
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGAATCGAT";  // Single base insertion after position 4

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty(), "Should accept a single base bulge");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = &results[0];
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences_sassy() {
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9

        let results = scan_contig_sassy(guide, target, 1, 1, 1, 0.75, false, true);
        assert!(results.is_empty());
    }

    #[test]
    fn test_perfect_match_with_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let target = create_flanked_sequence(&mut rng, guide, 500);

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty(), "Should match perfectly even with flanks");

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
        let core = b"ATCGTTCGAT";  // Single mismatch at position 5
        let target = create_flanked_sequence(&mut rng, core, 500);

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty(), "Should accept a single mismatch with flanks");

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
        let core = b"ATCGAATCGAT";  // Single base insertion after position 4
        let target = create_flanked_sequence(&mut rng, core, 500);

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty(), "Should accept a single base bulge with flanks");

        // Find the match at position 500 (SASSY may find other matches in random flanks)
        let match_at_500 = results.iter().find(|(_, _, _, _, _, pos)| *pos == 500);
        assert!(match_at_500.is_some(), "Should find match at position 500");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _pos) = match_at_500.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences_with_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        let target = create_flanked_sequence(&mut rng, core, 500);

        let results = scan_contig_sassy(guide, &target, 1, 1, 1, 0.75, false, true);

        // Should not find a match at position 500 (too many mismatches)
        // Note: SASSY may find accidental matches in the random flanks, so we only check position 500
        let match_at_500 = results.iter().find(|(_, _, _, _, _, pos)| *pos == 500);
        assert!(match_at_500.is_none(), "Should reject sequence with too many mismatches at position 500");
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

    #[test]
    fn test_handles_ns_in_target() {
        // Test that scan_contig_sassy handles N's in the target sequence (via Iupac profile)
        let guide = b"ATCGATCGAT";
        let target_with_n = b"ATCGATCNATCGATCGAT";  // Has an N in the middle

        // Should not panic - Iupac profile handles N's natively
        let results = scan_contig_sassy(guide, target_with_n, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty(), "Should handle N's in target and find match");
        assert!(results.iter().any(|(_, _, _, _, _, pos)| alignment_overlaps_ambiguous(target_with_n, *pos, guide.len())));

        // By default, ambiguous hits should be dropped
        let filtered = scan_contig_sassy(guide, target_with_n, 1, 1, 1, 0.75, false, false);
        assert!(filtered.iter().all(|(_, _, _, _, _, pos)| !alignment_overlaps_ambiguous(target_with_n, *pos, guide.len())));
    }

    #[test]
    fn test_handles_ns_in_guide() {
        // Test that scan_contig_sassy handles N's in the guide sequence (via Iupac profile)
        let guide_with_n = b"ATCNATCGAT";  // Has an N at position 3
        let target = b"ATCGATCGAT";

        // Should not panic - Iupac profile handles N's natively (N matches any base)
        let results = scan_contig_sassy(guide_with_n, target, 1, 1, 1, 0.75, false, true);
        assert!(!results.is_empty(), "Should handle N's in guide");

        let filtered = scan_contig_sassy(guide_with_n, target, 1, 1, 1, 0.75, false, false);
        assert!(filtered.is_empty(), "Ambiguous guide bases should be ignored unless requested");
    }
}

#[derive(Parser)]
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

    /// Include hits that overlap ambiguous bases (N/R/Y etc.)
    #[arg(long)]
    include_ambiguous: bool,
}


fn convert_to_minimap2_cigar(cigar: &str) -> String {
    // Remove the debug print line
    if cigar.is_empty() {
        return "".to_string();
    }

    cigar.to_string()  // Just return the CIGAR as-is
}

// Thread-local SASSY searcher to avoid repeated allocation
// Using Iupac profile for native N handling (faster than Dna + N->A conversion)
thread_local! {
    static SEARCHER: RefCell<Searcher<Iupac>> = RefCell::new(Searcher::new(false, None));
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
    let matches = SEARCHER.with(|searcher| {
        searcher.borrow_mut().search(guide, &contig, max_errors)
    });

    if matches.is_empty() {
        return Vec::new();
    }

    // Process ALL matches, not just the best one
    matches.into_iter()
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

            if no_filter || (
                matches_count >= 1 &&
                match_percentage >= min_match_fraction * 100.0 &&
                mismatches <= max_mismatches &&
                gaps <= max_bulges &&
                max_gap_size <= max_bulge_size
            ) {
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
    consolidated.iter()
        .map(|(count, op)| format!("{}{}", count, op))
        .collect::<String>()
}

fn reverse_cigar(cigar: &str) -> String {
    if cigar.is_empty() {
        return String::new();
    }

    let mut ops: Vec<(String, char)> = Vec::new();
    let mut chars = cigar.chars().peekable();
    while let Some(&ch) = chars.peek() {
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
    contig[start..end]
        .iter()
        .any(|&b| is_ambiguous_base(b))
}

fn is_ambiguous_base(b: u8) -> bool {
    matches!(
        b,
        b'N' | b'n'
            | b'R' | b'r'
            | b'Y' | b'y'
            | b'S' | b's'
            | b'W' | b'w'
            | b'K' | b'k'
            | b'M' | b'm'
            | b'B' | b'b'
            | b'D' | b'd'
            | b'H' | b'h'
            | b'V' | b'v'
    )
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

fn read_guides_from_file(path: &PathBuf) -> std::io::Result<Vec<Vec<u8>>> {
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut guides = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            guides.push(trimmed.as_bytes().to_vec());
        }
    }

    Ok(guides)
}

fn main() {
    let args = Args::parse();
    let perf_enabled = std::env::var("CRISPRAPIDO_TRACE_PERF").is_ok();
    let perf_counters = if perf_enabled {
        Some(Arc::new(PerfCounters::default()))
    } else {
        None
    };
    let perf_stop = perf_counters.as_ref().map(|_| Arc::new(AtomicBool::new(false)));
    let perf_handle = perf_counters.as_ref().zip(perf_stop.as_ref()).map(|(perf, stop)| {
        spawn_perf_monitor(perf.clone(), stop.clone())
    });

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

    // Load guides from either -g or --guides-file
    let guides: Vec<Vec<u8>> = if let Some(ref guide_str) = args.guide {
        vec![guide_str.as_bytes().to_vec()]
    } else if let Some(ref guides_path) = args.guides_file {
        read_guides_from_file(guides_path)
            .expect("Failed to read guides file")
    } else {
        eprintln!("Error: Must provide either --guide (-g) or --guides-file");
        std::process::exit(1);
    };

    if guides.is_empty() {
        eprintln!("Error: No guides provided");
        std::process::exit(1);
    }

    let guide_rcs: Vec<Vec<u8>> = guides
        .iter()
        .map(|g| reverse_complement(g))
        .collect();

    eprintln!("Loaded {} guide(s)", guides.len());

    // Set thread pool size if specified
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to initialize thread pool");
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
        fasta_reader.records()
            .map(|r| r.expect("Error during FASTA record parsing"))
            .map(ContigData::from_record)
            .collect()
    );

    eprintln!("Loaded {} contig(s)", contigs.len());

    // Create channel for sending hits from workers to output thread
    // Use bounded channel to provide backpressure and limit queued hits
    let (tx, rx) = bounded::<OutputHit>(2048);
    let perf_for_writer = perf_counters.clone();

    // Spawn single output consumer thread with buffered output
    let output_thread = thread::spawn(move || {
        let stdout = std::io::stdout();
        let writer = stdout.lock();
        let mut writer = CountingWriter::new(BufWriter::with_capacity(256 * 1024, writer), perf_for_writer);
        for hit in rx {
            hit.write_output(&mut writer);
            if let Some(perf) = &writer.perf {
                perf.record_hit_written();
            }
        }
        writer.flush().expect("Failed to flush output");
    });

    // Create flat list of all (guide_idx, contig_idx) task pairs
    // This eliminates nested parallelism and barriers
    let tasks: Vec<(usize, usize)> = (0..guides.len())
        .flat_map(|g_idx| {
            (0..contigs.len()).map(move |c_idx| (g_idx, c_idx))
        })
        .collect();

    eprintln!("Processing {} tasks (guides × contigs)", tasks.len());

    // Process all tasks in parallel - flat, no barriers!
    let perf_for_workers = perf_counters.clone();
    tasks.par_iter().for_each(|(guide_idx, contig_idx)| {
        let guide = &guides[*guide_idx];
        let rev_guide = &guide_rcs[*guide_idx];
        let contig = &contigs[*contig_idx];
        let guide_len = guide.len();

        let seq = contig.seq();
        let seq_len = seq.len();
        let record_id = contig.id();

        // Skip contigs shorter than guide
        if seq_len < guide_len {
            return;
        }

        // Scan entire contig with SASSY - returns ALL matches
        // Reserve buffers outside of scan loops to limit per-hit allocations
        let mut target_seq_buf = Vec::with_capacity(guide_len);
        let mut rc_target_seq_buf = Vec::with_capacity(guide_len);
        let mut rc_cache = Vec::with_capacity(guide_len);

        let worker_start = perf_for_workers.as_ref().map(|_| Instant::now());
        let matches = scan_contig_sassy(
            guide,
            seq,
            args.max_mismatches,
            args.max_bulges,
            args.max_bulge_size,
            args.min_match_fraction,
            args.no_filter,
            args.include_ambiguous,
        );

        for (score, cigar, _mismatches, _gaps, _max_gap_size, pos) in matches {
            // Extract target sequence slice into reusable buffer
            target_seq_buf.clear();
            let start = pos;
            let end = (pos + guide_len).min(seq_len);
            target_seq_buf.extend_from_slice(&seq[start..end]);

            if !args.include_ambiguous && target_seq_buf.iter().any(|&b| is_ambiguous_base(b)) {
                continue;
            }

            let output_hit = OutputHit {
                ref_id: record_id.to_string(),
                pos,
                guide_len,
                strand: '+',
                score,
                cigar,
                guide: guide.clone(),
                target_len: seq_len,
                max_mismatches: args.max_mismatches,
                max_bulges: args.max_bulges,
                max_bulge_size: args.max_bulge_size,
                target_seq: target_seq_buf.clone(),
                pam: args.pam.clone(),
            };

            if let Some(perf) = &perf_for_workers {
                let start = Instant::now();
                tx.send(output_hit).expect("Failed to send hit to output thread");
                let blocked = start.elapsed();
                perf.record_hit_sent(blocked);
            } else {
                tx.send(output_hit).expect("Failed to send hit to output thread");
            }
        }

        let rev_matches = scan_contig_sassy(
            rev_guide,
            seq,
            args.max_mismatches,
            args.max_bulges,
            args.max_bulge_size,
            args.min_match_fraction,
            args.no_filter,
            args.include_ambiguous,
        );

        for (score, cigar, _mismatches, _gaps, _max_gap_size, pos) in rev_matches {
            rc_target_seq_buf.clear();
            rc_cache.clear();
            let start = pos;
            let end = (pos + guide_len).min(seq_len);
            rc_cache.extend_from_slice(&seq[start..end]);
            rc_target_seq_buf.extend(reverse_complement(&rc_cache));

            if !args.include_ambiguous && rc_target_seq_buf.iter().any(|&b| is_ambiguous_base(b)) {
                continue;
            }

            let adjusted_cigar = reverse_cigar(&cigar);

            let output_hit = OutputHit {
                ref_id: record_id.to_string(),
                pos,
                guide_len,
                strand: '-',
                score,
                cigar: adjusted_cigar,
                guide: guide.clone(),
                target_len: seq_len,
                max_mismatches: args.max_mismatches,
                max_bulges: args.max_bulges,
                max_bulge_size: args.max_bulge_size,
                target_seq: rc_target_seq_buf.clone(),
                pam: args.pam.clone(),
            };

            if let Some(perf) = &perf_for_workers {
                let start = Instant::now();
                tx.send(output_hit).expect("Failed to send reverse-strand hit");
                let blocked = start.elapsed();
                perf.record_hit_sent(blocked);
            } else {
                tx.send(output_hit).expect("Failed to send reverse-strand hit");
            }
        }

        if let (Some(perf), Some(worker_start)) = (perf_for_workers.as_ref(), worker_start) {
            perf.record_task_duration(worker_start.elapsed());
        }
    });

    // Drop the original sender to signal completion
    drop(tx);

    // Wait for output thread to finish writing all results
    output_thread.join().expect("Output thread panicked");

    if let Some(stop) = perf_stop.as_ref() {
        stop.store(true, Ordering::Relaxed);
    }
    if let Some(handle) = perf_handle {
        let _ = handle.join();
    }
    if let Some(perf) = perf_counters {
        perf.print_summary();
    }
}
