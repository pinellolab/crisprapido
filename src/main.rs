use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;

use bio::io::fasta;
use clap::Parser;
use flate2::read::MultiGzDecoder;
use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentSpan};
use rayon::prelude::*;

mod cfd_score;
mod columba;
mod hit;
mod reporting;
mod verification;

use columba::{parse_columba_sam_file, run_columba, ColumbaRunConfig};
use hit::Hit;
use reporting::report_filtered_hits;
use verification::{
    build_verified_hit, load_reference_records, scan_window, uppercase_ascii_sequence,
    verify_columba_candidates,
};

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            _ => b'N',
        })
        .collect()
}

#[derive(Parser)]
#[command(author, version, about = "CRISPR guide RNA off-target scanner")]
pub(crate) struct Args {
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

    /// Existing Columba SAM file to parse for candidate locations
    #[arg(long)]
    columba_sam: Option<PathBuf>,

    /// Path to Columba executable
    #[arg(long)]
    columba_bin: Option<PathBuf>,

    /// Existing Columba FM-index prefix
    #[arg(long)]
    columba_index: Option<PathBuf>,

    /// Preserve temporary Columba SAM output
    #[arg(long)]
    keep_columba_sam: bool,
}

fn open_reference_reader(path: &PathBuf) -> fasta::Reader<BufReader<Box<dyn BufRead>>> {
    let file = File::open(path).expect("Failed to open reference file");
    let reader: Box<dyn BufRead> = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    fasta::Reader::new(reader)
}

fn main() {
    let args = Args::parse();

    if let Err(e) = cfd_score::init_score_matrices(
        args.mismatch_scores
            .to_str()
            .unwrap_or("mismatch_scores.txt"),
        args.pam_scores.to_str().unwrap_or("pam_scores.txt"),
    ) {
        eprintln!("Warning: CFD scoring disabled - {}", e);
    }

    let mut aligner = AffineWavefronts::with_penalties(0, 3, 5, 1);
    aligner.set_alignment_span(AlignmentSpan::EndsFree {
        pattern_begin_free: 1,
        pattern_end_free: 1,
        text_begin_free: 1,
        text_end_free: 1,
    });

    let guide_fwd = Arc::new(args.guide.as_bytes().to_vec());
    let guide_rc = Arc::new(reverse_complement(&guide_fwd));
    let guide_len = guide_fwd.len();

    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to initialize thread pool");
    }

    if args.columba_sam.is_some() || args.columba_bin.is_some() || args.columba_index.is_some() {
        let reference_records = load_reference_records(&args.reference).unwrap_or_else(|e| {
            eprintln!("{}", e);
            std::process::exit(1);
        });

        let candidates = if let Some(path) = &args.columba_sam {
            parse_columba_sam_file(path).unwrap_or_else(|e| {
                eprintln!("Failed to parse Columba SAM file: {}", e);
                std::process::exit(1);
            })
        } else {
            let columba_bin = args.columba_bin.as_ref().unwrap_or_else(|| {
                eprintln!("--columba-bin is required when --columba-index is supplied");
                std::process::exit(1);
            });
            let columba_index = args.columba_index.as_ref().unwrap_or_else(|| {
                eprintln!("--columba-index is required when --columba-bin is supplied");
                std::process::exit(1);
            });
            let columba_output = run_columba(&ColumbaRunConfig {
                columba_bin,
                index_prefix: columba_index,
                guide: &args.guide,
                max_mismatches: args.max_mismatches,
                threads: args.threads,
                keep_sam: args.keep_columba_sam,
            })
            .unwrap_or_else(|e| {
                eprintln!("{}", e);
                std::process::exit(1);
            });
            parse_columba_sam_file(&columba_output.sam_path).unwrap_or_else(|e| {
                eprintln!("Failed to parse generated Columba SAM file: {}", e);
                std::process::exit(1);
            })
        };

        let hits = verify_columba_candidates(&candidates, &reference_records, &guide_fwd, &args);
        report_filtered_hits(hits, &args.pam);
        return;
    }

    let reader = open_reference_reader(&args.reference);
    for result in reader.records() {
        let record = result.expect("Error during FASTA record parsing");
        let seq = record.seq().to_vec();
        let record_id = record.id().to_string();

        let window_size = args.window_size.unwrap_or(guide_len * 4);
        let step_size = window_size / 2;
        let windows: Vec<_> = (0..seq.len())
            .step_by(step_size)
            .map(|i| (i, (i + window_size).min(seq.len())))
            .collect();

        let hits: Vec<Hit> = windows
            .into_par_iter()
            .map_init(
                || AffineWavefronts::with_penalties(0, 3, 5, 1),
                |aligner, (i, end)| {
                    let window = &seq[i..end];
                    if window.len() < guide_len {
                        return None;
                    }
                    let normalized_window = uppercase_ascii_sequence(window);

                    if let Some((score, cigar, _mismatches, _gaps, _max_gap_size, leading_dels)) =
                        scan_window(
                            aligner,
                            &guide_fwd,
                            &normalized_window,
                            args.max_mismatches,
                            args.max_bulges,
                            args.max_bulge_size,
                            args.min_match_fraction,
                            args.no_filter,
                        )
                    {
                        return Some(build_verified_hit(
                            record_id.clone(),
                            &seq,
                            i + leading_dels,
                            '+',
                            score,
                            cigar.clone(),
                            Arc::clone(&guide_fwd),
                            &normalized_window,
                            leading_dels,
                            &args,
                        ));
                    }

                    if let Some((score, cigar, _mismatches, _gaps, _max_gap_size, leading_dels)) =
                        scan_window(
                            aligner,
                            &guide_rc,
                            &normalized_window,
                            args.max_mismatches,
                            args.max_bulges,
                            args.max_bulge_size,
                            args.min_match_fraction,
                            args.no_filter,
                        )
                    {
                        return Some(build_verified_hit(
                            record_id.clone(),
                            &seq,
                            i + leading_dels,
                            '-',
                            score,
                            cigar.clone(),
                            Arc::clone(&guide_rc),
                            &normalized_window,
                            leading_dels,
                            &args,
                        ));
                    }

                    None
                },
            )
            .filter_map(|x| x)
            .collect();

        report_filtered_hits(hits, &args.pam);
    }
}
