//! Benchmark comparing Dna vs Iupac profile performance
//! Run with: cargo bench --bench dna_vs_iupac

use rand::{rngs::SmallRng, RngCore, SeedableRng};
use sassy::profiles::{Dna, Iupac};
use sassy::Searcher;
use std::time::Instant;

fn generate_random_seq_with_ns(rng: &mut SmallRng, length: usize, n_frequency: f64) -> Vec<u8> {
    let bases = b"ACGT";
    (0..length)
        .map(|_| {
            if (rng.next_u32() as f64 / u32::MAX as f64) < n_frequency {
                b'N'
            } else {
                bases[rng.next_u32() as usize % 4]
            }
        })
        .collect()
}

fn replace_ns_with_a(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| if b == b'N' { b'A' } else { b })
        .collect()
}

fn main() {
    let mut rng = SmallRng::seed_from_u64(42);

    // Generate test data: 10MB sequence with ~1% Ns (realistic for genomic data)
    let seq_len = 10_000_000;
    let n_freq = 0.01;
    let target_with_ns = generate_random_seq_with_ns(&mut rng, seq_len, n_freq);
    let target_no_ns = replace_ns_with_a(&target_with_ns);

    let n_count = target_with_ns.iter().filter(|&&b| b == b'N').count();
    println!("Sequence length: {} bp", seq_len);
    println!(
        "N count: {} ({:.2}%)",
        n_count,
        100.0 * n_count as f64 / seq_len as f64
    );

    // Generate guide (no Ns)
    let guide = b"ATCGATCGATCGATCGATCG"; // 20bp guide

    let max_errors = 4;
    let iterations = 10;

    println!(
        "\nBenchmarking {} iterations with max_errors={}...\n",
        iterations, max_errors
    );

    // Benchmark Dna profile (requires N->A conversion)
    let mut dna_searcher: Searcher<Dna> = Searcher::new_fwd();
    let start = Instant::now();
    let mut dna_matches = 0;
    for _ in 0..iterations {
        let results = dna_searcher.search(guide, &target_no_ns, max_errors);
        dna_matches = results.len();
    }
    let dna_time = start.elapsed();
    println!("Dna profile (with N->A conversion):");
    println!("  Total time: {:?}", dna_time);
    println!("  Per iteration: {:?}", dna_time / iterations);
    println!("  Matches found: {}", dna_matches);

    // Benchmark Iupac profile (native N support)
    let mut iupac_searcher: Searcher<Iupac> = Searcher::new_fwd();
    let start = Instant::now();
    let mut iupac_matches = 0;
    for _ in 0..iterations {
        let results = iupac_searcher.search(guide, &target_with_ns, max_errors);
        iupac_matches = results.len();
    }
    let iupac_time = start.elapsed();
    println!("\nIupac profile (native N handling):");
    println!("  Total time: {:?}", iupac_time);
    println!("  Per iteration: {:?}", iupac_time / iterations);
    println!("  Matches found: {}", iupac_matches);

    // Compare
    let speedup = iupac_time.as_secs_f64() / dna_time.as_secs_f64();
    println!("\n--- Comparison ---");
    if speedup > 1.0 {
        println!("Dna is {:.2}x faster than Iupac", speedup);
    } else {
        println!("Iupac is {:.2}x faster than Dna", 1.0 / speedup);
    }

    // Also benchmark the N->A conversion overhead
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = replace_ns_with_a(&target_with_ns);
    }
    let conversion_time = start.elapsed();
    println!("\nN->A conversion overhead:");
    println!("  Total time: {:?}", conversion_time);
    println!("  Per iteration: {:?}", conversion_time / iterations);
}
