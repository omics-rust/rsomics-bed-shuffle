//! Correctness tests for rsomics-bed-shuffle.
//!
//! bedtools shuffle uses a different RNG so we cannot byte-match its output.
//! Instead we verify: count, length preservation, bounds, and determinism.

use rsomics_bed_shuffle::{ShuffleOptions, parse_genome, shuffle};
use std::collections::HashMap;
use std::io::Cursor;
use std::path::Path;

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn load_genome() -> HashMap<String, u64> {
    let f = std::fs::File::open(golden("genome.sizes")).unwrap();
    parse_genome(f).unwrap()
}

fn run(bed: &str, same_chrom: bool, seed: u64) -> Vec<String> {
    let g = load_genome();
    let mut out = Vec::new();
    shuffle(
        Cursor::new(bed),
        &mut out,
        &g,
        &ShuffleOptions {
            same_chrom,
            allow_beyond_chrom_end: false,
            seed: Some(seed),
        },
    )
    .unwrap();
    String::from_utf8(out)
        .unwrap()
        .lines()
        .map(str::to_owned)
        .collect()
}

#[test]
fn output_count_equals_input_count() {
    let bed = std::fs::read_to_string(golden("intervals.bed")).unwrap();
    let out = run(&bed, false, 1);
    let input_count = bed.lines().count();
    assert_eq!(
        out.len(),
        input_count,
        "output count must equal input count"
    );
}

#[test]
fn feature_lengths_preserved() {
    let bed = std::fs::read_to_string(golden("intervals.bed")).unwrap();
    let input_lens: Vec<u64> = bed
        .lines()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            c[2].parse::<u64>().unwrap() - c[1].parse::<u64>().unwrap()
        })
        .collect();
    let out = run(&bed, false, 42);
    let output_lens: Vec<u64> = out
        .iter()
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            c[2].parse::<u64>().unwrap() - c[1].parse::<u64>().unwrap()
        })
        .collect();
    assert_eq!(
        input_lens, output_lens,
        "feature lengths must be preserved after shuffle"
    );
}

#[test]
fn coords_within_genome_bounds() {
    let bed = std::fs::read_to_string(golden("intervals.bed")).unwrap();
    let genome = load_genome();
    for seed in [1u64, 42, 99, 777, 12345] {
        let out = run(&bed, false, seed);
        for line in &out {
            let c: Vec<&str> = line.split('\t').collect();
            let chrom = c[0];
            let end: u64 = c[2].parse().unwrap();
            let chrom_len = *genome.get(chrom).unwrap();
            assert!(
                end <= chrom_len,
                "end {end} > chrom_len {chrom_len} (chrom={chrom}, seed={seed})"
            );
        }
    }
}

#[test]
fn deterministic_with_seed() {
    let bed = std::fs::read_to_string(golden("intervals.bed")).unwrap();
    let a = run(&bed, false, 999);
    let b = run(&bed, false, 999);
    assert_eq!(a, b, "same seed must produce same output");
}
