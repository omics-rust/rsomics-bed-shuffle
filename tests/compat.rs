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

/// A 300 kb feature fits chr1 (1000 kb) and chr2 (500 kb) but not chr3
/// (250 kb); it must be retried onto a chromosome that fits and always emit
/// exactly one line — never silently dropped. bedtools v2.31.1 emits exactly
/// one line for this input on every seed.
#[test]
fn feature_fitting_only_some_chromosomes_is_always_placed() {
    let genome = load_genome();
    let feat_len = 300_000u64;
    let bed = "chr1\t0\t300000\n";
    for seed in [1u64, 2, 3, 4, 5, 99, 777, 12345] {
        let out = run(bed, false, seed);
        assert_eq!(out.len(), 1, "300 kb feature must be placed (seed={seed})");
        let c: Vec<&str> = out[0].split('\t').collect();
        let chrom = c[0];
        let start: u64 = c[1].parse().unwrap();
        let end: u64 = c[2].parse().unwrap();
        let chrom_len = *genome.get(chrom).unwrap();
        assert_eq!(end - start, feat_len, "length preserved (seed={seed})");
        assert!(end <= chrom_len, "in bounds (seed={seed})");
        assert!(
            chrom == "chr1" || chrom == "chr2",
            "must land on a chromosome that fits, got {chrom} (seed={seed})"
        );
    }
}

/// A feature longer than every chromosome cannot be placed: bedtools warns to
/// stderr and drops it (rc=0, 0 lines). We do the same — no panic, no output.
#[test]
fn feature_longer_than_all_chromosomes_is_dropped() {
    let bed = "chr1\t0\t2000000\n";
    let out = run(bed, false, 1);
    assert!(out.is_empty(), "unplaceable feature emits no line");
}

/// An empty / zero-total-length genome plus a feature must fail loud, not panic
/// with a divide-by-zero. (bedtools itself does not terminate cleanly here; a
/// loud error is the correct fail-fast behavior.)
#[test]
fn empty_genome_errors_loud() {
    let empty: HashMap<String, u64> = HashMap::new();
    let mut out = Vec::new();
    let res = shuffle(
        Cursor::new("chr1\t0\t100\n"),
        &mut out,
        &empty,
        &ShuffleOptions {
            same_chrom: false,
            allow_beyond_chrom_end: false,
            seed: Some(1),
        },
    );
    assert!(res.is_err(), "empty genome must error, not panic");
}

/// A genome file of only comments / blank lines parses to an empty map and must
/// likewise error loud rather than divide by zero.
#[test]
fn comment_only_genome_errors_loud() {
    let genome = parse_genome(Cursor::new("# just a comment\n\n")).unwrap();
    assert!(genome.is_empty());
    let mut out = Vec::new();
    let res = shuffle(
        Cursor::new("chr1\t0\t100\n"),
        &mut out,
        &genome,
        &ShuffleOptions {
            same_chrom: false,
            allow_beyond_chrom_end: false,
            seed: Some(1),
        },
    );
    assert!(res.is_err(), "comment-only genome must error, not panic");
}
