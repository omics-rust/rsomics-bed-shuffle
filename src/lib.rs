//! Randomly relocate BED intervals within a genome — bedtools shuffle equivalent.
//!
//! Chromosomes are weighted by length when picking a random target. Use `--seed`
//! for reproducible output.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read, Write};

use rsomics_common::{Result, RsomicsError};

pub fn parse_genome<R: Read>(r: R) -> Result<HashMap<String, u64>> {
    let mut map = HashMap::new();
    for line in BufReader::new(r).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut cols = t.splitn(2, '\t');
        let chrom = cols
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput(format!("genome parse error: {t}")))?
            .to_owned();
        let len: u64 = cols
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput(format!("genome missing size: {t}")))?
            .trim()
            .parse()
            .map_err(|e| RsomicsError::InvalidInput(format!("genome size parse error: {e}")))?;
        map.insert(chrom, len);
    }
    Ok(map)
}

pub struct ShuffleOptions {
    pub same_chrom: bool,
    pub allow_beyond_chrom_end: bool,
    pub seed: Option<u64>,
}

/// Features too long to fit any chromosome are skipped (matches bedtools).
#[allow(clippy::implicit_hasher)]
pub fn shuffle<R: Read, W: Write>(
    r: R,
    w: W,
    genome: &HashMap<String, u64>,
    opts: &ShuffleOptions,
) -> Result<()> {
    let mut chroms: Vec<(&str, u64)> = genome.iter().map(|(k, &v)| (k.as_str(), v)).collect();
    chroms.sort_unstable_by_key(|(c, _)| *c);

    let total_len: u64 = chroms.iter().map(|(_, l)| l).sum();
    let prefix: Vec<u64> = {
        let mut acc = 0u64;
        chroms
            .iter()
            .map(|(_, l)| {
                acc += l;
                acc
            })
            .collect()
    };

    let mut rng = LcgRng::new(opts.seed.unwrap_or_else(|| {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(42, |d| d.subsec_nanos().into())
    }));

    let mut bw = std::io::BufWriter::new(w);
    let mut rdr = BufReader::new(r);
    let mut line = String::new();

    loop {
        line.clear();
        if rdr.read_line(&mut line).map_err(RsomicsError::Io)? == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }

        let cols: Vec<&str> = trimmed.splitn(4, '\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let chrom = cols[0];
        let start: u64 = cols[1]
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("invalid start: {}", cols[1])))?;
        let end: u64 = cols[2]
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("invalid end: {}", cols[2])))?;
        if end < start {
            return Err(RsomicsError::InvalidInput(format!(
                "end < start for feature {chrom}:{start}-{end}"
            )));
        }
        let feat_len = end - start;

        let (target_chrom, chrom_len) = if opts.same_chrom {
            let len = genome.get(chrom).copied().ok_or_else(|| {
                RsomicsError::InvalidInput(format!("chromosome not in genome: {chrom}"))
            })?;
            (chrom, len)
        } else {
            let pick = rng.next_u64() % total_len;
            let idx = prefix.partition_point(|&p| p <= pick);
            (
                chroms[idx.min(chroms.len() - 1)].0,
                chroms[idx.min(chroms.len() - 1)].1,
            )
        };

        if chrom_len < feat_len {
            continue;
        }

        let max_start = chrom_len - feat_len;
        let new_start = rng.next_u64() % (max_start + 1);
        let new_end = new_start + feat_len;

        let extra = if cols.len() > 3 {
            format!("\t{}", cols[3])
        } else {
            String::new()
        };
        writeln!(bw, "{target_chrom}\t{new_start}\t{new_end}{extra}").map_err(RsomicsError::Io)?;
    }

    bw.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

/// Same constants as bed-sample (Knuth's MMIX).
struct LcgRng {
    state: u64,
}

impl LcgRng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn genome() -> HashMap<String, u64> {
        let mut m = HashMap::new();
        m.insert("chr1".to_owned(), 1000);
        m.insert("chr2".to_owned(), 500);
        m
    }

    fn run_shuffle(bed: &str, same_chrom: bool, seed: u64) -> Vec<String> {
        let g = genome();
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
    fn output_count_equals_input() {
        let bed = "chr1\t0\t100\nchr1\t200\t300\nchr2\t0\t50\n";
        let out = run_shuffle(bed, false, 1);
        assert_eq!(out.len(), 3, "one output per input record");
    }

    #[test]
    fn feature_length_preserved() {
        let bed = "chr1\t0\t100\nchr2\t50\t150\n";
        let out = run_shuffle(bed, false, 42);
        for line in &out {
            let cols: Vec<&str> = line.split('\t').collect();
            let end: u64 = cols[2].parse().unwrap();
            let start: u64 = cols[1].parse().unwrap();
            let len = end - start;
            assert!(len == 100, "feature length must be preserved (got {len})");
        }
    }

    #[test]
    fn same_chrom_stays_on_chrom() {
        let bed = "chr1\t0\t100\nchr1\t200\t300\n";
        let out = run_shuffle(bed, true, 7);
        for line in &out {
            assert!(
                line.starts_with("chr1\t"),
                "same_chrom must keep chromosome"
            );
        }
    }

    #[test]
    fn coords_within_chrom_bounds() {
        let bed = "chr1\t0\t100\n";
        // Run 100 times with different seeds to cover range.
        let g = genome();
        for seed in 0..100u64 {
            let mut out = Vec::new();
            shuffle(
                Cursor::new(bed),
                &mut out,
                &g,
                &ShuffleOptions {
                    same_chrom: false,
                    allow_beyond_chrom_end: false,
                    seed: Some(seed),
                },
            )
            .unwrap();
            let s = String::from_utf8(out).unwrap();
            let line = s.trim();
            let cols: Vec<&str> = line.split('\t').collect();
            let chrom = cols[0];
            let start: u64 = cols[1].parse().unwrap();
            let end: u64 = cols[2].parse().unwrap();
            let chrom_len = *g.get(chrom).unwrap();
            assert!(
                end <= chrom_len,
                "end {end} > chrom_len {chrom_len} for seed {seed}"
            );
            assert!(start <= end, "start > end for seed {seed}");
        }
    }

    #[test]
    fn deterministic_with_seed() {
        let bed = "chr1\t0\t100\nchr2\t0\t50\n";
        let a = run_shuffle(bed, false, 999);
        let b = run_shuffle(bed, false, 999);
        assert_eq!(a, b, "same seed must produce same output");
    }

    #[test]
    fn skips_headers_and_blanks() {
        let bed = "# comment\ntrack name=test\n\nchr1\t0\t100\n";
        let out = run_shuffle(bed, false, 1);
        assert_eq!(out.len(), 1, "only data records output");
    }

    #[test]
    fn extra_columns_preserved() {
        let bed = "chr1\t0\t100\tname\tscore\n";
        let out = run_shuffle(bed, true, 5);
        assert_eq!(out.len(), 1);
        let cols: Vec<&str> = out[0].split('\t').collect();
        // Extra cols after col 3 are preserved as a single string.
        assert!(cols.len() >= 4, "extra columns preserved");
    }
}
