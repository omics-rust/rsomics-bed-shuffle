use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use rsomics_bed_shuffle::{ShuffleOptions, parse_genome, shuffle};
use std::io::Cursor;
use std::path::Path;

fn bench_shuffle(c: &mut Criterion) {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let gf = std::fs::File::open(manifest.join("tests/golden/genome.sizes")).unwrap();
    let genome = parse_genome(gf).unwrap();

    let bed_line = "chr1\t100\t200\nchr1\t500\t600\nchr2\t0\t100\nchr3\t50\t150\n";
    let bed: String = bed_line.repeat(250);
    let bed_bytes = bed.len() as u64;

    let mut group = c.benchmark_group("bed-shuffle");
    group.throughput(Throughput::Bytes(bed_bytes));
    group.bench_function("shuffle_1000_records", |b| {
        b.iter(|| {
            let mut out = Vec::new();
            shuffle(
                Cursor::new(bed.as_bytes()),
                &mut out,
                &genome,
                &ShuffleOptions {
                    same_chrom: false,
                    allow_beyond_chrom_end: false,
                    seed: Some(42),
                },
            )
            .unwrap();
        });
    });
    group.finish();
}

criterion_group!(benches, bench_shuffle);
criterion_main!(benches);
