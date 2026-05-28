# rsomics-bed-shuffle

Randomly relocate BED intervals within a genome — `bedtools shuffle` equivalent.

Reads BED intervals and relocates each to a random genomic position, preserving
feature length. The genome file is a two-column TSV of chromosome name and size.

## Usage

```
rsomics-bed-shuffle -g <GENOME> [BED]
```

## Origin

Independent Rust reimplementation of `bedtools shuffle` based on:
- The bedtools documentation and man page
- Black-box behaviour testing against bedtools 2.31.1

No GPL source was used as reference.

License: MIT OR Apache-2.0
Upstream credit: bedtools <https://bedtools.readthedocs.io/> (MIT)
