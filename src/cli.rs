use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_shuffle::{ShuffleOptions, parse_genome, shuffle};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-shuffle", disable_help_flag = true)]
pub struct Cli {
    /// Input BED file (default: stdin).
    #[arg(short = 'i', long = "input")]
    pub input: Option<PathBuf>,

    /// Genome (chromosome sizes) file: two-column TSV chrom<tab>size.
    #[arg(short = 'g', long = "genome", required = true)]
    pub genome: PathBuf,

    /// Only relocate within the same chromosome as the input feature.
    #[arg(long = "chrom")]
    pub same_chrom: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let gf = File::open(&self.genome).map_err(RsomicsError::Io)?;
        let genome = parse_genome(gf)?;
        let opts = ShuffleOptions {
            same_chrom: self.same_chrom,
            allow_beyond_chrom_end: false,
            seed: self.common.seed,
        };
        let mut out: Box<dyn Write> = if self.common.json {
            Box::new(io::sink())
        } else {
            Box::new(io::stdout().lock())
        };
        if let Some(ref p) = self.input {
            let f = File::open(p).map_err(RsomicsError::Io)?;
            shuffle(f, &mut out, &genome, &opts)
        } else {
            let stdin = io::stdin();
            shuffle(stdin.lock(), &mut out, &genome, &opts)
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Randomly relocate BED intervals within a genome (bedtools shuffle equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["-g <GENOME> [OPTIONS] [-i <BED>]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdin"),
                description: "Input BED file",
                why_default: None,
            },
            FlagSpec {
                short: Some('g'),
                long: "genome",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Chromosome sizes file (chrom<tab>size)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "chrom",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Keep shuffled features on the same chromosome",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Shuffle intervals to random genome positions",
            command: "rsomics-bed-shuffle -g genome.sizes -i intervals.bed",
        },
        Example {
            description: "Reproducible shuffle, same chromosome only",
            command: "rsomics-bed-shuffle -g genome.sizes --chrom --seed 42 -i intervals.bed",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
