mod depth;
mod gc;

use clap::Parser;
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[clap(about, version, author)]
struct Opts {
    /// A level of verbosity, and can be used multiple times
    #[clap(short, long, parse(from_occurrences))]
    verbose: i32,
    #[clap(subcommand)]
    subcmd: SubCommand,
}

#[derive(Parser)]
enum SubCommand {
    GC(GC),
    Depth(Depth),
}

#[derive(Parser)]
struct GC {
    #[clap(short = 'b', long = "bam", help = "input bam file..", validator = file_path_validation)]
    bam: PathBuf,
    #[clap(short = 'f', long = "fasta", help = "input fa file..", validator = file_path_validation)]
    fa: PathBuf,
    #[clap(short = 'r', long = "region", help = "input bam file..", validator = file_path_validation)]
    region: PathBuf,
    #[clap(
        short = 'w',
        long = "window",
        help = "Set window size for GC calculation",
        default_value = "100"
    )]
    window: usize,
    #[clap(
        short = 's',
        long = "step",
        help = "Set step size for GC calculation",
        default_value = "10000"
    )]
    step: usize,
}

#[derive(Parser)]
struct Depth {
    #[clap(short, long, help = "debug")]
    debug: bool,
}

impl SubCommand {}

fn main() {
    let opts: Opts = Opts::parse();

    // Vary the output based on how many times the user used the "verbose" flag
    // (i.e. 'myprog -v -v -v' or 'myprog -vvv' vs 'myprog -v'
    match opts.verbose {
        0 => print!(""),
        1 => println!("Some verbose info"),
        2 => println!("Tons of verbose info"),
        _ => println!("Don't be ridiculous"),
    }

    // You can handle information about subcommands by requesting their matches by name
    // (as below), requesting just the name used, or both at the same time
    match opts.subcmd {
        SubCommand::GC(o) => {
            gc::run(o.bam, o.fa, o.region, o.window, o.step);
        }
        SubCommand::Depth(o) => {
            depth::run(o.debug);
        }
    }
}

fn file_path_validation(path: &str) -> Result<(), String> {
    let path = Path::new(path);
    if !path.exists() {
        Err(format!("{path:?} file doesn't exists"))
    } else if !path.is_file() {
        Err(format!("{path:?} is not a file"))
    } else {
        Ok(())
    }
}
