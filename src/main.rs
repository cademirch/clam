mod utils;

use clap::{Args, Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(name = "clam")]
#[command(about = "Callable Loci and More")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Parser)]
enum Commands {
    Loci {
        file: String
    },
    Stat {
        file: String
    }
}

fn main() {
    let args = Cli::parse();

    println!("Hello, world!");
}
