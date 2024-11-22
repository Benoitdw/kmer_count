use clap::{error::ErrorKind, CommandFactory, Parser};
use std::path::PathBuf;
use std::collections::HashMap;
use std::fs::File;
use rayon::prelude::*;
use count_kmer::{pattern_count, get_kmers};
use bio::io::fasta;



/// Program to compute the number of kmers apparing in a fasta
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path the the fasta
    #[arg(short, long)]
    fasta: PathBuf,

    /// Export path
    #[arg(short, long)]
    export: PathBuf,

    /// Fix number of g in the kmer (can be ignore)
    #[arg(short, long)]
    gtime: Option<usize>,

    /// len of the kmer
    #[arg(short, long, default_value_t = 4)]
    size: usize,

    ///inclue alternate contig (default True)
    #[arg(short, long, action, default_value_t=false)]
    alt : bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    if !args.fasta.exists() {
        let mut cmd = Args::command();
        cmd.error(
            ErrorKind::ValueValidation,
            format!(
                "input fasta `{}` doesn't exist",
                args.fasta.display()
            ),
        )
        .exit();
    } 
    let fasta_reader = fasta::Reader::from_file(args.fasta)?;
    // TODO
        // Load FASTA file
        
    let mut kmer_count: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let kmers = get_kmers(args.size, args.gtime);
    for record in fasta_reader.records() {
        let seq = record? ;
        if !args.alt && seq.id().contains('_') {
            continue;
        }
        let genome_sequence = String::from_utf8(seq.seq().to_vec())?;
        let kmer_count_in_contig: HashMap<String, usize> = kmers
        .par_iter()
        .map(|kmer| {
            println!("Searching for k-mer: {} in {}", kmer, seq.id());
            let count = pattern_count(&genome_sequence, kmer);
            (kmer.clone(), count)
        })
        .collect();
        kmer_count.insert(seq.id().to_string(), kmer_count_in_contig);
    }
    
    let mut writer = csv::Writer::from_writer(File::create(args.export)?);

    // Write the header row
    let mut header = vec!["Nucleotide"];
    for chr in kmer_count.keys() {
        header.push(chr);
    }
    writer.write_record(&header)?;

    // Write the data rows
    let mut all_nucleotides: Vec<_> = kmer_count
        .values()
        .flat_map(|inner_map| inner_map.keys().cloned())
        .collect();
    all_nucleotides.sort();
    all_nucleotides.dedup();

    for nucleotide in all_nucleotides {
        let mut row = vec![nucleotide.clone()];
        for (_chr, inner_map) in kmer_count.iter() {
            row.push(inner_map.get(&nucleotide).map(|count| count.to_string()).unwrap_or_else(|| "0".to_string()));
        }
        writer.write_record(&row)?;
    }
    
        Ok(())
    }
    