use bio::io::fasta;
use rayon::prelude::*;
use std::collections::HashMap;

// Function to count occurrences of a pattern in a text
fn pattern_count(text: &str, pattern: &str) -> usize {
    text.match_indices(pattern).count()
}

// Generate all k-mers of length `k` with an optional constraint on the number of 'G's
fn get_kmers(k: usize, number_of_g: Option<usize>) -> Vec<String> {
    let nucleotides = ['A', 'T', 'C', 'G'];
    let kmers: Vec<String> = (0..4_usize.pow(k as u32))
        .map(|i| {
            (0..k)
                .map(|j| nucleotides[(i >> (2 * j)) & 3])
                .rev()
                .collect()
        })
        .collect();

    if let Some(g_count) = number_of_g {
        kmers.into_iter().filter(|kmer| kmer.matches('G').count() == g_count).collect()
    } else {
        kmers
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load FASTA file
    let fasta_reader = fasta::Reader::from_file("/home/ben/genomes/hg38/genome.fasta")?;
    let mut genome_sequence = String::new();
    for result in fasta_reader.records() {
        let record = result?;
        if record.id() == "chr20" {
            genome_sequence = String::from_utf8(record.seq().to_vec())?;
            break;
        }
    }

    if genome_sequence.is_empty() {
        eprintln!("Chromosome 'chr20' not found in the FASTA file.");
        return Ok(());
    }

    // Generate k-mers with 2 'G's
    let kmers = get_kmers(4, Some(2));

    // Count occurrences for each k-mer
    let kmer_count: HashMap<String, usize> = kmers
        .par_iter()
        .map(|kmer| {
            println!("Searching for k-mer: {}", kmer);
            let count = pattern_count(&genome_sequence, kmer);
            (kmer.clone(), count)
        })
        .collect();

    // Print results
    for (kmer, count) in kmer_count.iter() {
        println!("{}: {}", kmer, count);
    }

    Ok(())
}
