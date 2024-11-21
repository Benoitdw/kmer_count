use bio::io::fasta;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;


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
    let mut kmer_count: HashMap<String, HashMap<String, usize>> = HashMap::new();
    let kmers = get_kmers(4, Some(2));
    for record in fasta_reader.records() {
        let seq = record? ;
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

    let mut writer = csv::Writer::from_writer(File::create("kmer_count_4N.tsv")?);

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
