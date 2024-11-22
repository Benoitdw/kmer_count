

// Function to count occurrences of a pattern in a text
pub fn pattern_count(text: &str, pattern: &str) -> usize {
    text.match_indices(pattern).count()
}

// Generate all k-mers of length `k` with an optional constraint on the number of 'G's
pub fn get_kmers(k: usize, number_of_g: Option<usize>) -> Vec<String> {
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
