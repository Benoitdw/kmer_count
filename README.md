# COUNT_KMER

```bash
Program to compute the number of kmers apparing in a fasta

Usage: count_kmer [OPTIONS] --fasta <FASTA> --export <EXPORT>

Options:
-f, --fasta <FASTA> Path the the fasta
-e, --export <EXPORT> Export path
-g, --gtime <GTIME> Fix number of g in the kmer (can be ignore)
-s, --size <SIZE> len of the kmer [default: 4]
-a, --alt inclue alternate contig (default True)
-h, --help Print help
-V, --version Print version
```

## EXAMPLES

```bash
./count_kmer -f /path_to_genome/genome.fasta -e /path_to_export/kmer_count_6N2G.csv -g 2 -s 6
```

Will compute and export all the 6mers count with a g fixed number of 2
