# Comparative Genomic Analysis of Virulent Acinetobacter Phages

This repository contains the Python scripts used in:

"Genomic Landscape of 340 Virulent Acinetobacter Bacteriophages Reveals Anti-CRISPR–Enriched Candidates for Therapeutic Prioritization"

The scripts reproduce analyses reported in the manuscript, including:

1. Whole-genome phylogenetic reconstruction using k-mer distances (k=13) and UPGMA clustering
2. Pangenome frequency analysis
3. Protein frequency distribution analysis
4. Anti-CRISPR candidate identification using a guilt-by-association approach
5. CRISPR spacer match visualization

## Requirements

Python 3.9+

Required packages:

pandas
matplotlib

## Scripts

### kmer_upgma_phylogeny.py
Builds pairwise k-mer genomic distances and constructs a UPGMA phylogenetic tree.

Input:
- genome FASTA file

Output:
- distances.tsv
- phage_tree_final.nwk

### anti_crispr_candidate_finder.py
Identifies putative Anti-CRISPR candidates based on annotation keywords and genomic neighborhood proximity to regulatory proteins.

Input:
- phage_database_master.faa

Output:
- acr_candidate_report.csv

### protein_frequency_analysis.py
Counts protein annotation frequencies across phage genomes.

Input:
- phage_database_master.faa

Output:
- all_protein_frequency.csv

### pangenome_frequency_distribution.py
Generates pangenome frequency distribution data.

Input:
- true_protein_frequency.csv

Output:
- pangenome_curve_data.csv

### crispr_top10_phage_plot.py
Plots top phages with CRISPR spacer matches.

Input:
- top10_phages.csv

Output:
- Top10_CRISPR_Phage_Matches.png
