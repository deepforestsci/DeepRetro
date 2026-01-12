# Molecular Sampling Based on Circular Fingerprint Clustering

This repository contains a Python script that samples 250 molecules from the USPTO-50k dataset based on clusters formed using RDKit circular fingerprints (Morgan fingerprints). The goal is to ensure diverse molecular sampling by clustering molecules with similar structural features and then sampling representatively from each cluster.

## Features

- **Morgan Fingerprint Generation**: Uses RDKit to generate circular fingerprints for molecular structures
- **Optimal Clustering**: Automatically finds the optimal number of clusters using silhouette score analysis
- **Multiple Sampling Strategies**: 
  - `proportional`: Sample proportionally to cluster size
  - `uniform`: Sample equally from each cluster
  - `weighted`: Inverse weighting (favor smaller clusters)
- **Visualization**: Generate plots showing cluster analysis and distribution
- **Comprehensive Logging**: Detailed logging of the entire process
- **Reproducible Results**: Set random seeds for consistent output

## Usage

### Command Line Interface

#### Basic Usage
```bash
python sample_molecules_by_clusters.py
```

#### Advanced Usage with Options
```bash
python sample_molecules_by_clusters.py \
    --input data/USPTO-50k_500.csv \
    --output my_sampled_molecules.csv \
    --n_samples 250 \
    --strategy proportional \
    --radius 2 \
    --n_bits 2048 \
    --n_clusters 20 \
    --mol_column input \
    --seed 42
```

#### Command Line Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--input` | Input CSV file path | `data/uspto_50k_test.csv` |
| `--output` | Output CSV file path | `sampled_molecules.csv` |
| `--n_samples` | Number of molecules to sample | `250` |
| `--strategy` | Sampling strategy (`proportional`, `uniform`, `weighted`) | `proportional` |
| `--radius` | Morgan fingerprint radius | `2` |
| `--n_bits` | Fingerprint size in bits | `2048` |
| `--n_clusters` | Number of clusters (auto if not specified) | `None` |
| `--mol_column` | Column containing SMILES strings | `input` |
| `--seed` | Random seed for reproducibility | `21` |




This will create three different samples using different strategies and compare their diversity.

## Input Data Format

The input CSV file should contain at least:
- A column with SMILES strings (default: `input`)
- Optional: `reaction_type` column for additional analysis
- Optional: Other molecular properties

Example:
```csv
mol_no,input,output,reaction_type
1,CCO,CC=O,8
2,CC(=O)O,CCO,7
...
```

## Output

The script generates:

1. **Sampled CSV file**: Contains the selected molecules with an additional `cluster_id` column
2. **Cluster analysis plot**: Shows cluster size distribution and silhouette score analysis
3. **Console output**: Detailed logging of the process and summary statistics

### Output CSV Format
```csv
mol_no,input,output,reaction_type,cluster_id
3,Cc1cc(C[C@@H]...),CS(=O)(=O)N1CCC...,2,0
14,COc1ccc(F)c(F)c1C...,COc1ccc(F)c(F)c1C...,6,1
...
```

## Sampling Strategies

### Proportional (Default)
- Samples molecules proportionally to cluster size
- Larger clusters contribute more molecules
- Maintains original distribution

### Uniform
- Samples equal number of molecules from each cluster
- Ensures all clusters are represented equally
- May oversample from small clusters

### Weighted
- Uses inverse weighting (favors smaller clusters)
- Helps discover rare molecular patterns
- Balances between cluster sizes

## Algorithm Overview

1. **Load Data**: Read CSV file containing molecular SMILES
2. **Generate Fingerprints**: Create Morgan fingerprints for each valid molecule
3. **Find Optimal Clusters**: Use silhouette score to determine best cluster count
4. **Perform Clustering**: Apply K-means clustering to fingerprints
5. **Sample Molecules**: Select molecules from clusters based on chosen strategy
6. **Save Results**: Export sampled molecules with cluster assignments

## Requirements

- Python 3.7+
- pandas
- numpy
- rdkit-pypi
- scikit-learn
- matplotlib
- seaborn

## Notes

- Invalid SMILES strings are automatically filtered out
- The script handles cases where requested samples exceed available molecules
- Clustering is performed using Euclidean distance on binary fingerprints
- Random seeds ensure reproducible results across runs
