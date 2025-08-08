#!/usr/bin/env python3
"""
Sample 250 molecules from USPTO-50k test set based on clusters formed using circular fingerprints.

This script:
1. Loads the USPTO-50k dataset
2. Generates Morgan (circular) fingerprints for input molecules
3. Performs clustering on the fingerprints
4. Samples molecules from clusters to ensure diversity
5. Saves the sampled molecules to a new CSV file
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import argparse
import logging

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class MolecularSampler:
    """Class to handle molecular sampling based on clustering."""

    def __init__(self,
                 csv_path='data/USPTO-50k_500.csv',
                 radius=2,
                 n_bits=2048):
        """
        Initialize the molecular sampler.

        Parameters
        ----------
        csv_path : str, optional
            Path to the CSV file containing molecular data.
        radius : int, optional
            Radius for Morgan fingerprints (default: 2).
        n_bits : int, optional
            Number of bits for fingerprints (default: 2048).
        """
        self.csv_path = csv_path
        self.radius = radius
        self.n_bits = n_bits
        self.df = None
        self.fingerprints = None
        self.valid_indices = None

    def load_data(self):
        """Load the dataset from CSV file."""
        logger.info(f"Loading data from {self.csv_path}")
        self.df = pd.read_csv(self.csv_path)
        logger.info(f"Loaded {len(self.df)} molecules")
        return self.df

    def generate_fingerprints(self, mol_column='input'):
        """
        Generate Morgan fingerprints for molecules.
        
        Parameters
        ----------
        mol_column : str, optional
            Column name containing SMILES strings
            
        Returns
        -------
        numpy.ndarray: Binary fingerprint matrix
        """
        logger.info("Generating Morgan fingerprints...")
        fingerprints = []
        valid_indices = []

        for idx, smiles in enumerate(self.df[mol_column]):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    # Generate Morgan fingerprint
                    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                        mol, radius=self.radius, nBits=self.n_bits)
                    # Convert to numpy array
                    fp_array = np.zeros((self.n_bits, ))
                    DataStructs.ConvertToNumpyArray(fp, fp_array)
                    fingerprints.append(fp_array)
                    valid_indices.append(idx)
                else:
                    logger.warning(f"Invalid SMILES at index {idx}: {smiles}")
            except Exception as e:
                logger.warning(
                    f"Error processing SMILES at index {idx}: {smiles}, Error: {e}"
                )

        self.fingerprints = np.array(fingerprints)
        self.valid_indices = valid_indices
        logger.info(
            f"Generated fingerprints for {len(self.fingerprints)} valid molecules"
        )

        return self.fingerprints

    def find_optimal_clusters(self, max_clusters=50, min_clusters=5):
        """
        Find optimal number of clusters using silhouette score.
        
        Parameters
        ----------
        max_clusters : int, optional
            Maximum number of clusters to test
        min_clusters : int, optional
            Minimum number of clusters to test
            
        Returns
        -------
        int: Optimal number of clusters
        """
        logger.info("Finding optimal number of clusters...")

        if len(self.fingerprints) < max_clusters:
            max_clusters = len(self.fingerprints) // 2

        silhouette_scores = []
        cluster_range = range(min_clusters,
                              min(max_clusters + 1, len(self.fingerprints)))

        for n_clusters in cluster_range:
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            cluster_labels = kmeans.fit_predict(self.fingerprints)
            score = silhouette_score(self.fingerprints, cluster_labels)
            silhouette_scores.append(score)
            logger.info(
                f"Clusters: {n_clusters}, Silhouette Score: {score:.4f}")

        optimal_clusters = list(cluster_range)[np.argmax(silhouette_scores)]
        logger.info(f"Optimal number of clusters: {optimal_clusters}")

        return optimal_clusters, silhouette_scores, list(cluster_range)

    def perform_clustering(self, n_clusters=None):
        """
        Perform K-means clustering on fingerprints.
        
        Parameters
        ----------
        n_clusters : int, optional
            Number of clusters. If None, will find optimal number.
            
        Returns
        -------
        numpy.ndarray: Cluster labels
        """
        if n_clusters is None:
            n_clusters, _, _ = self.find_optimal_clusters()

        logger.info(
            f"Performing K-means clustering with {n_clusters} clusters...")
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(self.fingerprints)

        # Log cluster distribution
        cluster_counts = Counter(cluster_labels)
        logger.info("Cluster distribution:")
        for cluster_id, count in sorted(cluster_counts.items()):
            logger.info(f"  Cluster {cluster_id}: {count} molecules")

        return cluster_labels, kmeans

    def sample_from_clusters(self,
                             cluster_labels,
                             n_samples=250,
                             strategy='proportional'):
        """
        Sample molecules from clusters.
        
        Parameters
        ----------
        cluster_labels : numpy.ndarray
            Cluster assignments
        n_samples : int, optional
            Total number of samples to select
        strategy : str, optional
            Sampling strategy ('proportional', 'uniform', or 'weighted')
            
        Returns
        -------
        pandas.DataFrame: Sampled molecules
        """
        logger.info(
            f"Sampling {n_samples} molecules using {strategy} strategy...")

        # Get cluster information
        unique_clusters = np.unique(cluster_labels)
        cluster_counts = Counter(cluster_labels)
        total_molecules = len(cluster_labels)

        if n_samples >= total_molecules:
            logger.warning(
                f"Requested {n_samples} samples, but only {total_molecules} available. Using all molecules."
            )
            sampled_indices = list(range(len(self.valid_indices)))
        else:
            sampled_indices = []

            if strategy == 'proportional':
                # Sample proportionally to cluster size
                for cluster_id in unique_clusters:
                    cluster_size = cluster_counts[cluster_id]
                    n_from_cluster = max(
                        1, int(n_samples * cluster_size / total_molecules))

                    cluster_molecule_indices = np.where(
                        cluster_labels == cluster_id)[0]

                    if len(cluster_molecule_indices) <= n_from_cluster:
                        selected = cluster_molecule_indices
                    else:
                        selected = np.random.choice(cluster_molecule_indices,
                                                    size=n_from_cluster,
                                                    replace=False)

                    sampled_indices.extend(selected)

            elif strategy == 'uniform':
                # Sample uniformly from each cluster
                samples_per_cluster = n_samples // len(unique_clusters)
                remaining_samples = n_samples % len(unique_clusters)

                for i, cluster_id in enumerate(unique_clusters):
                    cluster_molecule_indices = np.where(
                        cluster_labels == cluster_id)[0]
                    n_from_cluster = samples_per_cluster

                    # Distribute remaining samples to first few clusters
                    if i < remaining_samples:
                        n_from_cluster += 1

                    if len(cluster_molecule_indices) <= n_from_cluster:
                        selected = cluster_molecule_indices
                    else:
                        selected = np.random.choice(cluster_molecule_indices,
                                                    size=n_from_cluster,
                                                    replace=False)

                    sampled_indices.extend(selected)

            elif strategy == 'weighted':
                # Inverse weight by cluster size (favor smaller clusters)
                weights = []
                for cluster_id in unique_clusters:
                    cluster_size = cluster_counts[cluster_id]
                    weight = 1.0 / cluster_size  # Inverse weighting
                    weights.append(weight)

                weights = np.array(weights)
                weights = weights / weights.sum()  # Normalize

                for i, cluster_id in enumerate(unique_clusters):
                    cluster_molecule_indices = np.where(
                        cluster_labels == cluster_id)[0]
                    n_from_cluster = max(1, int(n_samples * weights[i]))

                    if len(cluster_molecule_indices) <= n_from_cluster:
                        selected = cluster_molecule_indices
                    else:
                        selected = np.random.choice(cluster_molecule_indices,
                                                    size=n_from_cluster,
                                                    replace=False)

                    sampled_indices.extend(selected)

            # Adjust if we have too many or too few samples
            if len(sampled_indices) > n_samples:
                sampled_indices = np.random.choice(sampled_indices,
                                                   size=n_samples,
                                                   replace=False)
            elif len(sampled_indices) < n_samples:
                # Add random samples from remaining molecules
                remaining_indices = list(
                    set(range(len(cluster_labels))) - set(sampled_indices))
                additional_needed = n_samples - len(sampled_indices)
                if len(remaining_indices) >= additional_needed:
                    additional_samples = np.random.choice(
                        remaining_indices,
                        size=additional_needed,
                        replace=False)
                    sampled_indices.extend(additional_samples)

        # Convert back to original dataframe indices
        original_indices = [self.valid_indices[i] for i in sampled_indices]
        sampled_df = self.df.iloc[original_indices].copy()

        # Add cluster information
        sampled_cluster_labels = [cluster_labels[i] for i in sampled_indices]
        sampled_df['cluster_id'] = sampled_cluster_labels

        logger.info(f"Successfully sampled {len(sampled_df)} molecules")
        return sampled_df

    def plot_cluster_analysis(self,
                              cluster_labels,
                              silhouette_scores=None,
                              cluster_range=None,
                              output_dir='./'):
        """
        Plot cluster analysis results.
        
        Parameters
        ----------
        cluster_labels : numpy.ndarray
            Cluster assignments
        silhouette_scores : list, optional
            Silhouette scores for different cluster numbers
        cluster_range : list, optional
            Range of cluster numbers tested
        output_dir : str, optional
            Directory to save plots
        """
        plt.style.use('default')

        # Plot 1: Cluster distribution
        plt.figure(figsize=(12, 5))

        plt.subplot(1, 2, 1)
        cluster_counts = Counter(cluster_labels)
        clusters, counts = zip(*sorted(cluster_counts.items()))
        plt.bar(clusters, counts, alpha=0.7, color='skyblue', edgecolor='navy')
        plt.xlabel('Cluster ID')
        plt.ylabel('Number of Molecules')
        plt.title('Cluster Size Distribution')
        plt.grid(True, alpha=0.3)

        # Plot 2: Silhouette scores (if available)
        if silhouette_scores and cluster_range:
            plt.subplot(1, 2, 2)
            plt.plot(cluster_range,
                     silhouette_scores,
                     'o-',
                     color='red',
                     linewidth=2,
                     markersize=6)
            plt.xlabel('Number of Clusters')
            plt.ylabel('Silhouette Score')
            plt.title('Silhouette Score vs Number of Clusters')
            plt.grid(True, alpha=0.3)

            # Highlight optimal point
            optimal_idx = np.argmax(silhouette_scores)
            plt.axvline(x=cluster_range[optimal_idx],
                        color='green',
                        linestyle='--',
                        alpha=0.7)
            plt.annotate(f'Optimal: {cluster_range[optimal_idx]}',
                         xy=(cluster_range[optimal_idx],
                             silhouette_scores[optimal_idx]),
                         xytext=(10, 10),
                         textcoords='offset points',
                         bbox=dict(boxstyle='round,pad=0.3',
                                   facecolor='yellow',
                                   alpha=0.7))

        plt.tight_layout()
        plt.savefig(f'{output_dir}/cluster_analysis_rad_4_n_bits_2048.png',
                    dpi=300,
                    bbox_inches='tight')
        plt.show()

        logger.info(
            f"Cluster analysis plot saved to {output_dir}/cluster_analysis_rad_4_n_bits_2048.png"
        )

    def save_results(self, sampled_df, output_path='sampled_molecules.csv'):
        """
        Save sampled molecules to CSV file.
        
        Parameters
        ----------
        sampled_df : pandas.DataFrame
            Sampled molecules
        output_path : str, optional
            Output file path
        """
        sampled_df.to_csv(output_path, index=False)
        logger.info(f"Sampled molecules saved to {output_path}")

        # Print summary statistics
        logger.info("\n=== SAMPLING SUMMARY ===")
        logger.info(f"Total molecules sampled: {len(sampled_df)}")
        logger.info(
            f"Number of clusters represented: {sampled_df['cluster_id'].nunique()}"
        )
        logger.info(f"Cluster distribution in sample:")
        cluster_dist = sampled_df['cluster_id'].value_counts().sort_index()
        for cluster_id, count in cluster_dist.items():
            logger.info(f"  Cluster {cluster_id}: {count} molecules")

        if 'reaction_type' in sampled_df.columns:
            logger.info(f"Reaction type distribution:")
            reaction_dist = sampled_df['reaction_type'].value_counts()
            for reaction_type, count in reaction_dist.items():
                logger.info(f"  Type {reaction_type}: {count} molecules")


def main():
    """Main function to run the molecular sampling pipeline."""
    parser = argparse.ArgumentParser(
        description='Sample molecules from USPTO-50k based on clustering')
    parser.add_argument('--input',
                        default='../data/uspto_50k_test.csv',
                        help='Input CSV file path')
    parser.add_argument('--output',
                        default='sampled_molecules_rad_4_n_bits_2048.csv',
                        help='Output CSV file path')
    parser.add_argument('--n_samples',
                        type=int,
                        default=250,
                        help='Number of molecules to sample')
    parser.add_argument('--strategy',
                        choices=['proportional', 'uniform', 'weighted'],
                        default='proportional',
                        help='Sampling strategy')
    parser.add_argument('--radius',
                        type=int,
                        default=2,
                        help='Morgan fingerprint radius')
    parser.add_argument('--n_bits',
                        type=int,
                        default=2048,
                        help='Fingerprint size in bits')
    parser.add_argument('--n_clusters',
                        type=int,
                        default=None,
                        help='Number of clusters (auto if not specified)')
    parser.add_argument('--mol_column',
                        default='input',
                        help='Column containing SMILES strings')
    parser.add_argument('--seed',
                        type=int,
                        default=21,
                        help='Random seed for reproducibility')

    args = parser.parse_args()

    # Set random seed for reproducibility
    np.random.seed(args.seed)

    # Initialize sampler
    sampler = MolecularSampler(csv_path=args.input,
                               radius=args.radius,
                               n_bits=args.n_bits)

    try:
        # Load data
        df = sampler.load_data()

        # Generate fingerprints
        fingerprints = sampler.generate_fingerprints(
            mol_column=args.mol_column)

        if len(fingerprints) == 0:
            logger.error(
                "No valid fingerprints generated. Check your input data.")
            return

        # Perform clustering
        if args.n_clusters:
            cluster_labels, kmeans_model = sampler.perform_clustering(
                n_clusters=args.n_clusters)
            silhouette_scores, cluster_range = None, None
        else:
            optimal_clusters, silhouette_scores, cluster_range = sampler.find_optimal_clusters(
            )
            cluster_labels, kmeans_model = sampler.perform_clustering(
                n_clusters=optimal_clusters)

        # Sample molecules
        sampled_df = sampler.sample_from_clusters(cluster_labels,
                                                  n_samples=args.n_samples,
                                                  strategy=args.strategy)

        # Plot analysis
        sampler.plot_cluster_analysis(cluster_labels, silhouette_scores,
                                      cluster_range)

        # Save results
        sampler.save_results(sampled_df, output_path=args.output)

        logger.info("Molecular sampling completed successfully!")

    except Exception as e:
        logger.error(f"Error during molecular sampling: {e}")
        raise


if __name__ == "__main__":
    main()
