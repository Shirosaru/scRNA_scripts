"""
Systematic Perturbation Analysis
================================

This script performs systematic analysis of multiple perturbations:
1. Analyzes each perturbation vs control
2. Compares perturbations to each other
3. Identifies common and unique effects
4. Creates comprehensive reports

Usage:
    python systematic_analysis.py
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import json
from datetime import datetime
from typing import List, Dict, Tuple
import random

warnings.filterwarnings('ignore')

# Set matplotlib parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.facecolor'] = 'white'
sc.settings.verbosity = 3

# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    """Configuration for systematic analysis."""
    
    DATA_FILE = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'
    OUTPUT_DIR = 'systematic_analysis_results'
    PERTURBATION_COLUMN = 'gene_target'
    CONTROL_PERTURBATION = 'Non-Targeting'
    
    # Analysis parameters
    MIN_CELLS_PER_PERTURBATION = 100  # Minimum cells to analyze a perturbation
    N_CELLS_TO_SUBSAMPLE = 200000  # Subsample size (None for all cells)
    N_TOP_PERTURBATIONS = 20  # Analyze top N perturbations (set to None for all)
    
    # QC parameters
    MIN_GENES_PER_CELL = 200
    MIN_CELLS_PER_GENE = 3
    MAX_MITO_PCT = 20.0
    
    # DE parameters
    DE_METHOD = 'wilcoxon'
    N_TOP_GENES = 100
    FDR_THRESHOLD = 0.05
    LOGFC_THRESHOLD = 0.5
    
    # Random seed
    RANDOM_SEED = 42

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def load_and_subsample(file_path: str, n_cells: int = None, random_seed: int = 42):
    """Load data and optionally subsample."""
    print(f"Loading data from: {file_path}")
    adata_backed = ad.read_h5ad(file_path, backed='r')
    print(f"Total cells: {adata_backed.n_obs:,}")
    
    if n_cells and n_cells < adata_backed.n_obs:
        print(f"Subsampling to {n_cells:,} cells...")
        random.seed(random_seed)
        np.random.seed(random_seed)
        all_indices = list(range(adata_backed.n_obs))
        subsample_indices = random.sample(all_indices, n_cells)
        adata = adata_backed[subsample_indices, :].to_memory()
    else:
        adata = adata_backed.to_memory()
    
    print(f"Loaded {adata.n_obs:,} cells")
    return adata

def preprocess_data(adata: ad.AnnData):
    """Standard preprocessing pipeline."""
    print("Preprocessing...")
    
    # QC
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Filter
    sc.pp.filter_cells(adata, min_genes=Config.MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=Config.MIN_CELLS_PER_GENE)
    if 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] < Config.MAX_MITO_PCT, :]
    
    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # HVG and PCA
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.pca(adata, use_highly_variable=True, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)
    
    print(f"After preprocessing: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    return adata

def get_perturbations_to_analyze(adata: ad.AnnData) -> List[str]:
    """Get list of perturbations to analyze."""
    perturbations = adata.obs[Config.PERTURBATION_COLUMN].value_counts()
    
    # Filter by minimum cells
    valid_perturbations = perturbations[
        perturbations >= Config.MIN_CELLS_PER_PERTURBATION
    ]
    
    # Exclude control
    valid_perturbations = valid_perturbations[
        valid_perturbations.index != Config.CONTROL_PERTURBATION
    ]
    
    # Get top N
    top_perturbations = valid_perturbations.head(Config.N_TOP_PERTURBATIONS)
    
    print(f"\nPerturbations to analyze: {len(top_perturbations)}")
    print("Top perturbations:")
    for pert, count in top_perturbations.items():
        print(f"  • {pert}: {count:,} cells")
    
    return list(top_perturbations.index)

def analyze_single_perturbation(adata: ad.AnnData, perturbation: str, 
                                output_dir: Path) -> Dict:
    """Analyze a single perturbation vs control."""
    print(f"\n{'='*60}")
    print(f"Analyzing: {perturbation}")
    print(f"{'='*60}")
    
    # Filter to perturbation and control
    pert_mask = adata.obs[Config.PERTURBATION_COLUMN].isin([perturbation, Config.CONTROL_PERTURBATION])
    adata_subset = adata[pert_mask, :].copy()
    
    print(f"  Cells: {adata_subset.n_obs:,}")
    print(f"    - {perturbation}: {(adata_subset.obs[Config.PERTURBATION_COLUMN] == perturbation).sum():,}")
    print(f"    - {Config.CONTROL_PERTURBATION}: {(adata_subset.obs[Config.PERTURBATION_COLUMN] == Config.CONTROL_PERTURBATION).sum():,}")
    
    # Create comparison group
    adata_subset.obs['comparison'] = adata_subset.obs[Config.PERTURBATION_COLUMN].copy()
    
    # Differential expression
    print("  Computing differential expression...")
    sc.tl.rank_genes_groups(
        adata_subset,
        'comparison',
        method=Config.DE_METHOD,
        groups=[perturbation],
        reference=Config.CONTROL_PERTURBATION,
        n_genes=Config.N_TOP_GENES
    )
    
    # Extract results
    de_results = pd.DataFrame({
        'gene': adata_subset.uns['rank_genes_groups']['names'][perturbation],
        'logfoldchange': adata_subset.uns['rank_genes_groups']['logfoldchanges'][perturbation],
        'pval': adata_subset.uns['rank_genes_groups']['pvals'][perturbation],
        'pval_adj': adata_subset.uns['rank_genes_groups']['pvals_adj'][perturbation],
        'scores': adata_subset.uns['rank_genes_groups']['scores'][perturbation]
    })
    
    # Filter significant genes
    significant = de_results[
        (de_results['pval_adj'] < Config.FDR_THRESHOLD) &
        (np.abs(de_results['logfoldchange']) > Config.LOGFC_THRESHOLD)
    ]
    
    n_up = (significant['logfoldchange'] > 0).sum()
    n_down = (significant['logfoldchange'] < 0).sum()
    
    print(f"  Significant genes: {len(significant)} (↑{n_up}, ↓{n_down})")
    
    # Save results
    output_file = output_dir / f'DE_{perturbation}.csv'
    de_results.to_csv(output_file, index=False)
    
    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Volcano plot
    de_results['-log10_pval'] = -np.log10(de_results['pval_adj'] + 1e-300)
    de_results['significant'] = (
        (de_results['pval_adj'] < Config.FDR_THRESHOLD) &
        (np.abs(de_results['logfoldchange']) > Config.LOGFC_THRESHOLD)
    )
    
    axes[0].scatter(
        de_results['logfoldchange'],
        de_results['-log10_pval'],
        c=de_results['significant'],
        cmap='RdYlBu_r',
        alpha=0.6,
        s=20
    )
    axes[0].axhline(-np.log10(Config.FDR_THRESHOLD), color='red', linestyle='--', alpha=0.5)
    axes[0].axvline(-Config.LOGFC_THRESHOLD, color='gray', linestyle='--', alpha=0.5)
    axes[0].axvline(Config.LOGFC_THRESHOLD, color='gray', linestyle='--', alpha=0.5)
    axes[0].set_xlabel('Log2 Fold Change')
    axes[0].set_ylabel('-Log10 Adjusted P-value')
    axes[0].set_title(f'Volcano Plot: {perturbation}')
    axes[0].grid(alpha=0.3)
    
    # Top genes bar plot
    top_20 = de_results.head(20)
    colors = ['red' if x > 0 else 'blue' for x in top_20['logfoldchange']]
    axes[1].barh(range(len(top_20)), top_20['logfoldchange'], color=colors)
    axes[1].set_yticks(range(len(top_20)))
    axes[1].set_yticklabels(top_20['gene'], fontsize=8)
    axes[1].set_xlabel('Log2 Fold Change')
    axes[1].set_title('Top 20 DE Genes')
    axes[1].axvline(0, color='black', linestyle='-', alpha=0.3)
    axes[1].grid(alpha=0.3, axis='x')
    
    plt.tight_layout()
    plot_file = output_dir / f'volcano_{perturbation}.pdf'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        'perturbation': perturbation,
        'n_cells': adata_subset.n_obs,
        'n_control_cells': (adata_subset.obs[Config.PERTURBATION_COLUMN] == Config.CONTROL_PERTURBATION).sum(),
        'n_pert_cells': (adata_subset.obs[Config.PERTURBATION_COLUMN] == perturbation).sum(),
        'n_significant_genes': len(significant),
        'n_upregulated': n_up,
        'n_downregulated': n_down,
        'top_upregulated': significant[significant['logfoldchange'] > 0].head(10)['gene'].tolist(),
        'top_downregulated': significant[significant['logfoldchange'] < 0].head(10)['gene'].tolist(),
        'de_file': str(output_file),
        'plot_file': str(plot_file)
    }

def create_summary_comparison(all_results: List[Dict], output_dir: Path):
    """Create summary comparison across all perturbations."""
    print(f"\n{'='*60}")
    print("Creating summary comparison")
    print(f"{'='*60}")
    
    # Create summary dataframe
    summary_data = []
    for result in all_results:
        summary_data.append({
            'Perturbation': result['perturbation'],
            'N_Cells': result['n_pert_cells'],
            'N_Significant_Genes': result['n_significant_genes'],
            'N_Upregulated': result['n_upregulated'],
            'N_Downregulated': result['n_downregulated']
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('N_Significant_Genes', ascending=False)
    
    # Save summary
    summary_file = output_dir / 'summary_all_perturbations.csv'
    summary_df.to_csv(summary_file, index=False)
    print(f"Summary saved to: {summary_file}")
    
    # Create visualizations
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Number of significant genes per perturbation
    axes[0, 0].barh(range(len(summary_df)), summary_df['N_Significant_Genes'])
    axes[0, 0].set_yticks(range(len(summary_df)))
    axes[0, 0].set_yticklabels(summary_df['Perturbation'], fontsize=8)
    axes[0, 0].set_xlabel('Number of Significant Genes')
    axes[0, 0].set_title('Significant Genes per Perturbation')
    axes[0, 0].grid(alpha=0.3, axis='x')
    axes[0, 0].invert_yaxis()
    
    # 2. Up vs Down regulated
    x_pos = np.arange(len(summary_df))
    width = 0.35
    axes[0, 1].barh(x_pos - width/2, summary_df['N_Upregulated'], width, 
                    label='Upregulated', color='red', alpha=0.7)
    axes[0, 1].barh(x_pos + width/2, summary_df['N_Downregulated'], width,
                    label='Downregulated', color='blue', alpha=0.7)
    axes[0, 1].set_yticks(x_pos)
    axes[0, 1].set_yticklabels(summary_df['Perturbation'], fontsize=8)
    axes[0, 1].set_xlabel('Number of Genes')
    axes[0, 1].set_title('Up vs Down Regulated Genes')
    axes[0, 1].legend()
    axes[0, 1].grid(alpha=0.3, axis='x')
    axes[0, 1].invert_yaxis()
    
    # 3. Scatter: cells vs significant genes
    axes[1, 0].scatter(summary_df['N_Cells'], summary_df['N_Significant_Genes'], 
                       s=100, alpha=0.6)
    for idx, row in summary_df.iterrows():
        axes[1, 0].annotate(row['Perturbation'], 
                           (row['N_Cells'], row['N_Significant_Genes']),
                           fontsize=7, alpha=0.7)
    axes[1, 0].set_xlabel('Number of Cells')
    axes[1, 0].set_ylabel('Number of Significant Genes')
    axes[1, 0].set_title('Cells vs Significant Genes')
    axes[1, 0].grid(alpha=0.3)
    
    # 4. Distribution of significant genes
    axes[1, 1].hist(summary_df['N_Significant_Genes'], bins=20, edgecolor='black', alpha=0.7)
    axes[1, 1].axvline(summary_df['N_Significant_Genes'].mean(), 
                       color='red', linestyle='--', 
                       label=f'Mean: {summary_df["N_Significant_Genes"].mean():.1f}')
    axes[1, 1].set_xlabel('Number of Significant Genes')
    axes[1, 1].set_ylabel('Number of Perturbations')
    axes[1, 1].set_title('Distribution of Significant Genes')
    axes[1, 1].legend()
    axes[1, 1].grid(alpha=0.3)
    
    plt.tight_layout()
    comparison_file = output_dir / 'summary_comparison.pdf'
    plt.savefig(comparison_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Comparison plot saved to: {comparison_file}")
    
    return summary_df

def find_common_genes(all_results: List[Dict], output_dir: Path):
    """Find genes commonly affected across perturbations."""
    print(f"\n{'='*60}")
    print("Finding common affected genes")
    print(f"{'='*60}")
    
    # Collect all significant genes per perturbation
    pert_to_genes = {}
    for result in all_results:
        pert = result['perturbation']
        # Load DE results
        de_df = pd.read_csv(result['de_file'])
        significant = de_df[
            (de_df['pval_adj'] < Config.FDR_THRESHOLD) &
            (np.abs(de_df['logfoldchange']) > Config.LOGFC_THRESHOLD)
        ]
        pert_to_genes[pert] = set(significant['gene'].tolist())
    
    # Find genes in multiple perturbations
    all_genes = set()
    for genes in pert_to_genes.values():
        all_genes.update(genes)
    
    gene_counts = {}
    for gene in all_genes:
        count = sum(1 for genes in pert_to_genes.values() if gene in genes)
        gene_counts[gene] = count
    
    # Sort by frequency
    common_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)
    
    # Save
    common_df = pd.DataFrame(common_genes, columns=['Gene', 'N_Perturbations'])
    common_file = output_dir / 'common_affected_genes.csv'
    common_df.to_csv(common_file, index=False)
    
    print(f"\nTop 20 most commonly affected genes:")
    for gene, count in common_genes[:20]:
        print(f"  • {gene}: affected in {count} perturbations")
    
    print(f"\nCommon genes saved to: {common_file}")
    
    return common_df

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Main systematic analysis pipeline."""
    print("=" * 80)
    print("SYSTEMATIC PERTURBATION ANALYSIS")
    print("=" * 80)
    
    # Set random seed
    random.seed(Config.RANDOM_SEED)
    np.random.seed(Config.RANDOM_SEED)
    
    # Create output directory
    output_dir = Path(Config.OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load and preprocess data
    adata = load_and_subsample(Config.DATA_FILE, Config.N_CELLS_TO_SUBSAMPLE, Config.RANDOM_SEED)
    adata = preprocess_data(adata)
    
    # Get perturbations to analyze
    perturbations = get_perturbations_to_analyze(adata)
    
    if not perturbations:
        print("No perturbations found to analyze!")
        return
    
    # Analyze each perturbation
    all_results = []
    for pert in perturbations:
        try:
            result = analyze_single_perturbation(adata, pert, output_dir)
            all_results.append(result)
        except Exception as e:
            print(f"Error analyzing {pert}: {e}")
            continue
    
    # Create summary comparison
    summary_df = create_summary_comparison(all_results, output_dir)
    
    # Find common genes
    common_genes_df = find_common_genes(all_results, output_dir)
    
    # Save final summary
    final_summary = {
        'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'n_perturbations_analyzed': len(all_results),
        'control_perturbation': Config.CONTROL_PERTURBATION,
        'summary_statistics': {
            'mean_significant_genes': float(summary_df['N_Significant_Genes'].mean()),
            'median_significant_genes': float(summary_df['N_Significant_Genes'].median()),
            'total_unique_genes_affected': len(common_genes_df),
        },
        'perturbations': [r['perturbation'] for r in all_results]
    }
    
    summary_json = output_dir / 'analysis_summary.json'
    with open(summary_json, 'w') as f:
        json.dump(final_summary, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_dir.absolute()}")
    print(f"\nSummary:")
    print(f"  • Analyzed {len(all_results)} perturbations")
    print(f"  • Mean significant genes per perturbation: {summary_df['N_Significant_Genes'].mean():.1f}")
    print(f"  • Total unique genes affected: {len(common_genes_df)}")
    print(f"\nKey files:")
    print(f"  • summary_all_perturbations.csv - Overview of all perturbations")
    print(f"  • common_affected_genes.csv - Genes affected across multiple perturbations")
    print(f"  • DE_<perturbation>.csv - DE results for each perturbation")
    print(f"  • summary_comparison.pdf - Visual comparison")
    print(f"{'='*80}")

if __name__ == '__main__':
    main()

