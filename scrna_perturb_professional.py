"""
Professional Single-Cell RNA Perturbation Analysis Pipeline
============================================================

This script provides a comprehensive analysis pipeline for single-cell RNA-seq
perturbation data, including quality control, dimensionality reduction, differential
expression analysis, and pathway enrichment.

Author: Professional Analysis Pipeline
Version: 1.0.0
Date: 2024
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import logging
import os
import json
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import random

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set scanpy settings for better reproducibility
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=300, facecolor='white', format='pdf')

# ============================================================================
# CONFIGURATION
# ============================================================================

class Config:
    """Configuration class for analysis parameters."""
    
    # Data paths
    DATA_FILE = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'
    OUTPUT_DIR = 'perturbation_analysis_results'
    
    # Perturbation settings
    PERTURBATION_COLUMN = 'gene_target'  # Column name for perturbation labels
    TARGET_PERTURBATION = 'TP53'
    CONTROL_PERTURBATION = 'AAV-control'
    
    # Subsampling (set to None to use all cells)
    N_CELLS_TO_SUBSAMPLE = 100000
    
    # Quality control parameters
    MIN_GENES_PER_CELL = 200
    MIN_CELLS_PER_GENE = 3
    MAX_MITO_PCT = 20.0  # Maximum mitochondrial percentage
    
    # Normalization
    TARGET_SUM = 1e4
    
    # Highly variable genes
    HVG_MIN_MEAN = 0.0125
    HVG_MAX_MEAN = 3
    HVG_MIN_DISP = 0.5
    
    # Dimensionality reduction
    N_NEIGHBORS = 15
    N_PCS = 50
    RESOLUTION = 0.5  # For clustering
    
    # Differential expression
    DE_METHOD = 'wilcoxon'  # 'wilcoxon', 't-test', or 'logreg'
    N_TOP_GENES = 50
    
    # Visualization
    FIGURE_FORMAT = 'pdf'
    DPI = 300
    
    # Random seed for reproducibility
    RANDOM_SEED = 42

# ============================================================================
# LOGGING SETUP
# ============================================================================

def setup_logging(output_dir: str) -> logging.Logger:
    """Set up logging to both file and console."""
    log_dir = Path(output_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = log_dir / f"analysis_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    logger = logging.getLogger(__name__)
    logger.info("=" * 80)
    logger.info("Single-Cell RNA Perturbation Analysis Pipeline")
    logger.info("=" * 80)
    logger.info(f"Log file: {log_file}")
    
    return logger

# ============================================================================
# DATA LOADING
# ============================================================================

def load_data(file_path: str, n_cells: Optional[int] = None, 
              logger: Optional[logging.Logger] = None) -> ad.AnnData:
    """
    Load AnnData object, optionally subsampling cells.
    
    Parameters:
    -----------
    file_path : str
        Path to .h5ad file
    n_cells : int, optional
        Number of cells to randomly subsample
    logger : logging.Logger, optional
        Logger instance
        
    Returns:
    --------
    ad.AnnData
        Loaded AnnData object
    """
    if logger:
        logger.info(f"Loading data from: {file_path}")
    
    # Load in backed mode for memory efficiency
    adata_backed = ad.read_h5ad(file_path, backed='r')
    
    if logger:
        logger.info(f"Total cells in dataset: {adata_backed.n_obs:,}")
        logger.info(f"Total genes in dataset: {adata_backed.n_vars:,}")
        logger.info(f"Available metadata columns: {list(adata_backed.obs.columns)}")
    
    # Check if perturbation column exists
    if Config.PERTURBATION_COLUMN not in adata_backed.obs.columns:
        raise ValueError(f"Perturbation column '{Config.PERTURBATION_COLUMN}' not found in data!")
    
    # Check if target and control exist
    if Config.TARGET_PERTURBATION not in adata_backed.obs[Config.PERTURBATION_COLUMN].cat.categories:
        raise ValueError(f"Target perturbation '{Config.TARGET_PERTURBATION}' not found!")
    if Config.CONTROL_PERTURBATION not in adata_backed.obs[Config.PERTURBATION_COLUMN].cat.categories:
        raise ValueError(f"Control perturbation '{Config.CONTROL_PERTURBATION}' not found!")
    
    # Subsample if requested
    if n_cells and n_cells < adata_backed.n_obs:
        if logger:
            logger.info(f"Subsampling to {n_cells:,} cells...")
        random.seed(Config.RANDOM_SEED)
        np.random.seed(Config.RANDOM_SEED)
        all_indices = list(range(adata_backed.n_obs))
        subsample_indices = random.sample(all_indices, n_cells)
        adata = adata_backed[subsample_indices, :].to_memory()
    else:
        adata = adata_backed.to_memory()
    
    if logger:
        logger.info(f"Loaded {adata.n_obs:,} cells and {adata.n_vars:,} genes")
    
    return adata

# ============================================================================
# QUALITY CONTROL
# ============================================================================

def calculate_qc_metrics(adata: ad.AnnData, logger: Optional[logging.Logger] = None) -> ad.AnnData:
    """
    Calculate quality control metrics.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    logger : logging.Logger, optional
        Logger instance
        
    Returns:
    --------
    ad.AnnData
        AnnData with QC metrics added
    """
    if logger:
        logger.info("Calculating quality control metrics...")
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    if logger:
        logger.info(f"Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.1f}")
        logger.info(f"Mean counts per cell: {adata.obs['total_counts'].mean():.1f}")
        logger.info(f"Mean mitochondrial %: {adata.obs['pct_counts_mt'].mean():.2f}%")
    
    return adata

def filter_data(adata: ad.AnnData, logger: Optional[logging.Logger] = None) -> ad.AnnData:
    """
    Filter cells and genes based on QC metrics.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    logger : logging.Logger, optional
        Logger instance
        
    Returns:
    --------
    ad.AnnData
        Filtered AnnData object
    """
    if logger:
        logger.info("Filtering cells and genes...")
        logger.info(f"Before filtering: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=Config.MIN_GENES_PER_CELL)
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=Config.MIN_CELLS_PER_GENE)
    
    # Filter by mitochondrial percentage
    if 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] < Config.MAX_MITO_PCT, :]
    
    if logger:
        logger.info(f"After filtering: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    return adata

def plot_qc_metrics(adata: ad.AnnData, output_dir: str, 
                   logger: Optional[logging.Logger] = None):
    """
    Create QC metric plots.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    output_dir : str
        Output directory for plots
    logger : logging.Logger, optional
        Logger instance
    """
    if logger:
        logger.info("Generating QC metric plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Total counts
    axes[0, 0].hist(adata.obs['total_counts'], bins=50, edgecolor='black')
    axes[0, 0].set_xlabel('Total counts per cell')
    axes[0, 0].set_ylabel('Number of cells')
    axes[0, 0].set_title('Distribution of Total Counts')
    axes[0, 0].axvline(adata.obs['total_counts'].mean(), color='red', 
                       linestyle='--', label=f'Mean: {adata.obs["total_counts"].mean():.0f}')
    axes[0, 0].legend()
    
    # Genes per cell
    axes[0, 1].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black')
    axes[0, 1].set_xlabel('Number of genes per cell')
    axes[0, 1].set_ylabel('Number of cells')
    axes[0, 1].set_title('Distribution of Genes per Cell')
    axes[0, 1].axvline(adata.obs['n_genes_by_counts'].mean(), color='red', 
                       linestyle='--', label=f'Mean: {adata.obs["n_genes_by_counts"].mean():.0f}')
    axes[0, 1].legend()
    
    # Mitochondrial percentage
    if 'pct_counts_mt' in adata.obs.columns:
        axes[1, 0].hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black')
        axes[1, 0].set_xlabel('Mitochondrial %')
        axes[1, 0].set_ylabel('Number of cells')
        axes[1, 0].set_title('Distribution of Mitochondrial %')
        axes[1, 0].axvline(adata.obs['pct_counts_mt'].mean(), color='red', 
                          linestyle='--', label=f'Mean: {adata.obs["pct_counts_mt"].mean():.2f}%')
        axes[1, 0].legend()
    
    # Scatter: counts vs genes
    axes[1, 1].scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], 
                      alpha=0.5, s=1)
    axes[1, 1].set_xlabel('Total counts')
    axes[1, 1].set_ylabel('Number of genes')
    axes[1, 1].set_title('Counts vs Genes per Cell')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'qc_metrics.pdf', dpi=Config.DPI, bbox_inches='tight')
    plt.close()
    
    if logger:
        logger.info("QC plots saved")

# ============================================================================
# PREPROCESSING
# ============================================================================

def preprocess_data(adata: ad.AnnData, logger: Optional[logging.Logger] = None) -> ad.AnnData:
    """
    Normalize and log-transform data.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    logger : logging.Logger, optional
        Logger instance
        
    Returns:
    --------
    ad.AnnData
        Preprocessed AnnData object
    """
    if logger:
        logger.info("Normalizing and log-transforming data...")
    
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=Config.TARGET_SUM)
    
    # Log transform
    sc.pp.log1p(adata)
    
    if logger:
        logger.info("Preprocessing complete")
    
    return adata

def find_hvg_and_reduce_dimensions(adata: ad.AnnData, 
                                   logger: Optional[logging.Logger] = None) -> ad.AnnData:
    """
    Find highly variable genes and perform dimensionality reduction.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    logger : logging.Logger, optional
        Logger instance
        
    Returns:
    --------
    ad.AnnData
        AnnData with PCA, neighbors, and UMAP computed
    """
    if logger:
        logger.info("Finding highly variable genes...")
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(
        adata, 
        min_mean=Config.HVG_MIN_MEAN, 
        max_mean=Config.HVG_MAX_MEAN, 
        min_disp=Config.HVG_MIN_DISP
    )
    
    n_hvg = adata.var['highly_variable'].sum()
    if logger:
        logger.info(f"Found {n_hvg:,} highly variable genes")
    
    # PCA
    if logger:
        logger.info("Computing PCA...")
    sc.pp.pca(adata, use_highly_variable=True, n_comps=Config.N_PCS)
    
    # Compute neighborhood graph
    if logger:
        logger.info("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=Config.N_NEIGHBORS, n_pcs=Config.N_PCS)
    
    # UMAP
    if logger:
        logger.info("Computing UMAP...")
    sc.tl.umap(adata)
    
    # Leiden clustering
    if logger:
        logger.info("Performing Leiden clustering...")
    sc.tl.leiden(adata, resolution=Config.RESOLUTION, key_added='leiden')
    
    if logger:
        logger.info("Dimensionality reduction complete")
    
    return adata

# ============================================================================
# PERTURBATION ANALYSIS
# ============================================================================

def prepare_perturbation_labels(adata: ad.AnnData, 
                                logger: Optional[logging.Logger] = None) -> ad.AnnData:
    """
    Prepare perturbation labels for analysis.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    logger : logging.Logger, optional
        Logger instance
        
    Returns:
    --------
    ad.AnnData
        AnnData with perturbation labels added
    """
    if logger:
        logger.info("Preparing perturbation labels...")
    
    # Create comparison groups
    adata.obs['perturbation_group'] = 'other'
    adata.obs.loc[
        adata.obs[Config.PERTURBATION_COLUMN] == Config.TARGET_PERTURBATION, 
        'perturbation_group'
    ] = Config.TARGET_PERTURBATION
    adata.obs.loc[
        adata.obs[Config.PERTURBATION_COLUMN] == Config.CONTROL_PERTURBATION, 
        'perturbation_group'
    ] = Config.CONTROL_PERTURBATION
    
    # Count cells per group
    group_counts = adata.obs['perturbation_group'].value_counts()
    if logger:
        logger.info(f"Cell counts per group:\n{group_counts}")
    
    return adata

def plot_perturbation_umap(adata: ad.AnnData, output_dir: str,
                          logger: Optional[logging.Logger] = None):
    """
    Create UMAP plots colored by perturbation status.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    output_dir : str
        Output directory for plots
    logger : logging.Logger, optional
        Logger instance
    """
    if logger:
        logger.info("Generating UMAP plots...")
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # UMAP colored by perturbation group
    sc.pl.umap(adata, color='perturbation_group', ax=axes[0], 
              title='Perturbation Groups', show=False, frameon=False)
    
    # UMAP colored by target perturbation
    adata.obs['is_target'] = (adata.obs['perturbation_group'] == Config.TARGET_PERTURBATION)
    sc.pl.umap(adata, color='is_target', ax=axes[1], 
              title=f'{Config.TARGET_PERTURBATION} vs Others', 
              show=False, frameon=False)
    
    # UMAP colored by control
    adata.obs['is_control'] = (adata.obs['perturbation_group'] == Config.CONTROL_PERTURBATION)
    sc.pl.umap(adata, color='is_control', ax=axes[2], 
              title=f'{Config.CONTROL_PERTURBATION} vs Others', 
              show=False, frameon=False)
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'umap_perturbation.pdf', dpi=Config.DPI, bbox_inches='tight')
    plt.close()
    
    if logger:
        logger.info("UMAP plots saved")

# ============================================================================
# DIFFERENTIAL EXPRESSION
# ============================================================================

def perform_differential_expression(adata: ad.AnnData, output_dir: str,
                                   logger: Optional[logging.Logger] = None) -> pd.DataFrame:
    """
    Perform differential expression analysis.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    output_dir : str
        Output directory for results
    logger : logging.Logger, optional
        Logger instance
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with differential expression results
    """
    if logger:
        logger.info("Performing differential expression analysis...")
        logger.info(f"Method: {Config.DE_METHOD}")
    
    # Perform DE analysis
    sc.tl.rank_genes_groups(
        adata, 
        'perturbation_group', 
        method=Config.DE_METHOD,
        groups=[Config.TARGET_PERTURBATION],
        reference=Config.CONTROL_PERTURBATION,
        n_genes=Config.N_TOP_GENES
    )
    
    # Extract results
    de_results = pd.DataFrame({
        'gene': adata.uns['rank_genes_groups']['names'][Config.TARGET_PERTURBATION],
        'logfoldchange': adata.uns['rank_genes_groups']['logfoldchanges'][Config.TARGET_PERTURBATION],
        'pval': adata.uns['rank_genes_groups']['pvals'][Config.TARGET_PERTURBATION],
        'pval_adj': adata.uns['rank_genes_groups']['pvals_adj'][Config.TARGET_PERTURBATION],
        'scores': adata.uns['rank_genes_groups']['scores'][Config.TARGET_PERTURBATION]
    })
    
    # Save results
    output_file = Path(output_dir) / f'differential_expression_{Config.TARGET_PERTURBATION}.csv'
    de_results.to_csv(output_file, index=False)
    
    if logger:
        logger.info(f"Differential expression results saved to {output_file}")
        n_sig = (de_results['pval_adj'] < 0.05).sum()
        logger.info(f"Significant genes (adj. p < 0.05): {n_sig}")
        logger.info(f"Top 10 upregulated genes:\n{de_results.head(10)[['gene', 'logfoldchange', 'pval_adj']]}")
    
    return de_results

def plot_de_results(adata: ad.AnnData, de_results: pd.DataFrame, output_dir: str,
                   logger: Optional[logging.Logger] = None):
    """
    Create visualizations for differential expression results.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    de_results : pd.DataFrame
        Differential expression results
    output_dir : str
        Output directory for plots
    logger : logging.Logger, optional
        Logger instance
    """
    if logger:
        logger.info("Generating DE visualization plots...")
    
    # Volcano plot
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Volcano plot
    de_results['-log10_pval'] = -np.log10(de_results['pval_adj'] + 1e-300)
    de_results['significant'] = de_results['pval_adj'] < 0.05
    
    scatter = axes[0].scatter(
        de_results['logfoldchange'], 
        de_results['-log10_pval'],
        c=de_results['significant'],
        cmap='RdYlBu_r',
        alpha=0.6,
        s=30
    )
    axes[0].axhline(-np.log10(0.05), color='red', linestyle='--', label='p_adj = 0.05')
    axes[0].axvline(0, color='black', linestyle='-', alpha=0.3)
    axes[0].set_xlabel('Log2 Fold Change')
    axes[0].set_ylabel('-Log10 Adjusted P-value')
    axes[0].set_title('Volcano Plot: Differential Expression')
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    
    # Top genes bar plot
    top_genes = de_results.head(20)
    colors = ['red' if x > 0 else 'blue' for x in top_genes['logfoldchange']]
    axes[1].barh(range(len(top_genes)), top_genes['logfoldchange'], color=colors)
    axes[1].set_yticks(range(len(top_genes)))
    axes[1].set_yticklabels(top_genes['gene'])
    axes[1].set_xlabel('Log2 Fold Change')
    axes[1].set_title('Top 20 Differentially Expressed Genes')
    axes[1].axvline(0, color='black', linestyle='-', alpha=0.3)
    axes[1].grid(alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'de_volcano_plot.pdf', dpi=Config.DPI, bbox_inches='tight')
    plt.close()
    
    # Violin plots for top genes
    top_n = 8
    top_genes_list = de_results.head(top_n)['gene'].tolist()
    
    sc.pl.rank_genes_groups_violin(
        adata, 
        groups=[Config.TARGET_PERTURBATION],
        n_genes=top_n,
        save=f'_top_{top_n}_genes.pdf'
    )
    
    if logger:
        logger.info("DE visualization plots saved")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

def generate_summary_report(adata: ad.AnnData, de_results: pd.DataFrame, 
                           output_dir: str, logger: Optional[logging.Logger] = None):
    """
    Generate a comprehensive summary report.
    
    Parameters:
    -----------
    adata : ad.AnnData
        AnnData object
    de_results : pd.DataFrame
        Differential expression results
    output_dir : str
        Output directory for report
    logger : logging.Logger, optional
        Logger instance
    """
    if logger:
        logger.info("Generating summary report...")
    
    report = {
        'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'dataset_info': {
            'n_cells': int(adata.n_obs),
            'n_genes': int(adata.n_vars),
            'n_highly_variable_genes': int(adata.var['highly_variable'].sum()),
            'n_clusters': len(adata.obs['leiden'].unique())
        },
        'perturbation_info': {
            'target': Config.TARGET_PERTURBATION,
            'control': Config.CONTROL_PERTURBATION,
            'n_target_cells': int((adata.obs['perturbation_group'] == Config.TARGET_PERTURBATION).sum()),
            'n_control_cells': int((adata.obs['perturbation_group'] == Config.CONTROL_PERTURBATION).sum())
        },
        'qc_metrics': {
            'mean_genes_per_cell': float(adata.obs['n_genes_by_counts'].mean()),
            'mean_counts_per_cell': float(adata.obs['total_counts'].mean()),
            'mean_mito_pct': float(adata.obs['pct_counts_mt'].mean()) if 'pct_counts_mt' in adata.obs.columns else None
        },
        'differential_expression': {
            'n_significant_genes': int((de_results['pval_adj'] < 0.05).sum()),
            'n_upregulated': int(((de_results['pval_adj'] < 0.05) & (de_results['logfoldchange'] > 0)).sum()),
            'n_downregulated': int(((de_results['pval_adj'] < 0.05) & (de_results['logfoldchange'] < 0)).sum()),
            'top_upregulated': de_results[de_results['logfoldchange'] > 0].head(10)['gene'].tolist(),
            'top_downregulated': de_results[de_results['logfoldchange'] < 0].head(10)['gene'].tolist()
        }
    }
    
    # Save JSON report
    report_file = Path(output_dir) / 'analysis_summary.json'
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    # Save text report
    report_text = f"""
Single-Cell RNA Perturbation Analysis Summary Report
{'=' * 60}

Analysis Date: {report['analysis_date']}

Dataset Information:
  - Total cells: {report['dataset_info']['n_cells']:,}
  - Total genes: {report['dataset_info']['n_genes']:,}
  - Highly variable genes: {report['dataset_info']['n_highly_variable_genes']:,}
  - Clusters identified: {report['dataset_info']['n_clusters']}

Perturbation Information:
  - Target perturbation: {report['perturbation_info']['target']}
  - Control perturbation: {report['perturbation_info']['control']}
  - Target cells: {report['perturbation_info']['n_target_cells']:,}
  - Control cells: {report['perturbation_info']['n_control_cells']:,}

Quality Control Metrics:
  - Mean genes per cell: {report['qc_metrics']['mean_genes_per_cell']:.1f}
  - Mean counts per cell: {report['qc_metrics']['mean_counts_per_cell']:.0f}
  - Mean mitochondrial %: {report['qc_metrics']['mean_mito_pct']:.2f}%

Differential Expression Results:
  - Significant genes (adj. p < 0.05): {report['differential_expression']['n_significant_genes']}
  - Upregulated genes: {report['differential_expression']['n_upregulated']}
  - Downregulated genes: {report['differential_expression']['n_downregulated']}
  
  Top 10 Upregulated Genes:
{chr(10).join(f'    {i+1}. {gene}' for i, gene in enumerate(report['differential_expression']['top_upregulated']))}
  
  Top 10 Downregulated Genes:
{chr(10).join(f'    {i+1}. {gene}' for i, gene in enumerate(report['differential_expression']['top_downregulated']))}

{'=' * 60}
"""
    
    report_text_file = Path(output_dir) / 'analysis_summary.txt'
    with open(report_text_file, 'w') as f:
        f.write(report_text)
    
    if logger:
        logger.info(f"Summary report saved to {report_file}")
        logger.info(f"Text report saved to {report_text_file}")
        print("\n" + report_text)

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    """Main analysis pipeline."""
    # Set random seeds for reproducibility
    random.seed(Config.RANDOM_SEED)
    np.random.seed(Config.RANDOM_SEED)
    
    # Create output directory
    output_dir = Path(Config.OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    logger = setup_logging(Config.OUTPUT_DIR)
    
    try:
        # 1. Load data
        adata = load_data(Config.DATA_FILE, Config.N_CELLS_TO_SUBSAMPLE, logger)
        
        # 2. Quality control
        adata = calculate_qc_metrics(adata, logger)
        plot_qc_metrics(adata, Config.OUTPUT_DIR, logger)
        adata = filter_data(adata, logger)
        
        # 3. Preprocessing
        adata = preprocess_data(adata, logger)
        
        # 4. Dimensionality reduction
        adata = find_hvg_and_reduce_dimensions(adata, logger)
        
        # 5. Perturbation analysis
        adata = prepare_perturbation_labels(adata, logger)
        plot_perturbation_umap(adata, Config.OUTPUT_DIR, logger)
        
        # 6. Differential expression
        de_results = perform_differential_expression(adata, Config.OUTPUT_DIR, logger)
        plot_de_results(adata, de_results, Config.OUTPUT_DIR, logger)
        
        # 7. Save processed data
        output_file = output_dir / 'processed_data.h5ad'
        adata.write(output_file)
        logger.info(f"Processed data saved to {output_file}")
        
        # 8. Generate summary report
        generate_summary_report(adata, de_results, Config.OUTPUT_DIR, logger)
        
        logger.info("=" * 80)
        logger.info("Analysis pipeline completed successfully!")
        logger.info(f"All results saved to: {output_dir.absolute()}")
        logger.info("=" * 80)
        
    except Exception as e:
        logger.error(f"Error in analysis pipeline: {e}", exc_info=True)
        raise

if __name__ == '__main__':
    main()

