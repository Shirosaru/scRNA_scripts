"""
Quick Visualization Script
==========================

Quickly visualize specific perturbations without running the full analysis.
Use this to quickly check how perturbations look on UMAP before running full analysis.
"""

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import random
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white')

# ============================================================================
# CONFIGURATION - EDIT THESE
# ============================================================================

DATA_FILE = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'
PERTURBATION_COLUMN = 'gene_target'

# Which perturbations to visualize (can be a list)
TARGET_PERTURBATIONS = ['TP53']  # Add more: ['TP53', 'BRCA1', 'RB1']
CONTROL_PERTURBATION = 'AAV-control'

# Subsampling (set to None to use all cells, or a number like 50000)
N_CELLS = 50000  # Smaller number = faster but less cells

# Random seed for reproducibility
RANDOM_SEED = 42

# ============================================================================
# MAIN SCRIPT
# ============================================================================

def quick_visualize():
    """Quick visualization of perturbations."""
    
    print("=" * 80)
    print("QUICK PERTURBATION VISUALIZATION")
    print("=" * 80)
    
    # Set random seed
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)
    
    # Load data
    print(f"\n📂 Loading data from: {DATA_FILE}")
    adata_backed = ad.read_h5ad(DATA_FILE, backed='r')
    print(f"   Total cells: {adata_backed.n_obs:,}")
    
    # Subsample if requested
    if N_CELLS and N_CELLS < adata_backed.n_obs:
        print(f"\n📊 Subsampling to {N_CELLS:,} cells...")
        all_indices = list(range(adata_backed.n_obs))
        subsample_indices = random.sample(all_indices, N_CELLS)
        adata = adata_backed[subsample_indices, :].to_memory()
    else:
        print(f"\n📊 Using all {adata_backed.n_obs:,} cells...")
        adata = adata_backed.to_memory()
    
    print(f"   Loaded {adata.n_obs:,} cells and {adata.n_vars:,} genes")
    
    # Check perturbations exist
    if PERTURBATION_COLUMN not in adata.obs.columns:
        print(f"\n❌ ERROR: Column '{PERTURBATION_COLUMN}' not found!")
        print(f"   Available columns: {list(adata.obs.columns)}")
        return
    
    # Check if perturbations exist
    available_perturbations = adata.obs[PERTURBATION_COLUMN].cat.categories
    print(f"\n🔍 Available perturbations: {len(available_perturbations)}")
    
    # Filter to only requested perturbations
    all_perturbations = list(TARGET_PERTURBATIONS) + [CONTROL_PERTURBATION]
    valid_perturbations = [p for p in all_perturbations if p in available_perturbations]
    
    if not valid_perturbations:
        print(f"\n❌ ERROR: None of the requested perturbations found!")
        print(f"   Requested: {all_perturbations}")
        print(f"   Available (first 20): {list(available_perturbations[:20])}")
        return
    
    # Filter data to only include requested perturbations
    print(f"\n🔬 Filtering to perturbations: {valid_perturbations}")
    adata = adata[adata.obs[PERTURBATION_COLUMN].isin(valid_perturbations), :]
    print(f"   Cells after filtering: {adata.n_obs:,}")
    
    # Show cell counts
    print(f"\n📊 Cell counts per perturbation:")
    counts = adata.obs[PERTURBATION_COLUMN].value_counts()
    for pert, count in counts.items():
        print(f"   • {pert}: {count:,} cells")
    
    # Quick preprocessing
    print(f"\n⚙️  Preprocessing...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print(f"   After filtering: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Dimensionality reduction
    print(f"\n📉 Computing UMAP...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.pca(adata, use_highly_variable=True, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)
    print("   ✓ UMAP computed")
    
    # Create visualizations
    print(f"\n🎨 Creating visualizations...")
    
    # Figure 1: All perturbations colored
    fig, axes = plt.subplots(1, len(TARGET_PERTURBATIONS) + 1, 
                            figsize=(5 * (len(TARGET_PERTURBATIONS) + 1), 5))
    
    if len(TARGET_PERTURBATIONS) == 1:
        axes = [axes]  # Make it iterable
    
    # Plot all perturbations together
    sc.pl.umap(adata, color=PERTURBATION_COLUMN, ax=axes[0], 
              title='All Perturbations', show=False, frameon=False, legend_loc='right margin')
    
    # Plot each target perturbation separately
    for idx, target in enumerate(TARGET_PERTURBATIONS, 1):
        if target in adata.obs[PERTURBATION_COLUMN].cat.categories:
            adata.obs[f'is_{target}'] = (adata.obs[PERTURBATION_COLUMN] == target)
            sc.pl.umap(adata, color=f'is_{target}', ax=axes[idx],
                      title=f'{target} vs Others', show=False, frameon=False)
    
    plt.tight_layout()
    output_file = 'quick_umap_visualization.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"   ✓ Saved: {output_file}")
    plt.close()
    
    # Figure 2: Side-by-side comparison
    if CONTROL_PERTURBATION in valid_perturbations:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Target vs Control
        for target in TARGET_PERTURBATIONS:
            if target in valid_perturbations:
                # Create comparison group
                adata.obs['comparison'] = 'other'
                adata.obs.loc[adata.obs[PERTURBATION_COLUMN] == target, 'comparison'] = target
                adata.obs.loc[adata.obs[PERTURBATION_COLUMN] == CONTROL_PERTURBATION, 'comparison'] = CONTROL_PERTURBATION
                
                # Plot
                sc.pl.umap(adata, color='comparison', ax=axes[0],
                          title=f'{target} vs {CONTROL_PERTURBATION}',
                          groups=[target, CONTROL_PERTURBATION], show=False, frameon=False)
                
                # Plot control
                adata.obs['is_control'] = (adata.obs[PERTURBATION_COLUMN] == CONTROL_PERTURBATION)
                sc.pl.umap(adata, color='is_control', ax=axes[1],
                          title=f'{CONTROL_PERTURBATION} (Control)', show=False, frameon=False)
                
                break
        
        plt.tight_layout()
        output_file2 = 'quick_comparison_visualization.png'
        plt.savefig(output_file2, dpi=300, bbox_inches='tight')
        print(f"   ✓ Saved: {output_file2}")
        plt.close()
    
    print(f"\n✅ Visualization complete!")
    print(f"\n📁 Output files:")
    print(f"   • {output_file}")
    if CONTROL_PERTURBATION in valid_perturbations:
        print(f"   • {output_file2}")
    
    print(f"\n💡 To analyze more perturbations, edit TARGET_PERTURBATIONS in this script")
    print(f"   or run the full analysis pipeline: python scrna_perturb_professional.py")


if __name__ == '__main__':
    quick_visualize()

