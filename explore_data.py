"""
Interactive Data Exploration Script
===================================

This script helps you explore your single-cell perturbation data:
- See what perturbations are available
- Check cell counts per perturbation
- Explore metadata columns
- Quick visualizations
- Browse the dataset structure

Run this BEFORE running the full analysis to understand your data!
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

DATA_FILE = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'
PERTURBATION_COLUMN = 'gene_target'  # Change this if your column has a different name

# ============================================================================
# EXPLORATION FUNCTIONS
# ============================================================================

def load_data_quick(file_path: str):
    """Load data in backed mode for quick exploration without loading everything into memory."""
    print("=" * 80)
    print("Loading data (backed mode - memory efficient)...")
    print("=" * 80)
    
    try:
        adata = ad.read_h5ad(file_path, backed='r')
        return adata
    except FileNotFoundError:
        print(f"\n❌ ERROR: File not found at: {file_path}")
        print("\nPlease update the DATA_FILE path in this script!")
        return None
    except Exception as e:
        print(f"\n❌ ERROR loading data: {e}")
        return None


def show_dataset_overview(adata):
    """Show basic information about the dataset."""
    print("\n" + "=" * 80)
    print("📊 DATASET OVERVIEW")
    print("=" * 80)
    print(f"\nTotal cells: {adata.n_obs:,}")
    print(f"Total genes: {adata.n_vars:,}")
    print(f"\nData shape: {adata.shape}")
    print(f"\nData type: {type(adata.X)}")
    
    # Show if data is sparse
    try:
        from scipy.sparse import issparse
        if issparse(adata.X):
            print("Data format: Sparse matrix (memory efficient)")
        else:
            print("Data format: Dense matrix")
    except:
        pass


def show_metadata_columns(adata):
    """Show all available metadata columns."""
    print("\n" + "=" * 80)
    print("📋 METADATA COLUMNS (Cell-level annotations)")
    print("=" * 80)
    
    obs_cols = list(adata.obs.columns)
    print(f"\nFound {len(obs_cols)} metadata columns:\n")
    
    for i, col in enumerate(obs_cols, 1):
        dtype = adata.obs[col].dtype
        n_unique = adata.obs[col].nunique()
        print(f"  {i:2d}. {col:30s} | Type: {str(dtype):15s} | Unique values: {n_unique:,}")
    
    print("\n" + "-" * 80)
    print("💡 TIP: Look for columns that might contain perturbation information")
    print("   Common names: 'gene_target', 'pert_name', 'perturbation', 'guide', etc.")
    print("-" * 80)
    
    return obs_cols


def explore_perturbation_column(adata, pert_col):
    """Explore a specific column that contains perturbation information."""
    print("\n" + "=" * 80)
    print(f"🔬 EXPLORING PERTURBATION COLUMN: '{pert_col}'")
    print("=" * 80)
    
    if pert_col not in adata.obs.columns:
        print(f"\n❌ Column '{pert_col}' not found in the dataset!")
        print(f"\nAvailable columns: {list(adata.obs.columns)}")
        return None
    
    # Get unique perturbations
    perturbations = adata.obs[pert_col].value_counts()
    
    print(f"\nFound {len(perturbations)} unique perturbations:\n")
    print(f"{'Rank':<6} {'Perturbation':<30} {'Cell Count':<15} {'Percentage':<15}")
    print("-" * 80)
    
    total_cells = perturbations.sum()
    for rank, (pert, count) in enumerate(perturbations.items(), 1):
        pct = (count / total_cells) * 100
        print(f"{rank:<6} {str(pert):<30} {count:>12,} {pct:>13.2f}%")
    
    print("-" * 80)
    print(f"{'TOTAL':<6} {'':<30} {total_cells:>12,} {'100.00%':>13}")
    
    # Show top and bottom perturbations
    print("\n" + "-" * 80)
    print("📈 TOP 10 PERTURBATIONS (most cells):")
    print("-" * 80)
    for pert, count in perturbations.head(10).items():
        print(f"  • {pert}: {count:,} cells")
    
    print("\n" + "-" * 80)
    print("📉 BOTTOM 10 PERTURBATIONS (fewest cells):")
    print("-" * 80)
    for pert, count in perturbations.tail(10).items():
        print(f"  • {pert}: {count:,} cells")
    
    return perturbations


def show_sample_metadata(adata, pert_col, n_samples=5):
    """Show sample rows of metadata."""
    print("\n" + "=" * 80)
    print(f"📄 SAMPLE METADATA (first {n_samples} cells)")
    print("=" * 80)
    
    sample_df = adata.obs.head(n_samples)
    print("\n", sample_df.to_string())
    
    print("\n" + "-" * 80)
    print("💡 This shows what the metadata looks like for individual cells")
    print("-" * 80)


def visualize_perturbation_distribution(adata, pert_col, output_file='perturbation_distribution.png'):
    """Create a visualization of perturbation distribution."""
    print("\n" + "=" * 80)
    print("📊 CREATING VISUALIZATION")
    print("=" * 80)
    
    if pert_col not in adata.obs.columns:
        print(f"❌ Column '{pert_col}' not found!")
        return
    
    perturbations = adata.obs[pert_col].value_counts()
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Bar plot of top 20 perturbations
    top_20 = perturbations.head(20)
    axes[0, 0].barh(range(len(top_20)), top_20.values)
    axes[0, 0].set_yticks(range(len(top_20)))
    axes[0, 0].set_yticklabels(top_20.index, fontsize=8)
    axes[0, 0].set_xlabel('Number of Cells', fontsize=10)
    axes[0, 0].set_title('Top 20 Perturbations by Cell Count', fontsize=12, fontweight='bold')
    axes[0, 0].grid(axis='x', alpha=0.3)
    axes[0, 0].invert_yaxis()
    
    # 2. Pie chart of top 10
    top_10 = perturbations.head(10)
    other_count = perturbations.iloc[10:].sum()
    if other_count > 0:
        pie_data = list(top_10.values) + [other_count]
        pie_labels = list(top_10.index) + ['Others']
    else:
        pie_data = top_10.values
        pie_labels = top_10.index
    
    axes[0, 1].pie(pie_data, labels=pie_labels, autopct='%1.1f%%', startangle=90)
    axes[0, 1].set_title('Distribution of Top 10 Perturbations', fontsize=12, fontweight='bold')
    
    # 3. Histogram of cell counts per perturbation
    axes[1, 0].hist(perturbations.values, bins=50, edgecolor='black', alpha=0.7)
    axes[1, 0].set_xlabel('Number of Cells per Perturbation', fontsize=10)
    axes[1, 0].set_ylabel('Number of Perturbations', fontsize=10)
    axes[1, 0].set_title('Distribution of Cell Counts per Perturbation', fontsize=12, fontweight='bold')
    axes[1, 0].grid(alpha=0.3)
    axes[1, 0].axvline(perturbations.median(), color='red', linestyle='--', 
                       label=f'Median: {perturbations.median():.0f}')
    axes[1, 0].legend()
    
    # 4. Log scale histogram
    axes[1, 1].hist(perturbations.values, bins=50, edgecolor='black', alpha=0.7)
    axes[1, 1].set_xlabel('Number of Cells per Perturbation (log scale)', fontsize=10)
    axes[1, 1].set_ylabel('Number of Perturbations (log scale)', fontsize=10)
    axes[1, 1].set_title('Distribution (Log Scale)', fontsize=12, fontweight='bold')
    axes[1, 1].set_xscale('log')
    axes[1, 1].set_yscale('log')
    axes[1, 1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✅ Visualization saved to: {output_file}")
    plt.close()


def find_control_perturbations(adata, pert_col):
    """Try to identify control perturbations."""
    print("\n" + "=" * 80)
    print("🔍 SEARCHING FOR CONTROL PERTURBATIONS")
    print("=" * 80)
    
    perturbations = adata.obs[pert_col].value_counts()
    
    # Common control names
    control_keywords = ['control', 'ctrl', 'non-targeting', 'nt', 'scramble', 
                       'empty', 'mock', 'aav-control', 'neg', 'negative']
    
    potential_controls = []
    for pert in perturbations.index:
        pert_lower = str(pert).lower()
        for keyword in control_keywords:
            if keyword in pert_lower:
                potential_controls.append(pert)
                break
    
    if potential_controls:
        print("\n✅ Found potential control perturbations:")
        for control in potential_controls:
            count = perturbations[control]
            print(f"  • {control}: {count:,} cells")
    else:
        print("\n⚠️  No obvious control perturbations found.")
        print("   You may need to manually identify controls based on your experimental design.")
    
    return potential_controls


def suggest_analysis_parameters(adata, pert_col):
    """Suggest parameters for analysis based on the data."""
    print("\n" + "=" * 80)
    print("💡 SUGGESTIONS FOR ANALYSIS")
    print("=" * 80)
    
    perturbations = adata.obs[pert_col].value_counts()
    total_cells = adata.n_obs
    
    print(f"\n📊 Dataset Statistics:")
    print(f"  • Total cells: {total_cells:,}")
    print(f"  • Total perturbations: {len(perturbations)}")
    print(f"  • Median cells per perturbation: {perturbations.median():.0f}")
    print(f"  • Mean cells per perturbation: {perturbations.mean():.0f}")
    
    # Suggest subsampling
    if total_cells > 200000:
        suggested_subsample = min(100000, int(total_cells * 0.5))
        print(f"\n💾 Memory Management:")
        print(f"  • Dataset is large ({total_cells:,} cells)")
        print(f"  • Suggested subsample size: {suggested_subsample:,} cells")
        print(f"  • This will speed up analysis while maintaining statistical power")
    
    # Suggest perturbations to analyze
    print(f"\n🎯 Perturbations with sufficient cells for analysis (>100 cells):")
    sufficient = perturbations[perturbations >= 100]
    print(f"  • {len(sufficient)} perturbations have ≥100 cells")
    print(f"  • Top candidates: {', '.join(sufficient.head(10).index.astype(str).tolist())}")
    
    print(f"\n⚠️  Perturbations with few cells (<100 cells):")
    few = perturbations[perturbations < 100]
    print(f"  • {len(few)} perturbations have <100 cells")
    print(f"  • These may have limited statistical power")


def interactive_explore_column(adata, col_name):
    """Interactively explore a specific column."""
    if col_name not in adata.obs.columns:
        print(f"❌ Column '{col_name}' not found!")
        return
    
    print(f"\n" + "=" * 80)
    print(f"🔍 DETAILED EXPLORATION: '{col_name}'")
    print("=" * 80)
    
    col_data = adata.obs[col_name]
    
    print(f"\nData type: {col_data.dtype}")
    print(f"Unique values: {col_data.nunique()}")
    print(f"Missing values: {col_data.isna().sum()}")
    
    if col_data.dtype == 'object' or col_data.dtype.name == 'category':
        # Categorical data
        value_counts = col_data.value_counts()
        print(f"\nValue distribution:")
        print(value_counts.head(20).to_string())
    else:
        # Numerical data
        print(f"\nStatistical summary:")
        print(col_data.describe())


# ============================================================================
# MAIN EXPLORATION SCRIPT
# ============================================================================

def main():
    """Main exploration function."""
    print("\n" + "=" * 80)
    print("🔬 SINGLE-CELL PERTURBATION DATA EXPLORATION")
    print("=" * 80)
    print("\nThis script will help you explore your dataset before running analysis.")
    print("It loads data in memory-efficient mode, so it's fast even for large datasets.\n")
    
    # Load data
    adata = load_data_quick(DATA_FILE)
    if adata is None:
        return
    
    # Show overview
    show_dataset_overview(adata)
    
    # Show metadata columns
    metadata_cols = show_metadata_columns(adata)
    
    # Explore perturbation column
    print("\n" + "=" * 80)
    print(f"🔬 EXPLORING PERTURBATION COLUMN: '{PERTURBATION_COLUMN}'")
    print("=" * 80)
    
    if PERTURBATION_COLUMN in adata.obs.columns:
        perturbations = explore_perturbation_column(adata, PERTURBATION_COLUMN)
        
        # Find controls
        controls = find_control_perturbations(adata, PERTURBATION_COLUMN)
        
        # Create visualization
        print("\nGenerating visualization...")
        visualize_perturbation_distribution(adata, PERTURBATION_COLUMN)
        
        # Show sample metadata
        show_sample_metadata(adata, PERTURBATION_COLUMN)
        
        # Suggestions
        suggest_analysis_parameters(adata, PERTURBATION_COLUMN)
        
    else:
        print(f"\n⚠️  Column '{PERTURBATION_COLUMN}' not found!")
        print("\nAvailable columns:")
        for col in metadata_cols:
            print(f"  • {col}")
        print(f"\n💡 Update PERTURBATION_COLUMN in this script to match your data!")
    
    # Summary
    print("\n" + "=" * 80)
    print("✅ EXPLORATION COMPLETE")
    print("=" * 80)
    print("\nNext steps:")
    print("  1. Review the perturbation list above")
    print("  2. Identify your target perturbation(s) and control(s)")
    print("  3. Update the configuration in 'scrna_perturb_professional.py'")
    print("  4. Run the full analysis pipeline")
    print("\n" + "=" * 80)


if __name__ == '__main__':
    main()

