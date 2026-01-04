"""
RUN ME FIRST - Simple Data Explorer
====================================

This is the simplest script to run first. It will show you:
1. What perturbations are in your data
2. How many cells each has
3. What columns are available

Just run: python RUN_ME_FIRST.py
"""

import anndata as ad
import pandas as pd
from pathlib import Path

# ============================================================================
# CONFIGURATION - EDIT THIS
# ============================================================================

DATA_FILE = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'
PERTURBATION_COLUMN = 'gene_target'  # Change if your column has different name

# ============================================================================
# MAIN SCRIPT
# ============================================================================

print("=" * 80)
print("SIMPLE DATA EXPLORER - RUN ME FIRST")
print("=" * 80)
print()

try:
    # Load data (memory efficient)
    print(f"📂 Loading: {DATA_FILE}")
    print("   (This may take a moment...)")
    adata = ad.read_h5ad(DATA_FILE, backed='r')
    
    print(f"\n✅ Data loaded successfully!")
    print(f"   Total cells: {adata.n_obs:,}")
    print(f"   Total genes: {adata.n_vars:,}")
    
    # Show metadata columns
    print(f"\n📋 Available metadata columns ({len(adata.obs.columns)}):")
    for i, col in enumerate(adata.obs.columns, 1):
        print(f"   {i:2d}. {col}")
    
    # Check perturbation column
    print(f"\n🔬 Checking perturbation column: '{PERTURBATION_COLUMN}'")
    if PERTURBATION_COLUMN not in adata.obs.columns:
        print(f"   ❌ NOT FOUND!")
        print(f"   Please update PERTURBATION_COLUMN in this script.")
        print(f"   Available columns are listed above.")
    else:
        # Show perturbations
        perturbations = adata.obs[PERTURBATION_COLUMN].value_counts()
        print(f"   ✅ Found {len(perturbations)} unique perturbations\n")
        
        print("   Top 20 Perturbations:")
        print("   " + "-" * 70)
        print(f"   {'Rank':<6} {'Perturbation':<40} {'Cells':<15}")
        print("   " + "-" * 70)
        
        for rank, (pert, count) in enumerate(perturbations.head(20).items(), 1):
            print(f"   {rank:<6} {str(pert):<40} {count:>12,}")
        
        print("   " + "-" * 70)
        print(f"   {'TOTAL':<6} {'':<40} {perturbations.sum():>12,}")
        
        # Find controls
        print(f"\n   🔍 Looking for controls...")
        control_keywords = ['control', 'ctrl', 'non-targeting', 'nt', 'scramble', 
                           'empty', 'mock', 'aav-control', 'neg', 'negative']
        found_controls = []
        for pert in perturbations.index:
            pert_lower = str(pert).lower()
            for keyword in control_keywords:
                if keyword in pert_lower:
                    found_controls.append(pert)
                    break
        
        if found_controls:
            print(f"   ✅ Potential controls found:")
            for ctrl in found_controls:
                print(f"      • {ctrl}: {perturbations[ctrl]:,} cells")
        else:
            print(f"   ⚠️  No obvious controls found")
    
    print("\n" + "=" * 80)
    print("NEXT STEPS:")
    print("=" * 80)
    print("1. Note down your target perturbation name (e.g., 'TP53')")
    print("2. Note down your control perturbation name (e.g., 'AAV-control')")
    print("3. Edit scrna_perturb_professional.py:")
    print("   - Set TARGET_PERTURBATION = 'your_target'")
    print("   - Set CONTROL_PERTURBATION = 'your_control'")
    print("4. Run: python scrna_perturb_professional.py")
    print("=" * 80)
    
except FileNotFoundError:
    print(f"\n❌ ERROR: File not found!")
    print(f"   Path: {DATA_FILE}")
    print(f"\n   Please update DATA_FILE in this script with the correct path.")
    
except Exception as e:
    print(f"\n❌ ERROR: {e}")
    print(f"\n   Please check:")
    print(f"   1. File path is correct")
    print(f"   2. You have read permissions")
    print(f"   3. File is a valid .h5ad file")

