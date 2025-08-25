import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os
import random
import pandas as pd

# --- CONFIGURATION ---
# Replace with the actual path to one of your downloaded .h5ad files.
file_path = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'

# You have correctly identified 'gene_target' as the perturbation column.
perturbation_column = 'gene_target'

# Define the perturbation you want to visualize.
target_pert_name = 'TP53'
control_pert_name = 'AAV-control'

# Define the number of cells to randomly subsample.
n_cells_to_subsample = 100000

# Name of the output file where the subsampled data will be saved.
output_file_name = 'subsampled_analysis_results.h5ad'

# --- SCRIPT START ---

try:
    print(f"Opening file in memory-safe 'backed' mode...")

    # Load the AnnData object in backed mode.
    adata_backed = ad.read_h5ad(file_path, backed='r')
    
    print(f"\nTotal number of cells in the dataset: {adata_backed.n_obs}")
    
    # NEW: Print all available metadata columns for inspection.
    print(f"\nAvailable metadata columns in your dataset: {list(adata_backed.obs.columns)}")
    print("\nPlease look at this list to confirm the correct column for your perturbations is 'gene_target'.")
    
    # Check if target and control names are present in the subsampled data's gene_target column
    if target_pert_name not in adata_backed.obs[perturbation_column].cat.categories:
        print(f"Warning: Target perturbation '{target_pert_name}' not found in '{perturbation_column}' column.")
    if control_pert_name not in adata_backed.obs[perturbation_column].cat.categories:
        print(f"Warning: Control perturbation '{control_pert_name}' not found in '{perturbation_column}' column.")

    # 1. Randomly sample the indices of cells to keep.
    print(f"\nRandomly selecting {n_cells_to_subsample} cells to analyze...")
    all_indices = list(range(adata_backed.n_obs))
    subsample_indices = random.sample(all_indices, n_cells_to_subsample)
    
    # 2. Slice the backed object and use .to_memory() to load the data into RAM.
    adata = adata_backed[subsample_indices, :].to_memory()
    
    print("\n--- Subsampled Dataset Info ---")
    print(adata)
    
    # 3. Preprocessing and Quality Control on the subsampled data.
    print("\nStarting preprocessing and quality control...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 4. Dimensionality Reduction
    print("\nPerforming dimensionality reduction...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)

    # 5. Visualization of Perturbation Effect (and saving the UMAP plot)
    print(f"\nGenerating UMAP plots for '{target_pert_name}' and '{control_pert_name}'...")

    # Create a new column to highlight the cells of interest.
    adata.obs['is_target_or_control'] = 'other'
    adata.obs.loc[adata.obs[perturbation_column] == target_pert_name, 'is_target_or_control'] = target_pert_name
    adata.obs.loc[adata.obs[perturbation_column] == control_pert_name, 'is_target_or_control'] = control_pert_name
    
    # Plot the UMAP colored by the perturbation status and save to file.
    sc.pl.umap(adata, color='is_target_or_or_control', title=f"UMAP of Subsampled Data",
               groups=[target_pert_name, control_pert_name], size=30,
               save="_umap_perturbation_effect.pdf")

    # 6. Differential Expression Analysis (and saving the results)
    print("\nStarting differential expression analysis...")
    sc.tl.rank_genes_groups(adata, 'is_target_or_control', method='t-test')
    
    # Save the differential gene expression results to a CSV file.
    results = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    results.to_csv('TP53_differential_genes.csv', index=False)

    # Save the violin plots of the top genes to a file.
    sc.pl.rank_genes_groups_violin(adata, groups=[target_pert_name, control_pert_name], n_genes=8,
                                   save="_DE_violin_plots.pdf")

    # 7. Final outputs saved to a file
    print(f"\nAll outputs have been successfully saved to your working directory:")
    print(f"- Processed AnnData object: '{output_file_name}'")
    print(f"- UMAP plot: 'umap_is_target_or_or_control_umap_perturbation_effect.pdf'")
    print(f"- Differential expression results: 'TP53_differential_genes.csv'")
    print(f"- Differential expression violin plots: 'rank_genes_groups_violin_DE_violin_plots.pdf'")
    
    print("\nScript completed successfully.")

except FileNotFoundError:
    print(f"Error: The file at '{file_path}' was not found.")
    print("Please make sure you have provided the correct and absolute path to your data file.")
except Exception as e:
    print(f"An error occurred: {e}")