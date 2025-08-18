import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
# Replace this with the actual path to one of your downloaded .h5ad files.
# Example: 'processed_data_lenti.h5ad' or 'processed_data_crispr.h5ad'
file_path = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'

# Name of the perturbation (gene) you want to visualize.
# You can change this to any gene in the dataset's 'pert_name' column.
# A common control group is 'AAV-control' or similar.
target_pert_name = 'TP53'
control_pert_name = 'AAV-control'

# Name of the perturbation (gene) you want to visualize.
target_pert_name = 'TP53'
control_pert_name = 'AAV-control'

# Define the number of cells to randomly subsample.
# 100,000 cells is a good starting point for a 64 GB machine.
n_cells_to_subsample = 100000

# Name of the output file where the subsampled data will be saved.
output_file_name = 'subsampled_analysis_results.h5ad'

# --- SCRIPT START ---

try:
    print(f"Loading data from {file_path}...")
    
    # Load the full dataset into memory.
    # This might take a moment but should not exceed your RAM.
    adata = sc.read_h5ad(file_path)

    # 1. Subsample the data to a manageable size.
    print(f"\nSubsampling to {n_cells_to_subsample} cells...")
    sc.pp.subsample(adata, n_obs=n_cells_to_subsample)
    
    print("\n--- Subsampled Dataset Info ---")
    print(adata)

    # 2. Preprocessing and Quality Control on the subsampled data.
    print("\nStarting preprocessing and quality control...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 3. Dimensionality Reduction
    print("\nPerforming dimensionality reduction...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.pca(adata, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)

    # 4. Save the processed data to a new file.
    print(f"\nSaving processed data to '{output_file_name}'...")
    adata.write(output_file_name, compression='gzip')
    print("File saved successfully.")

    # 5. Visualization of Perturbation Effect
    print(f"\nGenerating UMAP plots for '{target_pert_name}' and '{control_pert_name}'...")

    # Create a new column to highlight the cells of interest
    adata.obs['is_target_or_control'] = 'other'
    adata.obs.loc[adata.obs['pert_name'] == target_pert_name, 'is_target_or_control'] = target_pert_name
    adata.obs.loc[adata.obs['pert_name'] == control_pert_name, 'is_target_or_control'] = control_pert_name
    
    # Plot the UMAP colored by the perturbation status
    sc.pl.umap(adata, color='is_target_or_control', title=f"UMAP of Subsampled Data",
               groups=[target_pert_name, control_pert_name], size=30)
    plt.show()

    print("\nScript completed successfully.")

except FileNotFoundError:
    print(f"Error: The file at '{file_path}' was not found.")
    print("Please make sure you have replaced the placeholder path with the correct path to your data file.")
except Exception as e:
    print(f"An error occurred: {e}")