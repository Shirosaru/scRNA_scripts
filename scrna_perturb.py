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

# --- SCRIPT START ---

try:
    # 1. Load the dataset
    print(f"Loading data from {file_path}...")
    adata = sc.read_h5ad(file_path)

    # Print some basic information about the loaded AnnData object
    print("\n--- Initial Dataset Info ---")
    print(adata)
    print(f"Available cell metadata columns: {list(adata.obs.columns)}")
    print(f"Available gene metadata columns: {list(adata.var.columns)}")

    # 2. Basic Preprocessing and Quality Control
    print("\nStarting preprocessing and quality control...")
    
    # Simple filtering: keep cells with at least 200 genes expressed
    sc.pp.filter_cells(adata, min_genes=200)

    # Simple filtering: keep genes expressed in at least 3 cells
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Normalize data to a total count per cell and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    print("\n--- Processed Dataset Info ---")
    print(adata)

    # 3. Dimensionality Reduction
    print("\nPerforming dimensionality reduction...")

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Run PCA on highly variable genes
    sc.pp.pca(adata, use_highly_variable=True)
    
    # Run UMAP to create a 2D visualization
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)

    # 4. Visualization of Perturbation Effect
    print(f"\nVisualizing the effect of the '{target_pert_name}' perturbation...")

    # Create a new column in the metadata to highlight our target and control cells.
    adata.obs['is_target'] = adata.obs['pert_name'] == target_pert_name
    adata.obs['is_control'] = adata.obs['pert_name'] == control_pert_name
    
    # Plot the UMAP colored by whether a cell received the perturbation
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # Plot all cells, highlighting the target perturbation
    sc.pl.umap(adata, color='is_target', title=f"Cells with {target_pert_name} Perturbation", ax=axes[0], show=False)
    
    # Plot all cells, highlighting the control perturbation
    sc.pl.umap(adata, color='is_control', title=f"Cells with {control_pert_name} (Control) Perturbation", ax=axes[1], show=False)

    plt.tight_layout()
    plt.show()
    
    print("\nScript completed successfully.")
    print("The UMAP plots show the distribution of your target and control cells.")
    print("If the perturbation has a strong effect, you may see the target cells clustered in a distinct region of the UMAP.")

except FileNotFoundError:
    print(f"Error: The file at '{file_path}' was not found.")
    print("Please make sure you have replaced 'path/to/your/downloaded/file.h5ad' with the correct path to your data file.")
except Exception as e:
    print(f"An error occurred: {e}")