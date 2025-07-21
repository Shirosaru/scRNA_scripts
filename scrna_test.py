import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import os

# --- Configuration ---
# Define a list of example cytokine genes and tumor target genes.
# You can modify these lists with the specific genes you are interested in.
# Note: Gene names should match those in your scRNA-seq dataset (e.g., official gene symbols).
CYTOKINE_GENES = ['IL6', 'TNF', 'IFNG', 'CCL2', 'CXCL8']
TUMOR_TARGET_GENES = ['PDCD1', 'CD274', 'CTLA4', 'VEGFA', 'EGFR'] # PD-1, PD-L1, CTLA-4, VEGF-A, EGFR

# --- Data Loading ---
# For demonstration purposes, we will use a small, built-in PBMC dataset from scanpy.
# For actual cancer data, you would typically load an AnnData object (.h5ad file)
# that you have downloaded from a public repository (e.g., GEO, ENA, Broad Institute's Single Cell Portal).
#
# Example of how you might load your own .h5ad cancer data:
# DATA_FILE = 'your_cancer_data.h5ad'
# if not os.path.exists(DATA_FILE):
#     # You would need to download your specific cancer data file here.
#     # For example, using requests for a direct download link:
#     # data_url = "https://example.com/path/to/your_cancer_data.h5ad"
#     # print(f"Downloading {DATA_FILE} from {data_url}...")
#     # response = requests.get(data_url, stream=True)
#     # response.raise_for_status() # Raise an HTTPError for bad responses (4xx or 5xx)
#     # with open(DATA_FILE, 'wb') as f:
#     #     for chunk in response.iter_content(chunk_size=8192):
#     #         f.write(chunk)
#     # print("Download complete.")
# else:
#     print(f"Using existing data file: {DATA_FILE}")

print("Loading example PBMC 3k dataset from scanpy.datasets...")
# This loads a pre-processed AnnData object.
# For real cancer data, replace this with sc.read_h5ad(DATA_FILE) after downloading.
adata = sc.datasets.pbmc3k()
print("Dataset loaded successfully.")
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")

# --- Basic Preprocessing (Essential for meaningful analysis) ---
# For raw count data, these steps are crucial. If your .h5ad is already normalized/log-transformed,
# you might skip some of these.

print("Performing basic preprocessing: Normalization and Log-transformation...")
# 1. Normalize total counts per cell
# This ensures that differences in sequencing depth between cells don't mislead analysis.
sc.pp.normalize_total(adata, target_sum=1e4)

# 2. Logarithmize the data
# This helps to stabilize variance and make gene expression distributions more symmetric.
sc.pp.log1p(adata)

print("Preprocessing complete.")

# --- Dimensionality Reduction and Clustering (for visualization context) ---
# These steps help to visualize cells in a lower-dimensional space and identify cell clusters.
print("Performing dimensionality reduction (PCA) and UMAP...")
sc.pp.pca(adata) # Principal Component Analysis
sc.pp.neighbors(adata) # Compute a neighborhood graph of observations
sc.tl.umap(adata) # Uniform Manifold Approximation and Projection for dimension reduction
sc.tl.leiden(adata) # Leiden clustering algorithm to find cell clusters
print("UMAP and clustering complete.")

# --- Gene Expression Observation ---

# Combine all genes of interest into one list
genes_to_observe = CYTOKINE_GENES + TUMOR_TARGET_GENES

# Filter for genes actually present in the dataset
# It's common that some genes from your list might not be detected in a specific dataset.
present_genes = [gene for gene in genes_to_observe if gene in adata.var_names]
missing_genes = [gene for gene in genes_to_observe if gene not in adata.var_names]

if missing_genes:
    print(f"\nWarning: The following genes were not found in the dataset and will be skipped: {missing_genes}")
if not present_genes:
    print("\nError: No specified genes were found in the dataset. Please check gene names or dataset.")
else:
    print(f"\nAnalyzing expression for present genes: {present_genes}")

    # --- Visualization ---
    # 1. UMAP plots showing gene expression overlay
    # This visualizes the expression level of each gene on the UMAP embedding.
    print("Generating UMAP plots for gene expression...")
    sc.pl.umap(adata, color=present_genes, cmap='viridis', size=20,
               title=[f"{gene} Expression" for gene in present_genes], # This title argument is fine for sc.pl.umap
               ncols=3, # Adjust number of columns for display
               show=False) # Set to False to combine plots or save them programmatically
    plt.suptitle('Gene Expression on UMAP Embedding', fontsize=16, y=1.02)
    plt.tight_layout(rect=[0, 0.03, 1, 0.98]) # Adjust layout to prevent title overlap
    plt.show()

    # 2. Violin plots showing gene expression across clusters
    # This shows the distribution of gene expression within each identified cell cluster.
    print("Generating Violin plots for gene expression across clusters...")
    # Removed the 'title' and 'ncols' arguments from sc.pl.violin as they cause an AttributeError.
    # When 'keys' is a list, sc.pl.violin automatically creates a separate subplot for each gene.
    # If you need to arrange these plots in a grid, you would typically iterate and plot them
    # on manually created subplots, or save them individually.
    sc.pl.violin(adata, keys=present_genes, groupby='leiden', rotation=90,
                 show=False)
    plt.suptitle('Gene Expression Distribution by Leiden Cluster', fontsize=16, y=1.02)
    plt.tight_layout(rect=[0, 0.03, 1, 0.98])
    plt.show()

    print("\nAnalysis complete. Visualizations displayed.")

# --- Further Exploration ---
print("\n--- Next Steps for Your Analysis ---")
print("1. Replace the example PBMC dataset with your actual cancer scRNA-seq data (.h5ad).")
print("   You can often find processed .h5ad files in supplementary data of publications,")
print("   or on dedicated single-cell data portals like:")
print("   - Broad Institute's Single Cell Portal (portals.broadinstitute.org/single_cell)")
print("   - GEO (Gene Expression Omnibus) - look for processed matrices or AnnData objects.")
print("   - European Nucleotide Archive (ENA) or NCBI SRA (for raw data, requires more extensive preprocessing).")
print("2. Adjust `CYTOKINE_GENES` and `TUMOR_TARGET_GENES` lists to your specific interests.")
print("3. For more in-depth analysis (e.g., differential expression, cell type annotation),")
print("   you would expand the preprocessing and analysis steps using `scanpy`'s extensive functionalities.")
print("   - `sc.tl.rank_genes_groups()` for marker genes.")
print("   - `sc.tl.paga()` for trajectory inference.")
print("   - Integration of multiple datasets (batch correction).")
