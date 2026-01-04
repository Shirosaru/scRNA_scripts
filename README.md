# Professional Single-Cell RNA Perturbation Analysis Pipeline

A comprehensive, production-ready analysis pipeline for single-cell RNA sequencing perturbation data.

## Features

- **Comprehensive Quality Control**: Automated QC metrics calculation and filtering
- **Advanced Preprocessing**: Normalization, log transformation, and highly variable gene detection
- **Dimensionality Reduction**: PCA, UMAP, and Leiden clustering
- **Differential Expression Analysis**: Statistical testing with multiple methods (Wilcoxon, t-test, logistic regression)
- **Professional Visualizations**: Publication-quality figures for QC, UMAP, and differential expression
- **Automated Reporting**: JSON and text summary reports
- **Reproducibility**: Configurable random seeds and detailed logging
- **Memory Efficient**: Supports large datasets with optional subsampling

## Installation

1. Install required packages:

```bash
pip install -r requirements.txt
```

2. Ensure you have the data file in the correct location (update `Config.DATA_FILE` in the script if needed).

## Usage

### Step 1: Explore Your Data First! 🔍

**Before running the full analysis, explore your dataset to see what perturbations are available:**

```bash
python explore_data.py
```

This will show you:
- All available perturbations and their cell counts
- Metadata columns in your dataset
- Distribution of perturbations
- Suggestions for analysis parameters
- Visualization of perturbation distribution

**Output:** `perturbation_distribution.png` - A comprehensive visualization of your perturbations

### Step 2: Quick Visualization (Optional) 🎨

**Quickly visualize specific perturbations on UMAP:**

```bash
# Edit quick_visualize.py to set your perturbations, then:
python quick_visualize.py
```

Edit the configuration in `quick_visualize.py`:
```python
TARGET_PERTURBATIONS = ['TP53', 'BRCA1']  # Your genes of interest
CONTROL_PERTURBATION = 'AAV-control'      # Your control
N_CELLS = 50000  # Number of cells to use (smaller = faster)
```

**Output:** Quick UMAP visualizations showing your perturbations

### Step 3: Run Full Analysis 📊

**After exploring, run the comprehensive analysis:**

```bash
python scrna_perturb_professional.py
```

### Example Scripts

See `example_usage.py` for examples of:
- Running with default configuration
- Customizing for different perturbation targets
- Batch analysis of multiple perturbations

### Configuration

Edit the `Config` class in `scrna_perturb_professional.py` to customize:

- **Data paths**: Set `DATA_FILE` and `OUTPUT_DIR`
- **Perturbation settings**: Configure `TARGET_PERTURBATION` and `CONTROL_PERTURBATION`
- **QC parameters**: Adjust filtering thresholds (`MIN_GENES_PER_CELL`, `MIN_CELLS_PER_GENE`, `MAX_MITO_PCT`)
- **Analysis parameters**: Modify dimensionality reduction and DE analysis settings
- **Subsampling**: Set `N_CELLS_TO_SUBSAMPLE` to None to use all cells, or specify a number

### Example Configuration

```python
class Config:
    DATA_FILE = '/path/to/your/data.h5ad'
    OUTPUT_DIR = 'results'
    TARGET_PERTURBATION = 'TP53'
    CONTROL_PERTURBATION = 'AAV-control'
    N_CELLS_TO_SUBSAMPLE = 50000  # or None for all cells
    DE_METHOD = 'wilcoxon'  # 'wilcoxon', 't-test', or 'logreg'
```

## Output Files

The pipeline generates the following outputs in the specified output directory:

1. **QC Metrics Plot** (`qc_metrics.pdf`): Distribution plots for counts, genes, and mitochondrial percentage
2. **UMAP Plots** (`umap_perturbation.pdf`): UMAP visualizations colored by perturbation status
3. **Differential Expression Results** (`differential_expression_<target>.csv`): Full DE results table
4. **DE Visualizations** (`de_volcano_plot.pdf`, `rank_genes_groups_violin_top_8_genes.pdf`): Volcano plot and violin plots
5. **Processed Data** (`processed_data.h5ad`): Processed AnnData object with all computed metrics
6. **Summary Reports** (`analysis_summary.json`, `analysis_summary.txt`): Comprehensive analysis summaries
7. **Log File** (`analysis_log_<timestamp>.log`): Detailed execution log

## Analysis Pipeline

The pipeline performs the following steps:

1. **Data Loading**: Loads AnnData object with optional subsampling for memory efficiency
2. **Quality Control**: Calculates QC metrics (counts, genes, mitochondrial percentage) and filters low-quality cells/genes
3. **Preprocessing**: Normalizes to 10,000 reads per cell and log-transforms
4. **Dimensionality Reduction**: 
   - Identifies highly variable genes
   - Computes PCA
   - Builds neighborhood graph
   - Computes UMAP embedding
   - Performs Leiden clustering
5. **Perturbation Analysis**: Prepares perturbation labels and visualizes on UMAP
6. **Differential Expression**: Performs statistical testing between target and control groups
7. **Visualization**: Generates publication-quality figures
8. **Reporting**: Creates comprehensive summary reports

## Key Improvements Over Basic Scripts

- **Modular Design**: Functions are well-organized and reusable
- **Error Handling**: Comprehensive error checking and informative error messages
- **Logging**: Detailed logging to both file and console
- **Reproducibility**: Configurable random seeds
- **Documentation**: Extensive docstrings and comments
- **Professional Visualizations**: Publication-ready figures
- **Statistical Rigor**: Multiple DE methods and proper multiple testing correction
- **Memory Management**: Efficient handling of large datasets

## Requirements

- Python 3.8+
- See `requirements.txt` for package versions

## Citation

If you use this pipeline in your research, please cite the underlying tools:

- Scanpy: Wolf et al., Genome Biology, 2018
- AnnData: Virshup et al., Nature Methods, 2021

## License

See LICENSE file for details.

## Support

For issues or questions, please check the log files in the output directory for detailed error messages.

