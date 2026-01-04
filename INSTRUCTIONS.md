# Step-by-Step Instructions for scRNA Perturbation Analysis

## 🎯 Quick Start (3 Steps)

### Step 1: Explore Your Data
**First, see what perturbations you have in your dataset**

```bash
python explore_data.py
```

**What to do:**
1. Open `explore_data.py` in a text editor
2. Check line 18 - make sure `DATA_FILE` path is correct:
   ```python
   DATA_FILE = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'
   ```
3. Check line 19 - make sure `PERTURBATION_COLUMN` matches your data:
   ```python
   PERTURBATION_COLUMN = 'gene_target'  # Change if your column has different name
   ```
4. Run the script:
   ```bash
   python explore_data.py
   ```

**What you'll get:**
- List of ALL perturbations with cell counts
- Visualization saved as `perturbation_distribution.png`
- Suggestions for which perturbations to analyze

**Example output:**
```
Found 150 unique perturbations:

Rank   Perturbation              Cell Count      Percentage
------------------------------------------------------------
1      TP53                      15,234         12.5%
2      AAV-control               14,891         12.2%
3      BRCA1                      8,234          6.7%
...
```

---

### Step 2: Choose Your Perturbations
**Based on Step 1, decide which perturbation(s) you want to analyze**

**Write down:**
- Your target perturbation (e.g., `TP53`)
- Your control perturbation (e.g., `AAV-control`)

---

### Step 3: Run the Analysis
**Edit and run the professional analysis script**

#### Option A: Quick Preview (Fast, ~5-10 minutes)
**Good for: Quick check before full analysis**

1. Open `quick_visualize.py`
2. Edit these lines (around line 26-31):
   ```python
   DATA_FILE = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'
   PERTURBATION_COLUMN = 'gene_target'
   
   TARGET_PERTURBATIONS = ['TP53']  # Change to your target gene(s)
   CONTROL_PERTURBATION = 'AAV-control'  # Change to your control
   
   N_CELLS = 50000  # Number of cells (smaller = faster)
   ```
3. Run:
   ```bash
   python quick_visualize.py
   ```
4. Check output: `quick_umap_visualization.png`

#### Option B: Full Professional Analysis (Complete, ~30-60 minutes)
**Good for: Publication-ready results**

1. Open `scrna_perturb_professional.py`
2. Find the `Config` class (around line 40)
3. Edit these settings:
   ```python
   class Config:
       # Data paths
       DATA_FILE = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'
       OUTPUT_DIR = 'perturbation_analysis_results'
       
       # Perturbation settings
       PERTURBATION_COLUMN = 'gene_target'  # Your perturbation column name
       TARGET_PERTURBATION = 'TP53'  # CHANGE THIS to your target
       CONTROL_PERTURBATION = 'AAV-control'  # CHANGE THIS to your control
       
       # Subsampling (set to None to use all cells)
       N_CELLS_TO_SUBSAMPLE = 100000  # Or None for all cells
   ```
4. Run:
   ```bash
   python scrna_perturb_professional.py
   ```
5. Check output folder: `perturbation_analysis_results/`

---

## 📋 Detailed Code Instructions

### For Multiple Perturbations

If you want to analyze multiple genes, use `example_usage.py`:

1. Open `example_usage.py`
2. Edit the perturbations list (around line 70):
   ```python
   perturbations = ['TP53', 'BRCA1', 'RB1']  # Your list of genes
   ```
3. Run:
   ```bash
   python example_usage.py
   ```

### Customizing Analysis Parameters

In `scrna_perturb_professional.py`, you can customize:

```python
class Config:
    # Quality control - adjust if needed
    MIN_GENES_PER_CELL = 200      # Minimum genes per cell
    MIN_CELLS_PER_GENE = 3        # Minimum cells per gene
    MAX_MITO_PCT = 20.0           # Maximum mitochondrial %
    
    # Differential expression method
    DE_METHOD = 'wilcoxon'        # Options: 'wilcoxon', 't-test', 'logreg'
    N_TOP_GENES = 50              # Number of top genes to report
```

---

## 🔧 Troubleshooting

### Problem: "File not found"
**Solution:**
1. Check the file path in the script
2. Make sure the path is correct (use absolute path)
3. Example:
   ```python
   DATA_FILE = '/home2/figshare_data/HEK293T_filtered_dual_guide_cells.h5ad'
   ```

### Problem: "Column not found"
**Solution:**
1. Run `explore_data.py` first to see available columns
2. Update `PERTURBATION_COLUMN` in your script
3. Example:
   ```python
   PERTURBATION_COLUMN = 'pert_name'  # Use the correct column name
   ```

### Problem: "Perturbation not found"
**Solution:**
1. Run `explore_data.py` to see available perturbations
2. Check spelling (case-sensitive!)
3. Example:
   ```python
   TARGET_PERTURBATION = 'TP53'  # Not 'tp53' or 'Tp53'
   ```

### Problem: Out of memory
**Solution:**
1. Reduce the number of cells:
   ```python
   N_CELLS_TO_SUBSAMPLE = 50000  # Smaller number
   ```
2. Or set to None and let it use all (if you have enough RAM):
   ```python
   N_CELLS_TO_SUBSAMPLE = None
   ```

### Problem: Script is too slow
**Solution:**
1. Use smaller subsample:
   ```python
   N_CELLS_TO_SUBSAMPLE = 50000  # Instead of 100000
   ```
2. Or use `quick_visualize.py` first for a preview

---

## 📊 Understanding the Output

### From `explore_data.py`:
- **`perturbation_distribution.png`**: Shows all perturbations and their cell counts

### From `quick_visualize.py`:
- **`quick_umap_visualization.png`**: UMAP plot with your perturbations colored
- **`quick_comparison_visualization.png`**: Target vs control comparison

### From `scrna_perturb_professional.py`:
All files saved in `perturbation_analysis_results/` folder:

1. **`qc_metrics.pdf`**: Quality control plots
   - Shows distribution of counts, genes, mitochondrial %

2. **`umap_perturbation.pdf`**: UMAP visualizations
   - Shows how perturbations cluster in 2D space

3. **`differential_expression_<GENE>.csv`**: DE results table
   - Columns: gene, logfoldchange, pval, pval_adj, scores
   - Sorted by significance

4. **`de_volcano_plot.pdf`**: Volcano plot
   - Shows fold change vs significance
   - Red = significant genes

5. **`rank_genes_groups_violin_top_8_genes.pdf`**: Violin plots
   - Shows expression of top DE genes

6. **`analysis_summary.txt`**: Human-readable summary
   - Quick overview of results

7. **`analysis_summary.json`**: Machine-readable summary
   - For programmatic access

8. **`processed_data.h5ad`**: Processed dataset
   - Can be loaded later for further analysis

---

## 🎯 Recommended Workflow

### First Time User:
```
1. python explore_data.py          # See what you have
2. python quick_visualize.py      # Quick preview
3. python scrna_perturb_professional.py  # Full analysis
```

### Experienced User:
```
1. Edit scrna_perturb_professional.py (set your perturbations)
2. python scrna_perturb_professional.py
3. Check results in perturbation_analysis_results/
```

### Batch Analysis (Multiple Genes):
```
1. Edit example_usage.py
2. python example_usage.py
3. Each gene gets its own output folder
```

---

## 💻 Code Snippets

### Check if a perturbation exists:
```python
import anndata as ad

adata = ad.read_h5ad('your_file.h5ad', backed='r')
perturbations = adata.obs['gene_target'].value_counts()
print(perturbations.head(20))  # See top 20
```

### Load and check your data:
```python
import anndata as ad

adata = ad.read_h5ad('your_file.h5ad', backed='r')
print(f"Cells: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")
print(f"Columns: {list(adata.obs.columns)}")
```

### Quick UMAP for one perturbation:
```python
import scanpy as sc
import anndata as ad

# Load
adata = ad.read_h5ad('your_file.h5ad', backed='r')
adata = adata[:50000].to_memory()  # Subsample

# Process
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot
sc.pl.umap(adata, color='gene_target', save='quick_look.pdf')
```

---

## ✅ Checklist Before Running

- [ ] Data file path is correct
- [ ] Perturbation column name is correct
- [ ] Target perturbation name is correct (case-sensitive!)
- [ ] Control perturbation name is correct
- [ ] Subsample size is appropriate for your memory
- [ ] Output directory exists or will be created

---

## 🆘 Need Help?

1. **Check the log file**: `perturbation_analysis_results/analysis_log_*.log`
2. **Run explore_data.py first**: Always helps identify issues
3. **Check file paths**: Use absolute paths
4. **Verify column names**: Run explore_data.py to see available columns
5. **Check perturbation names**: Must match exactly (case-sensitive)

---

## 📝 Example: Complete Workflow

```bash
# Step 1: Explore
python explore_data.py

# Step 2: Edit quick_visualize.py with your perturbations
# Then run:
python quick_visualize.py

# Step 3: Edit scrna_perturb_professional.py with your perturbations
# Then run:
python scrna_perturb_professional.py

# Step 4: Check results
ls perturbation_analysis_results/
cat perturbation_analysis_results/analysis_summary.txt
```

---

**That's it! Follow these steps and you'll have professional analysis results.**

