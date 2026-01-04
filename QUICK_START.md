# Quick Start Guide

## How to Look at Your Perturbation Data

### 🚀 Three Simple Steps

### 1️⃣ **Explore What's in Your Dataset**

First, see what perturbations you have:

```bash
python explore_data.py
```

**What this does:**
- Shows you ALL available perturbations
- Tells you how many cells each perturbation has
- Identifies potential control groups
- Creates a visualization of the distribution
- Suggests which perturbations are good for analysis

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

### 2️⃣ **Quick Visual Check (Optional)**

Want to quickly see how a specific perturbation looks? Use the quick visualization:

```bash
# 1. Edit quick_visualize.py - change these lines:
TARGET_PERTURBATIONS = ['TP53']  # Your gene of interest
CONTROL_PERTURBATION = 'AAV-control'  # Your control

# 2. Run it:
python quick_visualize.py
```

**What this does:**
- Loads a subset of cells (fast!)
- Computes UMAP
- Shows your perturbation colored on the UMAP
- Saves quick visualization images

**Output:** `quick_umap_visualization.png` - See your perturbation on UMAP

### 3️⃣ **Run Full Professional Analysis**

Once you know what perturbations you want to analyze:

```bash
# 1. Edit scrna_perturb_professional.py - change these in the Config class:
TARGET_PERTURBATION = 'TP53'      # Your target gene
CONTROL_PERTURBATION = 'AAV-control'  # Your control
N_CELLS_TO_SUBSAMPLE = 100000     # Or None for all cells

# 2. Run it:
python scrna_perturb_professional.py
```

**What this does:**
- Full quality control analysis
- Comprehensive preprocessing
- Dimensionality reduction (PCA, UMAP, clustering)
- Differential expression analysis
- Professional visualizations
- Summary reports

**Output:** Everything saved in `perturbation_analysis_results/` folder

---

## 📋 Common Questions

### Q: How do I know which column has my perturbations?

**A:** Run `explore_data.py` - it shows all columns. Look for columns like:
- `gene_target`
- `pert_name`
- `perturbation`
- `guide`

Then update `PERTURBATION_COLUMN` in the scripts.

### Q: How do I find controls?

**A:** Run `explore_data.py` - it automatically searches for common control names like:
- "control", "ctrl", "AAV-control", "non-targeting", etc.

### Q: My dataset is huge, will it be slow?

**A:** All scripts support subsampling:
- `explore_data.py` - Uses backed mode (memory efficient, fast)
- `quick_visualize.py` - Set `N_CELLS = 50000` for speed
- `scrna_perturb_professional.py` - Set `N_CELLS_TO_SUBSAMPLE = 100000`

### Q: Can I analyze multiple perturbations?

**A:** Yes! See `example_usage.py` for batch analysis examples.

### Q: What if my perturbation column has a different name?

**A:** Edit the `PERTURBATION_COLUMN` variable in each script:
```python
PERTURBATION_COLUMN = 'your_column_name'
```

---

## 🎯 Recommended Workflow

1. **First time?** → Run `explore_data.py` to understand your data
2. **Want quick preview?** → Run `quick_visualize.py` with your favorite perturbation
3. **Ready for full analysis?** → Run `scrna_perturb_professional.py`
4. **Multiple perturbations?** → Use `example_usage.py` as a template

---

## 📁 Output Files Explained

### From `explore_data.py`:
- `perturbation_distribution.png` - Overview of all perturbations

### From `quick_visualize.py`:
- `quick_umap_visualization.png` - Quick UMAP view
- `quick_comparison_visualization.png` - Target vs control

### From `scrna_perturb_professional.py`:
- `qc_metrics.pdf` - Quality control plots
- `umap_perturbation.pdf` - UMAP visualizations
- `differential_expression_<gene>.csv` - DE results table
- `de_volcano_plot.pdf` - Volcano plot
- `analysis_summary.json` - Full summary (JSON)
- `analysis_summary.txt` - Full summary (human-readable)
- `processed_data.h5ad` - Processed dataset

---

## 💡 Pro Tips

1. **Start small**: Use subsampling first to test your analysis
2. **Check cell counts**: Perturbations with <100 cells may have limited power
3. **Save time**: Use `quick_visualize.py` to preview before full analysis
4. **Batch processing**: Use `example_usage.py` to analyze multiple genes at once

---

## 🆘 Troubleshooting

**"Column not found" error?**
→ Run `explore_data.py` to see available columns, then update `PERTURBATION_COLUMN`

**"Perturbation not found" error?**
→ Run `explore_data.py` to see available perturbations, check spelling

**Out of memory?**
→ Reduce `N_CELLS_TO_SUBSAMPLE` or set it to a smaller number

**Slow analysis?**
→ Use smaller subsample size, or run `quick_visualize.py` first to test

