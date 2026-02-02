# MK4 Biomarker - Dataset Documentation

This directory contains information about datasets used in MK4 analysis.

## Current Status

**Important:** Current version uses **synthetic data** for demonstration.
All analysis scripts are ready for real data integration.

## Dataset Overview

| Dataset | GEO ID | Disease | Tissue | n (D/C) | Balance | Status |
|---------|--------|---------|--------|---------|---------|--------|
| Cancer | TBD | Multiple cancers | Blood | 500/42 | 8.4% | Synthetic |
| MS | GSE21942 | Multiple Sclerosis | Blood PBMCs | 14/15 | 48.3% | Synthetic |
| T2D-1 | GSE38642 | Type 2 Diabetes | Islets | 10/53 | 15.9% | Synthetic |
| T2D-2 | GSE25724 | Type 2 Diabetes | Islets | 6/7 | 46.2% | Synthetic |
| T2D-3 | GSE76894 | Type 2 Diabetes | Islets | 19/84 | 18.4% | Synthetic |
| PD | GSE49036 | Parkinson's | Brain SN | 15/8 | 65.2% | Synthetic |

## How to Add Real Data

### Step 1: Download from GEO

```bash
# Example for MS dataset (GSE21942)
# Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21942
# Download: Series Matrix File(s)
```

### Step 2: Convert to CSV

Place files in this directory:
- `cancer_expression.csv` - Expression matrix (samples x genes)
- `cancer_labels.csv` - Sample labels (sample_id, label)
- `ms_expression.csv`
- `ms_labels.csv`
- etc.

### Step 3: Update Analysis Scripts

In each `experiments/*_analysis.py`, replace synthetic data loading:

```python
# Current (synthetic):
expression_matrix, labels = generate_synthetic_data()

# Replace with:
from mk4.utils import load_expression_data, load_labels
expression_matrix = load_expression_data('data/cancer_expression.csv')
labels = load_labels('data/cancer_labels.csv')
```

## Data Format

### Expression Matrix CSV

```
gene_id,sample_1,sample_2,sample_3,...
GENE001,5.23,4.87,5.91,...
GENE002,8.12,7.94,8.45,...
```

### Labels CSV

```
sample_id,label,group
sample_1,0,control
sample_2,1,disease
sample_3,0,control
```

Where: `label` 0 = control, 1 = disease.

## Quality Control

Before analysis, data should meet:

- **Minimum genes:** >1,000 (recommend >10,000)
- **Minimum samples:** >5 per group
- **No missing values:** NaN handled or removed
- **Balanced data:** 30-70% disease samples (optimal: 40-60%)

## Citation

When using these datasets, please cite the original study authors.
See [CITATIONS.md](../CITATIONS.md) for complete references.

---

**Questions?** Open an issue: https://github.com/Alexandria-dynamics/mk4-biomarker/issues
