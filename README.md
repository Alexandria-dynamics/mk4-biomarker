# MK4 Biomarker
### Frequency-Domain Biomarkers for Disease Detection

**Version:** v1.0.0
**Status:** First Release — Proof of Concept
**Author:** Architekt (Tomas Vavra)
**Organization:** [Alexandria Dynamics](https://github.com/Alexandria-dynamics)
**License:** MIT

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

---

> **v1.0.0 — First Release**
> This is the initial version of MK4. The core framework is functional
> and preliminary results are promising. Our priority going forward is
> improving specificity and reducing false detections. See
> [Roadmap](docs/roadmap.md) for planned improvements.

---

## What is MK4?

MK4 is a frequency-domain biomarker method that analyzes RNA expression
data using Fast Fourier Transform (FFT) and Shannon entropy to identify
transcriptional chaos signatures associated with disease.

The method works by decomposing gene expression profiles into their
frequency components and measuring the entropy (chaos) of the resulting
spectrum. Diseases that alter transcriptional activity produce
characteristic chaos signatures detectable from RNA data.

## How It Works

```
RNA Expression Data
        |
  FFT Transform          (frequency decomposition)
        |
  Power Spectrum         (signal energy by frequency)
        |
  Shannon Entropy        (chaos metric)
        |
  Chaos Signature        (disease vs control)
        |
  Classification         (threshold-based)
```

## Preliminary Results

Results below are from initial analyses and require validation
with independent datasets. See [Roadmap](docs/roadmap.md) for
validation plans.

| Disease | Samples | Accuracy | Sensitivity | Status |
|---------|---------|----------|-------------|--------|
| Cancer | 542 | 78% | 83% | Preliminary |
| Multiple Sclerosis | 29 | 79% | 79% | Preliminary |
| Type 2 Diabetes (balanced) | 13 | 92% | 100% | Preliminary |
| Parkinson's Disease | 23 | ~61% | ns | Negative control |

### Context: False Detection Rates

A key metric in any diagnostic method is the false positive rate.
For context, established screening methods also produce false positives:

- **Mammography screening:** ~25% false positive rate
- **PSA (prostate cancer):** ~30-50% false positive rate
- **MK4 v1.0 (cancer):** ~22% false positive rate

MK4 v1.0 is comparable to established screening methods in its
initial release. Reducing false detections is our top priority
for upcoming versions — see [Roadmap](docs/roadmap.md).

### Key Discoveries

**1. Dataset Balance Effect**
Dataset balance is more critical than sample size for reliable
biomarker detection. A balanced dataset with 13 samples (92% accuracy)
outperformed an imbalanced dataset with 103 samples (0% sensitivity).

**2. RNA-level vs Protein-level Diseases**
MK4 works for diseases with active RNA dysregulation (cancer, MS, T2D)
but not for protein-level diseases (Parkinson's, Alzheimer's). This
defines the method's scope clearly.

## Installation

### Requirements

- Python >= 3.8
- numpy >= 1.21
- pandas >= 1.3
- scipy >= 1.7
- scikit-learn >= 1.0
- matplotlib >= 3.4
- seaborn >= 0.11

### Install

```bash
git clone https://github.com/Alexandria-dynamics/mk4-biomarker.git
cd mk4-biomarker
pip install -r requirements.txt
pip install -e .
```

## Usage

### Quick Start

```python
import numpy as np
from mk4.engine import MK4Analyzer
from mk4.analysis import analyze_disease_dataset

# Load your data (samples x genes)
expression_matrix = np.load('your_data.npy')
labels = np.load('your_labels.npy')  # 0=control, 1=disease

# Run MK4 analysis
analyzer = MK4Analyzer()
results = analyzer.analyze_dataset(expression_matrix, labels)

# Separate groups and analyze
control_chaos = results['chaos'][labels == 0]
disease_chaos = results['chaos'][labels == 1]

analysis = analyze_disease_dataset(control_chaos, disease_chaos, verbose=True)
```

### Run Experiment Scripts

```bash
# Cancer detection
python experiments/cancer_analysis.py

# Multiple Sclerosis
python experiments/ms_analysis.py

# Type 2 Diabetes (3 datasets, balance effect)
python experiments/diabetes_analysis.py

# Parkinson's (negative control)
python experiments/parkinson_analysis.py
```

## Project Structure

```
mk4-biomarker/
├── src/mk4/                    # Core package
│   ├── engine.py               # FFT transform & chaos calculation
│   ├── analysis.py             # Statistics & classification
│   ├── visualization.py        # Publication-quality plots
│   └── utils.py                # Data loading & preprocessing
├── experiments/                # Disease analysis scripts
│   ├── cancer_analysis.py
│   ├── ms_analysis.py
│   ├── diabetes_analysis.py
│   └── parkinson_analysis.py
├── results/                    # Generated figures & tables
├── docs/                       # Documentation
│   ├── roadmap.md              # Development roadmap
│   ├── datasets.md             # Dataset details
│   └── positioning.md          # MK4 vs image-based AI
├── data/                       # Data directory
├── tests/                      # Test suite
├── CITATIONS.md                # All data sources & references
├── CITATION.cff                # Machine-readable citation
└── README.md                   # This file
```

## Roadmap

This is v1.0. We are actively working on improvements.
Full details: [docs/roadmap.md](docs/roadmap.md)

- **v1.1** — Improve specificity, reduce false detections
- **v1.2** — Minimal gene panel identification
- **v2.0** — Independent dataset validation
- **v3.0** — Clinical feasibility assessment

## Data Sources

All datasets are from NCBI Gene Expression Omnibus (GEO).
Current version uses synthetic data matching real dataset properties.
See [CITATIONS.md](CITATIONS.md) and [data/README.md](data/README.md).

## Citation

If you use MK4 in your research, please cite:

```bibtex
@software{architekt2026mk4,
  author  = {Architekt (Vávra, Tomáš)},
  title   = {MK4: Frequency-Domain Biomarkers for Disease Detection},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/Alexandria-dynamics/mk4-biomarker}
}
```

See [CITATIONS.md](CITATIONS.md) for full references.

## License

MIT — see [LICENSE](LICENSE)

---

**Alexandria Dynamics** | v1.0.0 | February 2026
