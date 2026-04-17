# mk4.eeg — Spectral Biomarker Pipeline for EEG

Companion module to the original MK4 (RNA expression) for EEG time-series analysis.

## Files

- `batch.py` — per-subject feature extraction (Welch PSD, 10 spectral features) and per-dataset orientation-vector optimization via grid search maximizing Youden's J.
- `stats.py` — statistical validation: 5-fold stratified cross-validation, permutation testing, bootstrap 95% confidence intervals. Strict no-data-leakage protocol.
- `figures.py` — publication figures (cosine heatmap, performance barchart, feature distributions, pipeline schematic, PCA embedding).

## Quickstart

```bash
pip install -r requirements.txt

# Download an OpenNeuro dataset with datalad (example: Alzheimer/FTD)
datalad clone https://github.com/OpenNeuroDatasets/ds004504.git
cd ds004504
datalad get sub-*/eeg/*.set sub-*/eeg/*.fdt
cd ..

# Feature extraction + initial weight training
python -m mk4.eeg.batch ./ds004504

# Statistical validation
python -m mk4.eeg.stats ds004504 ./mk4e_results_ds004504.json --cv --perm 1000 --boot 1000

# Figures
export MK4_EEG_DATA=.
export MK4_EEG_FIGURES=./figures
python -m mk4.eeg.figures
```

## Datasets used in the companion paper

- ds002778 (Parkinson)
- ds003478 (MDD)
- ds003523 (mild TBI, task-based — used as domain-shift probe)
- ds004504 (Alzheimer + FTD)

## Reproducibility

All random seeds are fixed to 42. Python 3.14, MNE ≥1.10, scikit-learn ≥1.0.

See validation outputs in `results/eeg/mk4e_stats_*.json`.
