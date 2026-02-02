"""
MK4 Disease Biomarker Package v1.0.0

Frequency-domain biomarkers for disease detection from RNA-seq data.

Version 1.0 — Proof of Concept
- Core framework: FFT + Shannon entropy chaos analysis
- Preliminary results: Cancer (78%), MS (79%), T2D (92%)
- Key discovery: Dataset balance > sample size
- Scope: RNA-level diseases (not protein-level)
- Priority: Reducing false detections in upcoming versions

See docs/roadmap.md for development plans.
"""

from .engine import MK4Analyzer, mk4_transform, calculate_chaos
from .analysis import analyze_disease_dataset, find_optimal_threshold
from .visualization import plot_chaos_distribution, plot_results_summary
from .utils import load_expression_data, preprocess_expression_data

__version__ = "1.0.0"
__author__ = "Architekt (Tomas Vavra)"
__organization__ = "Alexandria Dynamics"

__all__ = [
    'MK4Analyzer',
    'mk4_transform',
    'calculate_chaos',
    'analyze_disease_dataset',
    'find_optimal_threshold',
    'plot_chaos_distribution',
    'plot_results_summary',
    'load_expression_data',
    'preprocess_expression_data',
]
