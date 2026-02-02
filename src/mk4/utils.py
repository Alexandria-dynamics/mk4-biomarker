"""
MK4 Utilities - Data loading, validation, preprocessing, and I/O.

Provides helper functions for handling expression data,
saving results, and printing summaries.
"""

import numpy as np
import json
import os


def load_expression_data(filepath):
    """
    Load expression matrix from CSV file.

    Expects rows = samples, columns = genes.
    First row may be a header (gene names).

    Parameters
    ----------
    filepath : str
        Path to CSV file.

    Returns
    -------
    np.ndarray
        Expression matrix (n_samples, n_genes).
    """
    import pandas as pd
    df = pd.read_csv(filepath, index_col=0)
    return df.values.astype(float)


def load_labels(filepath):
    """
    Load sample labels from CSV file.

    Parameters
    ----------
    filepath : str
        Path to CSV or text file with one label per line.

    Returns
    -------
    np.ndarray
        Binary labels (0=control, 1=disease).
    """
    import pandas as pd
    df = pd.read_csv(filepath)
    # Try common column names
    for col in ['label', 'Label', 'class', 'Class', 'group', 'Group']:
        if col in df.columns:
            labels = df[col].values
            break
    else:
        labels = df.iloc[:, 0].values

    # Convert to binary if needed
    unique = np.unique(labels)
    if len(unique) == 2:
        if set(unique) != {0, 1}:
            mapping = {unique[0]: 0, unique[1]: 1}
            labels = np.array([mapping[v] for v in labels])
    return labels.astype(int)


def validate_expression_data(expression_matrix):
    """
    Validate expression matrix for common issues.

    Parameters
    ----------
    expression_matrix : np.ndarray
        Expression matrix (n_samples, n_genes).

    Returns
    -------
    bool
        True if valid.

    Raises
    ------
    ValueError
        If data has critical issues.
    """
    if expression_matrix.ndim != 2:
        raise ValueError(f"Expected 2D matrix, got {expression_matrix.ndim}D")

    n_samples, n_genes = expression_matrix.shape
    print(f"  Shape: {n_samples} samples x {n_genes} genes")

    if n_samples < 2:
        raise ValueError("Need at least 2 samples")

    # Check for NaN/Inf
    n_nan = np.sum(np.isnan(expression_matrix))
    n_inf = np.sum(np.isinf(expression_matrix))
    if n_nan > 0:
        print(f"  Warning: {n_nan} NaN values found")
    if n_inf > 0:
        print(f"  Warning: {n_inf} Inf values found")
    if n_nan == 0 and n_inf == 0:
        print(f"  No NaN/Inf values")

    # Check variance
    variances = np.var(expression_matrix, axis=0)
    zero_var = np.sum(variances == 0)
    if zero_var > 0:
        print(f"  Warning: {zero_var} genes with zero variance")

    print(f"  Value range: [{expression_matrix.min():.2f}, {expression_matrix.max():.2f}]")
    print(f"  Validation passed")
    return True


def preprocess_expression_data(expression_matrix, remove_low_variance=True,
                               variance_percentile=10.0):
    """
    Preprocess expression matrix.

    Parameters
    ----------
    expression_matrix : np.ndarray
        Expression matrix (n_samples, n_genes).
    remove_low_variance : bool
        Whether to remove low-variance genes.
    variance_percentile : float
        Percentile threshold for variance filtering.

    Returns
    -------
    np.ndarray
        Preprocessed expression matrix.
    """
    data = expression_matrix.copy()

    # Replace NaN with gene means
    nan_mask = np.isnan(data)
    if np.any(nan_mask):
        gene_means = np.nanmean(data, axis=0)
        for j in range(data.shape[1]):
            data[nan_mask[:, j], j] = gene_means[j]
        print(f"  Replaced {np.sum(nan_mask)} NaN values with gene means")

    # Remove low-variance genes
    if remove_low_variance:
        variances = np.var(data, axis=0)
        threshold = np.percentile(variances, variance_percentile)
        keep = variances > threshold
        n_removed = np.sum(~keep)
        data = data[:, keep]
        print(f"  Removed {n_removed} low-variance genes (threshold: {threshold:.4f})")
        print(f"  Remaining: {data.shape[1]} genes")

    return data


def create_balanced_dataset(expression_matrix, labels, target_ratio=1.0, random_state=42):
    """
    Create a balanced dataset by downsampling the majority class.

    Parameters
    ----------
    expression_matrix : np.ndarray
        Expression matrix.
    labels : np.ndarray
        Binary labels.
    target_ratio : float
        Target ratio of minority/majority (1.0 = equal).
    random_state : int
        Random seed.

    Returns
    -------
    tuple
        (balanced_expression, balanced_labels)
    """
    rng = np.random.RandomState(random_state)
    classes, counts = np.unique(labels, return_counts=True)
    minority_class = classes[np.argmin(counts)]
    majority_class = classes[np.argmax(counts)]
    n_minority = int(np.min(counts))
    n_majority_target = int(n_minority / target_ratio)

    minority_idx = np.where(labels == minority_class)[0]
    majority_idx = np.where(labels == majority_class)[0]
    majority_sample = rng.choice(majority_idx, size=n_majority_target, replace=False)

    idx = np.concatenate([minority_idx, majority_sample])
    rng.shuffle(idx)

    print(f"  Balanced: {n_minority} minority + {n_majority_target} majority = {len(idx)} total")

    return expression_matrix[idx], labels[idx]


def save_results(results, output_dir="results/tables", prefix="analysis"):
    """
    Save analysis results to JSON.

    Parameters
    ----------
    results : dict
        Results dictionary.
    output_dir : str
        Output directory.
    prefix : str
        Filename prefix.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Convert numpy types for JSON serialization
    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [convert(v) for v in obj]
        return obj

    filepath = os.path.join(output_dir, f'{prefix}_results.json')
    with open(filepath, 'w') as f:
        json.dump(convert(results), f, indent=2)
    print(f"  Results saved: {filepath}")


def print_summary(expression_matrix, labels):
    """
    Print dataset summary.

    Parameters
    ----------
    expression_matrix : np.ndarray
        Expression matrix.
    labels : np.ndarray
        Sample labels.
    """
    n_samples, n_genes = expression_matrix.shape
    classes, counts = np.unique(labels, return_counts=True)

    print(f"\n--- Dataset Summary ---")
    print(f"  Samples: {n_samples}")
    print(f"  Genes: {n_genes}")
    for cls, cnt in zip(classes, counts):
        name = "Control" if cls == 0 else "Disease"
        print(f"  {name} (class {cls}): {cnt} samples")
    print(f"  Class ratio: {counts.max()/counts.min():.1f}:1")
