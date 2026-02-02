"""
MK4 Analysis - Statistical analysis and classification.

Provides functions for comparing groups, finding optimal thresholds,
computing classification metrics, and running complete disease analyses.
"""

import numpy as np
from scipy import stats
from sklearn.metrics import roc_curve, auc, confusion_matrix


def compare_groups(control_scores, disease_scores):
    """
    Statistical comparison between control and disease groups.

    Parameters
    ----------
    control_scores : np.ndarray
        Chaos scores for control samples.
    disease_scores : np.ndarray
        Chaos scores for disease samples.

    Returns
    -------
    dict
        't_stat': t-statistic,
        'p_value': p-value (two-sided),
        'effect_size': Cohen's d,
        'control_mean': mean of control,
        'disease_mean': mean of disease,
        'control_std': std of control,
        'disease_std': std of disease.
    """
    control = np.asarray(control_scores)
    disease = np.asarray(disease_scores)

    t_stat, p_value = stats.ttest_ind(control, disease, equal_var=False)

    # Cohen's d
    pooled_std = np.sqrt(
        ((len(control) - 1) * np.var(control, ddof=1) +
         (len(disease) - 1) * np.var(disease, ddof=1)) /
        (len(control) + len(disease) - 2)
    )
    effect_size = (np.mean(disease) - np.mean(control)) / pooled_std if pooled_std > 0 else 0.0

    # Mann-Whitney U (non-parametric)
    u_stat, u_pvalue = stats.mannwhitneyu(control, disease, alternative='two-sided')

    return {
        't_stat': t_stat,
        'p_value': p_value,
        'effect_size': effect_size,
        'control_mean': np.mean(control),
        'disease_mean': np.mean(disease),
        'control_std': np.std(control, ddof=1),
        'disease_std': np.std(disease, ddof=1),
        'u_stat': u_stat,
        'u_pvalue': u_pvalue,
    }


def find_optimal_threshold(control_scores, disease_scores):
    """
    Find the optimal classification threshold using Youden's J.

    Parameters
    ----------
    control_scores : np.ndarray
        Chaos scores for control samples.
    disease_scores : np.ndarray
        Chaos scores for disease samples.

    Returns
    -------
    dict
        'threshold': optimal threshold,
        'fpr': false positive rates,
        'tpr': true positive rates,
        'thresholds': all thresholds from ROC,
        'auc': area under ROC curve.
    """
    scores = np.concatenate([control_scores, disease_scores])
    labels = np.concatenate([
        np.zeros(len(control_scores)),
        np.ones(len(disease_scores))
    ])

    fpr, tpr, thresholds = roc_curve(labels, scores)
    roc_auc = auc(fpr, tpr)

    # Youden's J statistic
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    optimal_threshold = thresholds[optimal_idx]

    return {
        'threshold': optimal_threshold,
        'fpr': fpr,
        'tpr': tpr,
        'thresholds': thresholds,
        'auc': roc_auc,
    }


def calculate_metrics(control_scores, disease_scores, threshold):
    """
    Calculate classification metrics at a given threshold.

    Samples with score >= threshold are classified as disease.

    Parameters
    ----------
    control_scores : np.ndarray
        Chaos scores for control samples.
    disease_scores : np.ndarray
        Chaos scores for disease samples.
    threshold : float
        Classification threshold.

    Returns
    -------
    dict
        'accuracy', 'sensitivity', 'specificity', 'ppv', 'npv',
        'tp', 'fp', 'tn', 'fn'.
    """
    control_pred = (np.asarray(control_scores) >= threshold).astype(int)
    disease_pred = (np.asarray(disease_scores) >= threshold).astype(int)

    tp = np.sum(disease_pred == 1)
    fn = np.sum(disease_pred == 0)
    fp = np.sum(control_pred == 1)
    tn = np.sum(control_pred == 0)

    total = tp + fn + fp + tn
    accuracy = (tp + tn) / total if total > 0 else 0.0
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0

    return {
        'accuracy': accuracy,
        'sensitivity': sensitivity,
        'specificity': specificity,
        'ppv': ppv,
        'npv': npv,
        'tp': int(tp),
        'fp': int(fp),
        'tn': int(tn),
        'fn': int(fn),
    }


def analyze_disease_dataset(control_scores, disease_scores, verbose=True):
    """
    Complete analysis of a disease dataset.

    Runs statistical comparison, finds optimal threshold,
    and computes classification metrics.

    Parameters
    ----------
    control_scores : np.ndarray
        Chaos scores for control samples.
    disease_scores : np.ndarray
        Chaos scores for disease samples.
    verbose : bool
        Whether to print results.

    Returns
    -------
    dict
        'statistics': group comparison results,
        'roc': ROC analysis results,
        'classification': classification metrics.
    """
    # Statistical comparison
    stat_results = compare_groups(control_scores, disease_scores)

    # ROC analysis
    roc_results = find_optimal_threshold(control_scores, disease_scores)

    # Classification at optimal threshold
    class_results = calculate_metrics(
        control_scores, disease_scores, roc_results['threshold']
    )

    if verbose:
        print(f"\n--- Group Comparison ---")
        print(f"  Control: {stat_results['control_mean']:.4f} +/- {stat_results['control_std']:.4f} "
              f"(n={len(control_scores)})")
        print(f"  Disease: {stat_results['disease_mean']:.4f} +/- {stat_results['disease_std']:.4f} "
              f"(n={len(disease_scores)})")
        print(f"  t-stat: {stat_results['t_stat']:.3f}")
        print(f"  p-value: {stat_results['p_value']:.6f}")
        print(f"  Cohen's d: {stat_results['effect_size']:.3f}")

        print(f"\n--- Classification ---")
        print(f"  Threshold: {roc_results['threshold']:.4f}")
        print(f"  AUC: {roc_results['auc']:.3f}")
        print(f"  Accuracy: {class_results['accuracy']*100:.1f}%")
        print(f"  Sensitivity: {class_results['sensitivity']*100:.1f}%")
        print(f"  Specificity: {class_results['specificity']*100:.1f}%")
        print(f"  PPV: {class_results['ppv']*100:.1f}%")
        print(f"  NPV: {class_results['npv']*100:.1f}%")

    return {
        'statistics': stat_results,
        'roc': roc_results,
        'classification': class_results,
    }
