"""
MK4 Visualization - Publication-quality plotting functions.

Generates figures for chaos distributions, group comparisons,
ROC curves, and comprehensive result summaries.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os


def plot_chaos_distribution(control_scores, disease_scores, disease_name="Disease",
                            ax=None, save_path=None):
    """
    Plot overlapping histograms of chaos scores.

    Parameters
    ----------
    control_scores : np.ndarray
        Chaos scores for control group.
    disease_scores : np.ndarray
        Chaos scores for disease group.
    disease_name : str
        Name of disease for labels.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.
    save_path : str, optional
        Path to save figure.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
        created_fig = True

    bins = np.linspace(
        min(np.min(control_scores), np.min(disease_scores)),
        max(np.max(control_scores), np.max(disease_scores)),
        30
    )

    ax.hist(control_scores, bins=bins, alpha=0.6, label='Control',
            color='#2196F3', edgecolor='white', density=True)
    ax.hist(disease_scores, bins=bins, alpha=0.6, label=disease_name,
            color='#F44336', edgecolor='white', density=True)

    ax.axvline(np.mean(control_scores), color='#1565C0', linestyle='--',
               linewidth=1.5, label=f'Control mean ({np.mean(control_scores):.3f})')
    ax.axvline(np.mean(disease_scores), color='#C62828', linestyle='--',
               linewidth=1.5, label=f'{disease_name} mean ({np.mean(disease_scores):.3f})')

    ax.set_xlabel('MK4 Chaos Score', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(f'MK4 Chaos Distribution: Control vs {disease_name}', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    if save_path and created_fig:
        fig.tight_layout()
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {save_path}")
        plt.close(fig)
        return None

    if created_fig:
        return fig
    return None


def plot_comparison_boxplot(control_scores, disease_scores, disease_name="Disease",
                           ax=None, save_path=None):
    """
    Plot box-and-whisker comparison of groups.

    Parameters
    ----------
    control_scores : np.ndarray
        Chaos scores for control group.
    disease_scores : np.ndarray
        Chaos scores for disease group.
    disease_name : str
        Name of disease for labels.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.
    save_path : str, optional
        Path to save figure.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
        created_fig = True

    bp = ax.boxplot(
        [control_scores, disease_scores],
        labels=['Control', disease_name],
        patch_artist=True,
        widths=0.5,
    )
    bp['boxes'][0].set_facecolor('#BBDEFB')
    bp['boxes'][1].set_facecolor('#FFCDD2')

    # Overlay points
    for i, data in enumerate([control_scores, disease_scores], 1):
        jitter = np.random.normal(0, 0.04, len(data))
        ax.scatter(np.full(len(data), i) + jitter, data,
                   alpha=0.4, s=15, color='black', zorder=3)

    ax.set_ylabel('MK4 Chaos Score', fontsize=12)
    ax.set_title(f'Group Comparison: {disease_name}', fontsize=13)
    ax.grid(True, alpha=0.3, axis='y')

    if save_path and created_fig:
        fig.tight_layout()
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {save_path}")
        plt.close(fig)
        return None

    if created_fig:
        return fig
    return None


def plot_roc_curve(roc_results, disease_name="Disease", ax=None, save_path=None):
    """
    Plot ROC curve with AUC.

    Parameters
    ----------
    roc_results : dict
        Output from find_optimal_threshold() containing 'fpr', 'tpr', 'auc'.
    disease_name : str
        Name of disease for title.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on.
    save_path : str, optional
        Path to save figure.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
        created_fig = True

    ax.plot(roc_results['fpr'], roc_results['tpr'],
            color='#F44336', linewidth=2,
            label=f"AUC = {roc_results['auc']:.3f}")
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random (AUC = 0.5)')

    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title(f'ROC Curve: {disease_name} Detection', fontsize=13)
    ax.legend(fontsize=11, loc='lower right')
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    if save_path and created_fig:
        fig.tight_layout()
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {save_path}")
        plt.close(fig)
        return None

    if created_fig:
        return fig
    return None


def plot_results_summary(control_scores, disease_scores, results, disease_name="Disease",
                         save_path=None):
    """
    Generate a 2x2 summary figure with all key plots.

    Parameters
    ----------
    control_scores : np.ndarray
        Chaos scores for control group.
    disease_scores : np.ndarray
        Chaos scores for disease group.
    results : dict
        Output from analyze_disease_dataset().
    disease_name : str
        Name of disease.
    save_path : str, optional
        Path to save figure.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Distribution
    plot_chaos_distribution(control_scores, disease_scores, disease_name, ax=axes[0, 0])

    # Boxplot
    plot_comparison_boxplot(control_scores, disease_scores, disease_name, ax=axes[0, 1])

    # ROC
    plot_roc_curve(results['roc'], disease_name, ax=axes[1, 0])

    # Metrics text
    ax = axes[1, 1]
    ax.axis('off')
    cl = results['classification']
    st = results['statistics']
    text = (
        f"{disease_name} Detection Results\n"
        f"{'='*35}\n\n"
        f"Accuracy:    {cl['accuracy']*100:.1f}%\n"
        f"Sensitivity: {cl['sensitivity']*100:.1f}%\n"
        f"Specificity: {cl['specificity']*100:.1f}%\n"
        f"PPV:         {cl['ppv']*100:.1f}%\n"
        f"NPV:         {cl['npv']*100:.1f}%\n"
        f"AUC:         {results['roc']['auc']:.3f}\n\n"
        f"p-value:     {st['p_value']:.2e}\n"
        f"Cohen's d:   {st['effect_size']:.3f}\n\n"
        f"Control:  n={len(control_scores)}\n"
        f"Disease:  n={len(disease_scores)}"
    )
    ax.text(0.1, 0.9, text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    fig.suptitle(f'MK4 Analysis: {disease_name}', fontsize=15, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {save_path}")
        plt.close(fig)

    return fig


def save_all_plots(control_scores, disease_scores, results, disease_name="Disease",
                   output_dir="results/figures"):
    """
    Generate and save all standard plots for a disease analysis.

    Parameters
    ----------
    control_scores : np.ndarray
        Chaos scores for control group.
    disease_scores : np.ndarray
        Chaos scores for disease group.
    results : dict
        Output from analyze_disease_dataset().
    disease_name : str
        Name of disease.
    output_dir : str
        Directory to save figures.
    """
    os.makedirs(output_dir, exist_ok=True)
    tag = disease_name.lower().replace(' ', '_')

    plot_chaos_distribution(
        control_scores, disease_scores, disease_name,
        save_path=os.path.join(output_dir, f'{tag}_distribution.png')
    )

    plot_comparison_boxplot(
        control_scores, disease_scores, disease_name,
        save_path=os.path.join(output_dir, f'{tag}_boxplot.png')
    )

    plot_roc_curve(
        results['roc'], disease_name,
        save_path=os.path.join(output_dir, f'{tag}_roc.png')
    )

    plot_results_summary(
        control_scores, disease_scores, results, disease_name,
        save_path=os.path.join(output_dir, f'{tag}_summary.png')
    )

    print(f"  All figures saved to {output_dir}/")
