#!/usr/bin/env python3
"""
mk4.eeg.figures — Publication figures for MK4-EEG paper.

Figures:
  1. Cosine similarity heatmap (disease orientation vectors)
  2. Classification performance (pre/post filter + CV) with bootstrap CI
  3. Feature distributions per disease (top-3 discriminative)
  4. Pipeline schematic
  5. PCA subject embedding + cosine graph overlay

Usage:
  export MK4_EEG_DATA=/path/to/results
  export MK4_EEG_FIGURES=/path/to/figures
  python -m mk4.eeg.figures
"""
import os

EEG_DIR = os.environ.get("MK4_EEG_DATA", "./results/eeg")
FIG_DIR = os.environ.get("MK4_EEG_FIGURES", "./results/eeg/figures")

os.makedirs(FIG_DIR, exist_ok=True)

import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

os.makedirs(FIG_DIR, exist_ok=True)

FEATS = ["band_delta", "band_theta", "band_alpha", "band_beta", "band_gamma",
         "alpha_theta_ratio", "spectral_slope", "coherence", "peak_ratio", "hurst"]

FEAT_SHORT = ["δ", "θ", "α", "β", "γ", "α/θ", "1/f", "coh", "peak", "H"]

DATASETS = [
    ("ds003478", "MDD", "Depression"),
    ("ds004504", "AD+FTD", "Alzheimer/FTD"),
    ("ds002778", "PD", "Parkinson"),
    ("ds003523", "TBI", "Mild TBI"),
]


def save_both(fig, name):
    fig.savefig(f"{FIG_DIR}/{name}.png", dpi=300, bbox_inches='tight')
    fig.savefig(f"{FIG_DIR}/{name}.pdf", bbox_inches='tight')
    print(f"  saved: {FIG_DIR}/{name}.{{png,pdf}}")


def get_trained_weights(dataset_id):
    """Pull trained weights from stats JSON full_data_fit (if present)
    or from base mk4e_results JSON trained_weights (dict form)."""
    stats_path = f"{EEG_DIR}/mk4e_stats_{dataset_id}.json"
    if os.path.exists(stats_path):
        s = json.load(open(stats_path))
        if 'full_data_fit' in s:
            return np.array(s['full_data_fit']['weights'])
    # Fallback: from mk4e_results
    r = json.load(open(f"{EEG_DIR}/mk4e_results_{dataset_id}.json"))
    tw = r.get('trained_weights', {})
    w = np.zeros(len(FEATS))
    for i, f in enumerate(FEATS):
        w[i] = tw.get(f, 0)
    return w


# ── Fig 1: Cosine similarity heatmap ──

def fig1_cosine_heatmap():
    print("\nFig 1: cosine similarity heatmap")
    W = np.stack([get_trained_weights(ds[0]) for ds in DATASETS])
    W_norm = W / (np.linalg.norm(W, axis=1, keepdims=True) + 1e-12)
    S = W_norm @ W_norm.T
    labels = [ds[1] for ds in DATASETS]

    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    # Blue-white-red diverging colormap
    cmap = plt.get_cmap('RdBu_r')
    im = ax.imshow(S, cmap=cmap, vmin=-1, vmax=1)

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_yticklabels(labels, fontsize=11)

    for i in range(len(labels)):
        for j in range(len(labels)):
            val = S[i, j]
            color = 'white' if abs(val) > 0.5 else 'black'
            ax.text(j, i, f"{val:+.2f}", ha='center', va='center',
                    color=color, fontsize=10, fontweight='bold')

    ax.set_title('Cosine similarity of learned weight vectors\n(spectral disease fingerprints)',
                 fontsize=12, pad=10)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Cosine similarity', fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    save_both(fig, "fig1_cosine_heatmap")
    plt.close(fig)
    return S


# ── Fig 2: Classification performance pre/post filter ──

def fig2_performance_barchart():
    print("\nFig 2: classification performance pre/post filter")
    labels, j_pre, j_post, j_cv_mean, j_cv_std, ci_lo, ci_hi = [], [], [], [], [], [], []
    for ds_id, short, _ in DATASETS:
        r = json.load(open(f"{EEG_DIR}/mk4e_results_{ds_id}.json"))
        analyses = r.get('analyses', {})
        all_analysis = analyses.get('CTRL_vs_ALL', {})
        if not all_analysis:
            # Fallback
            b = r.get('before', {}); a = r.get('after', {}) or b
        else:
            b = all_analysis.get('before', {}); a = all_analysis.get('after', {}) or b
        j_pre.append(b.get('youden_j', 0))
        j_post.append(a.get('youden_j', 0))
        labels.append(short)

        stats_path = f"{EEG_DIR}/mk4e_stats_{ds_id}.json"
        if os.path.exists(stats_path):
            s = json.load(open(stats_path))
            if 'bootstrap' in s:
                b_stats = s['bootstrap']
                ci_lo.append(b_stats['j_ci_lower'])
                ci_hi.append(b_stats['j_ci_upper'])
            else:
                ci_lo.append(a.get('youden_j', 0))
                ci_hi.append(a.get('youden_j', 0))
            if 'cv' in s:
                agg = s['cv']['aggregate']
                j_cv_mean.append(agg['j_mean'])
                j_cv_std.append(agg['j_std'])
            else:
                j_cv_mean.append(None)
                j_cv_std.append(None)
        else:
            ci_lo.append(a.get('youden_j', 0))
            ci_hi.append(a.get('youden_j', 0))
            j_cv_mean.append(None)
            j_cv_std.append(None)

    x = np.arange(len(labels))
    width = 0.28

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.bar(x - width, j_pre, width, label='Pre-filter', color='#b0c4de', edgecolor='k', linewidth=0.8)

    # Post-filter with bootstrap CI as error bars
    err_lo = [p - l for p, l in zip(j_post, ci_lo)]
    err_hi = [u - p for p, u in zip(j_post, ci_hi)]
    ax.bar(x, j_post, width, label='Post-filter (95% CI)', color='#4a8fc2', edgecolor='k', linewidth=0.8,
           yerr=[err_lo, err_hi], capsize=4, error_kw={'ecolor': 'k', 'elinewidth': 1})

    # CV J as separate bars where available
    cv_x, cv_mean, cv_std = [], [], []
    for i, (m, s) in enumerate(zip(j_cv_mean, j_cv_std)):
        if m is not None:
            cv_x.append(x[i] + width)
            cv_mean.append(m)
            cv_std.append(s or 0)
    if cv_mean:
        ax.bar(cv_x, cv_mean, width, label='5-fold CV (train-only sanitize)',
               color='#d89a50', edgecolor='k', linewidth=0.8,
               yerr=cv_std, capsize=4, error_kw={'ecolor': 'k', 'elinewidth': 1})

    ax.axhline(0, color='k', linewidth=0.6)
    ax.axhline(0.5, color='gray', linestyle=':', linewidth=0.6, alpha=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("Youden's J", fontsize=11)
    ax.set_title('Classification performance: pre-filter vs. post-filter vs. 5-fold CV',
                 fontsize=11, pad=10)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(axis='y', alpha=0.3, linestyle=':')

    save_both(fig, "fig2_performance")
    plt.close(fig)


# ── Fig 3: Feature distributions per disease ──

def fig3_feature_distributions():
    print("\nFig 3: feature distributions (top 3 discriminative per disease)")
    fig, axes = plt.subplots(len(DATASETS), 3, figsize=(10, 10))

    for row, (ds_id, short, full) in enumerate(DATASETS):
        r = json.load(open(f"{EEG_DIR}/mk4e_results_{ds_id}.json"))
        subjects = r.get('subjects', [])
        # Build feature matrix
        X, y = [], []
        for s in subjects:
            if s.get('label') not in ('CONTROL', 'DISEASE'):
                continue
            try:
                X.append([float(s[k]) for k in FEATS])
                y.append(1 if s['label'] == 'DISEASE' else 0)
            except Exception:
                continue
        X = np.array(X)
        y = np.array(y)

        # Use trained weights magnitude as "importance"
        w = get_trained_weights(ds_id)
        top3 = np.argsort(-np.abs(w))[:3]

        for col, feat_idx in enumerate(top3):
            ax = axes[row, col]
            c_vals = X[y == 0, feat_idx]
            d_vals = X[y == 1, feat_idx]
            bins = np.linspace(min(c_vals.min(), d_vals.min()),
                               max(c_vals.max(), d_vals.max()), 20)
            ax.hist(c_vals, bins=bins, alpha=0.5, color='#4a8fc2', label='CTRL', density=True)
            ax.hist(d_vals, bins=bins, alpha=0.5, color='#c24a4a', label=short, density=True)
            ax.set_title(f"{FEAT_SHORT[feat_idx]} (w={w[feat_idx]:+.1f})", fontsize=10)
            ax.tick_params(labelsize=8)
            if col == 0:
                ax.set_ylabel(full, fontsize=10)
            if row == 0 and col == 2:
                ax.legend(fontsize=8, loc='best')

    fig.suptitle('Top 3 discriminative features per diagnosis (by |weight|)',
                 fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    save_both(fig, "fig3_features")
    plt.close(fig)


# ── Fig 4: Pipeline schematic (simple) ──

def fig4_pipeline():
    print("\nFig 4: pipeline schematic")
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.axis('off')

    stages = [
        (0.04, "Raw EEG\n(.set/.bdf/.edf)", '#aec8e0'),
        (0.22, "Preprocessing\n(0.5-45 Hz, notch)", '#aec8e0'),
        (0.40, "Welch PSD\n(per channel, avg)", '#aec8e0'),
        (0.60, "10 features\n(δ/θ/α/β/γ +\nH, coh, slope...)", '#d0a840'),
        (0.80, "MK4E score\nΣ wᵢ · featᵢ", '#c2a04a'),
        (0.96, "Youden\nthreshold", '#c26060'),
    ]
    for x, label, color in stages:
        ax.add_patch(plt.Rectangle((x - 0.07, 0.55), 0.14, 0.3, facecolor=color,
                                   edgecolor='black', linewidth=1.2))
        ax.text(x, 0.70, label, ha='center', va='center', fontsize=9)

    # Arrows between stages
    for i in range(len(stages) - 1):
        x1 = stages[i][0] + 0.07
        x2 = stages[i + 1][0] - 0.07
        ax.annotate('', xy=(x2, 0.70), xytext=(x1, 0.70),
                    arrowprops=dict(arrowstyle='->', lw=1.2, color='black'))

    # CTRL sanitize side-branch
    ax.add_patch(plt.Rectangle((0.53, 0.15), 0.14, 0.22,
                               facecolor='#d08a50', edgecolor='black', linewidth=1.2))
    ax.text(0.60, 0.26, "CTRL sanitize\nmedian ± MAD\n(strict/soft/bias)",
            ha='center', va='center', fontsize=8)

    ax.annotate('', xy=(0.60, 0.37), xytext=(0.60, 0.55),
                arrowprops=dict(arrowstyle='->', lw=1.0, color='black', linestyle='--'))
    ax.text(0.68, 0.46, "flag suspect\nCTRL samples", fontsize=8, style='italic')

    # Weight training loop
    ax.add_patch(plt.Rectangle((0.73, 0.15), 0.14, 0.22,
                               facecolor='#9fc2a0', edgecolor='black', linewidth=1.2))
    ax.text(0.80, 0.26, "Grid search\nmaximize J\n(train only)",
            ha='center', va='center', fontsize=8)
    ax.annotate('', xy=(0.80, 0.37), xytext=(0.80, 0.55),
                arrowprops=dict(arrowstyle='->', lw=1.0, color='black', linestyle='--'))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("MK4-EEG pipeline: feature extraction → weight learning → Youden classification",
                 fontsize=11, pad=10)

    save_both(fig, "fig4_pipeline")
    plt.close(fig)


# ── Fig 5: PCA embedding of subjects + cosine overlay ──

def fig5_pca_embedding():
    print("\nFig 5: PCA of subject feature vectors + disease centroid cosine overlay")
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Collect per-subject feature vectors from all 4 datasets, labeled
    X_all, meta_all = [], []
    for ds_id, short, _ in DATASETS:
        r = json.load(open(f"{EEG_DIR}/mk4e_results_{ds_id}.json"))
        for s in r.get('subjects', []):
            lbl = s.get('label')
            if lbl not in ('CONTROL', 'DISEASE'):
                continue
            try:
                X_all.append([float(s[k]) for k in FEATS])
                meta_all.append((ds_id, short, lbl))
            except Exception:
                continue
    X_all = np.array(X_all)
    # Standardize then PCA to 2D
    X_std = StandardScaler().fit_transform(X_all)
    pca = PCA(n_components=2, random_state=42)
    XY = pca.fit_transform(X_std)
    var_explained = pca.explained_variance_ratio_

    # Disease centroids in PCA space (DIS subjects only per dataset)
    centroids = {}
    for ds_id, short, _ in DATASETS:
        mask = np.array([m[0] == ds_id and m[2] == 'DISEASE' for m in meta_all])
        if mask.any():
            centroids[short] = XY[mask].mean(axis=0)

    # Cosine similarity from weights (for overlay)
    W = np.stack([get_trained_weights(ds[0]) for ds in DATASETS])
    Wn = W / (np.linalg.norm(W, axis=1, keepdims=True) + 1e-12)
    S = Wn @ Wn.T
    labels = [ds[1] for ds in DATASETS]

    fig, (ax_pca, ax_graph) = plt.subplots(1, 2, figsize=(13, 5.5))

    # --- Left panel: PCA scatter ---
    colors = {'MDD': '#c24a4a', 'AD+FTD': '#d89a50', 'PD': '#4a8fc2', 'TBI': '#6a6a6a'}
    for ds_id, short, _ in DATASETS:
        for lbl, alpha, marker in [('CONTROL', 0.35, 'o'), ('DISEASE', 0.85, 's')]:
            mask = np.array([m[0] == ds_id and m[2] == lbl for m in meta_all])
            if mask.any():
                ax_pca.scatter(XY[mask, 0], XY[mask, 1],
                               c=colors[short], alpha=alpha,
                               marker=marker, s=28 if lbl == 'DISEASE' else 18,
                               edgecolor='black', linewidth=0.3,
                               label=f"{short} {lbl}" if lbl == 'DISEASE' else None)

    # Disease centroids with labels
    for short, c in centroids.items():
        ax_pca.scatter(c[0], c[1], c=colors[short], s=220, marker='*',
                       edgecolor='black', linewidth=1.2, zorder=5)
        ax_pca.annotate(short, c, xytext=(8, 8), textcoords='offset points',
                        fontsize=10, fontweight='bold')

    ax_pca.set_xlabel(f'PC1 ({var_explained[0]*100:.1f}%)', fontsize=10)
    ax_pca.set_ylabel(f'PC2 ({var_explained[1]*100:.1f}%)', fontsize=10)
    ax_pca.set_title('Subject feature vectors in shared PCA space\n(CTRL = faded circles, DIS = dark squares, ★ = disease centroid)',
                     fontsize=10.5, pad=8)
    ax_pca.legend(loc='best', fontsize=8, framealpha=0.9)
    ax_pca.grid(alpha=0.3, linestyle=':')
    ax_pca.axhline(0, color='k', linewidth=0.4, alpha=0.4)
    ax_pca.axvline(0, color='k', linewidth=0.4, alpha=0.4)

    # --- Right panel: cosine graph overlay ---
    ax_graph.set_xlim(-1.3, 1.3)
    ax_graph.set_ylim(-1.3, 1.3)
    ax_graph.set_aspect('equal')
    ax_graph.axis('off')

    # Place 4 nodes on a square
    pos = {'MDD': (1.0, 0.0), 'AD+FTD': (0.0, 1.0),
           'PD': (-1.0, 0.0), 'TBI': (0.0, -1.0)}

    # Edges (colored by cosine sign/magnitude)
    for i, l1 in enumerate(labels):
        for j, l2 in enumerate(labels):
            if j <= i:
                continue
            cos = S[i, j]
            # Color: red positive, blue negative
            color = plt.cm.RdBu_r((cos + 1) / 2)
            lw = 0.5 + 4.5 * abs(cos)
            ax_graph.plot([pos[l1][0], pos[l2][0]],
                          [pos[l1][1], pos[l2][1]],
                          color=color, linewidth=lw, alpha=0.85, zorder=1)
            mx, my = (pos[l1][0] + pos[l2][0]) / 2, (pos[l1][1] + pos[l2][1]) / 2
            ax_graph.text(mx, my, f"{cos:+.2f}", fontsize=9, ha='center', va='center',
                          bbox=dict(facecolor='white', edgecolor='none',
                                    boxstyle='round,pad=0.15', alpha=0.9), zorder=3)

    # Nodes
    for lbl, (x, y) in pos.items():
        ax_graph.scatter(x, y, s=1200, c=colors[lbl], edgecolor='black',
                         linewidth=1.5, zorder=2)
        ax_graph.text(x, y, lbl, ha='center', va='center',
                      fontsize=10, fontweight='bold', color='white', zorder=4)

    ax_graph.set_title('Orientation-vector cosine graph\n(edge color & width ∝ similarity)',
                       fontsize=10.5, pad=8)

    fig.suptitle('Shared spectral embedding + disease-vector cosine geometry',
                 fontsize=12, y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_both(fig, "fig5_embedding")
    plt.close(fig)


if __name__ == "__main__":
    print(f"MK4-EEG figures → {FIG_DIR}/")
    fig1_cosine_heatmap()
    fig2_performance_barchart()
    fig3_feature_distributions()
    fig4_pipeline()
    fig5_pca_embedding()
    print("\nDone.")
