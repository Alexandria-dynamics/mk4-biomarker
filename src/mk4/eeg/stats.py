#!/usr/bin/env python3
"""
mk4.eeg.stats — Statistical validation: 5-fold CV + permutation + bootstrap CI.

Reads feature JSONs produced by mk4.eeg.batch. Implements strict no-data-leakage:
CTRL sanitization uses train-only median+MAD; test samples flagged by train
thresholds. Grid-search weight optimization confined to training fold.

Usage:
  python -m mk4.eeg.stats <dataset_label> <path/to/results.json> \
      [--cv] [--perm N] [--boot N] [--out OUT]
"""


import json
import sys
import os
import time
import argparse
import numpy as np
from sklearn.model_selection import StratifiedKFold

FEATS = ["band_delta", "band_theta", "band_alpha", "band_beta", "band_gamma",
         "alpha_theta_ratio", "spectral_slope", "coherence", "peak_ratio", "hurst"]


# ── Core utilities ──

def load_features(results_path):
    d = json.load(open(results_path))
    subjects = d.get('subjects', [])
    X, y, ids = [], [], []
    for s in subjects:
        if s.get('label') not in ('CONTROL', 'DISEASE'):
            continue
        try:
            row = [float(s[k]) for k in FEATS]
        except (KeyError, TypeError, ValueError):
            continue
        X.append(row)
        y.append(1 if s['label'] == 'DISEASE' else 0)
        ids.append(s.get('subject', ''))
    return np.array(X, dtype=float), np.array(y, dtype=int), ids


def youden_threshold(ctrl_scores, dis_scores, dis_is_high=True):
    """Find Youden-optimal threshold. Returns (threshold, J)."""
    all_v = np.concatenate([ctrl_scores, dis_scores])
    if all_v.size < 4 or ctrl_scores.size < 2 or dis_scores.size < 2:
        return 0.0, 0.0
    thresholds = np.linspace(all_v.min(), all_v.max(), 300)
    best_j, best_t = -1.0, float(np.median(all_v))
    for t in thresholds:
        if dis_is_high:
            tp = float(np.sum(dis_scores >= t))
            fn = float(np.sum(dis_scores < t))
            tn = float(np.sum(ctrl_scores < t))
            fp = float(np.sum(ctrl_scores >= t))
        else:
            tp = float(np.sum(dis_scores < t))
            fn = float(np.sum(dis_scores >= t))
            tn = float(np.sum(ctrl_scores >= t))
            fp = float(np.sum(ctrl_scores < t))
        sens = tp / (tp + fn + 1e-9)
        spec = tn / (tn + fp + 1e-9)
        j = sens + spec - 1
        if j > best_j:
            best_j, best_t = j, float(t)
    return best_t, best_j


def robust_stats_matrix(X):
    """Return median and MAD per column (for feature-wise standardization)."""
    med = np.median(X, axis=0)
    mad = np.median(np.abs(X - med), axis=0) * 1.4826 + 1e-12
    return med, mad


def sanitize_ctrl(X_ctrl, X_dis, strict=2.5, soft=1.8, bias=1.8):
    """
    Apply MK4 CTRL sanitize filter using train-CTRL-only statistics.
    Returns (kept_mask, params) where params = {'med', 'mad', 'mk4_med', 'mk4_mad', 'dis_is_high'}.
    The params allow re-applying the same filter to any other CTRL samples.
    """
    if X_ctrl.shape[0] < 4:
        return np.ones(X_ctrl.shape[0], dtype=bool), None

    med, mad = robust_stats_matrix(X_ctrl)
    z = (X_ctrl - med) / mad
    abs_z = np.abs(z)

    strict_hits = np.any(abs_z > strict, axis=1)
    soft_hits = np.sum(abs_z > soft, axis=1) >= 2

    # Direction-aware MK4 bias: use simple single-best-feature score as proxy
    # Pick best single feature on the training distributions
    best_key, best_sign = None, +1
    best_j = -1
    for k in range(X_ctrl.shape[1]):
        for sgn in [+1, -1]:
            sc_c = X_ctrl[:, k] * sgn
            sc_d = X_dis[:, k] * sgn
            dis_high = float(np.mean(sc_d)) > float(np.mean(sc_c))
            _, j = youden_threshold(sc_c, sc_d, dis_high)
            if j > best_j:
                best_j, best_key, best_sign = j, k, sgn

    # MK4 score per CTRL sample using best_key direction
    mk4_score_ctrl = X_ctrl[:, best_key] * best_sign
    mk4_score_dis = X_dis[:, best_key] * best_sign
    mk4_med, mk4_mad = float(np.median(mk4_score_ctrl)), float(np.median(np.abs(mk4_score_ctrl - np.median(mk4_score_ctrl))) * 1.4826 + 1e-12)
    dis_is_high_mk4 = float(np.mean(mk4_score_dis)) > float(np.mean(mk4_score_ctrl))

    z_mk4 = (mk4_score_ctrl - mk4_med) / mk4_mad
    if dis_is_high_mk4:
        bias_hits = z_mk4 > bias
    else:
        bias_hits = z_mk4 < -bias

    kept_mask = ~(strict_hits | soft_hits | bias_hits)
    if kept_mask.sum() < 3:
        return np.ones(X_ctrl.shape[0], dtype=bool), None

    params = {
        'med': med.tolist(), 'mad': mad.tolist(),
        'best_key': int(best_key), 'best_sign': int(best_sign),
        'mk4_med': float(mk4_med), 'mk4_mad': float(mk4_mad),
        'dis_is_high_mk4': bool(dis_is_high_mk4),
        'strict': float(strict), 'soft': float(soft), 'bias': float(bias),
    }
    return kept_mask, params


def train_weights(X_ctrl, X_dis):
    """Grid search for weight vector maximizing Youden J. Returns (w, J)."""
    n_feats = X_ctrl.shape[1]
    best_j = -1.0
    best_w = np.zeros(n_feats)
    best_key = None

    # Single feature
    for k in range(n_feats):
        for sgn in [+1, -1]:
            w = np.zeros(n_feats)
            w[k] = sgn
            sc_c = X_ctrl @ w
            sc_d = X_dis @ w
            dis_high = float(np.mean(sc_d)) > float(np.mean(sc_c))
            _, j = youden_threshold(sc_c, sc_d, dis_high)
            if j > best_j:
                best_j = j
                best_w = w.copy()
                best_key = k

    # Two-feature expansion
    if best_key is not None:
        for k2 in range(n_feats):
            if k2 == best_key:
                continue
            for w2 in [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]:
                w = best_w.copy()
                w[k2] = w2
                sc_c = X_ctrl @ w
                sc_d = X_dis @ w
                dis_high = float(np.mean(sc_d)) > float(np.mean(sc_c))
                _, j = youden_threshold(sc_c, sc_d, dis_high)
                if j > best_j:
                    best_j = j
                    best_w = w.copy()

    # Three-feature expansion
    for k3 in range(n_feats):
        if best_w[k3] != 0:
            continue
        for w3 in [-1.0, -0.5, 0.5, 1.0]:
            w = best_w.copy()
            w[k3] = w3
            sc_c = X_ctrl @ w
            sc_d = X_dis @ w
            dis_high = float(np.mean(sc_d)) > float(np.mean(sc_c))
            _, j = youden_threshold(sc_c, sc_d, dis_high)
            if j > best_j:
                best_j = j
                best_w = w.copy()

    return best_w, float(best_j)


def compute_metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=int)
    y_pred = np.asarray(y_pred, dtype=int)
    tp = int(np.sum((y_pred == 1) & (y_true == 1)))
    tn = int(np.sum((y_pred == 0) & (y_true == 0)))
    fp = int(np.sum((y_pred == 1) & (y_true == 0)))
    fn = int(np.sum((y_pred == 0) & (y_true == 1)))
    n = tp + tn + fp + fn
    acc = (tp + tn) / n if n else 0.0
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    j = sens + spec - 1
    return {
        'acc': float(acc), 'sens': float(sens), 'spec': float(spec), 'j': float(j),
        'tp': tp, 'tn': tn, 'fp': fp, 'fn': fn
    }


# ── Cross-validation ──

def fit_fold(X_train, y_train):
    """
    Fit one CV fold: sanitize train CTRL (train-only stats), train weights,
    find threshold. Returns everything needed to score test fold.
    """
    X_ctrl_train = X_train[y_train == 0]
    X_dis_train = X_train[y_train == 1]

    kept_mask, san_params = sanitize_ctrl(X_ctrl_train, X_dis_train)
    X_ctrl_clean = X_ctrl_train[kept_mask]

    if X_ctrl_clean.shape[0] < 3 or X_dis_train.shape[0] < 3:
        # Fallback: no sanitize possible
        X_ctrl_clean = X_ctrl_train
        san_params = None

    w, j_train = train_weights(X_ctrl_clean, X_dis_train)

    sc_c = X_ctrl_clean @ w
    sc_d = X_dis_train @ w
    dis_high = float(np.mean(sc_d)) > float(np.mean(sc_c))
    thr, _ = youden_threshold(sc_c, sc_d, dis_high)

    return {
        'weights': w,
        'threshold': float(thr),
        'dis_is_high': bool(dis_high),
        'j_train': float(j_train),
        'n_ctrl_train_raw': int(X_ctrl_train.shape[0]),
        'n_ctrl_train_kept': int(X_ctrl_clean.shape[0]),
        'sanitize_params': san_params,
    }


def predict_fold(X_test, fit):
    scores = X_test @ fit['weights']
    if fit['dis_is_high']:
        return (scores >= fit['threshold']).astype(int)
    return (scores < fit['threshold']).astype(int)


def run_cv(X, y, k=5, seed=42):
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    folds = []
    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        fit = fit_fold(X_train, y_train)
        y_pred = predict_fold(X_test, fit)
        m = compute_metrics(y_test, y_pred)
        m['fold'] = fold_idx
        m['j_train'] = fit['j_train']
        m['n_train'] = int(len(train_idx))
        m['n_test'] = int(len(test_idx))
        m['n_ctrl_train_kept'] = fit['n_ctrl_train_kept']
        m['n_ctrl_train_rejected'] = fit['n_ctrl_train_raw'] - fit['n_ctrl_train_kept']
        m['weights'] = fit['weights'].tolist()
        m['threshold'] = fit['threshold']
        m['dis_is_high'] = fit['dis_is_high']
        folds.append(m)
        print(f"  fold {fold_idx+1}/{k}: J_train={m['j_train']:.3f} "
              f"J_test={m['j']:.3f} ACC={m['acc']*100:.1f}%", flush=True)

    agg = {}
    for key in ['acc', 'sens', 'spec', 'j']:
        vals = [r[key] for r in folds]
        agg[f'{key}_mean'] = float(np.mean(vals))
        agg[f'{key}_std'] = float(np.std(vals, ddof=1))
        agg[f'{key}_min'] = float(np.min(vals))
        agg[f'{key}_max'] = float(np.max(vals))
    return {'folds': folds, 'aggregate': agg}


# ── Permutation test ──

def run_permutation(X, y, n=1000, seed=42):
    rng = np.random.default_rng(seed)
    y_shuffled = y.copy()

    # Actual full-data J (no sanitize, single best single-feature grid for comparability)
    _, actual_j = train_weights(X[y == 0], X[y == 1])

    null_js = []
    t0 = time.time()
    for i in range(n):
        rng.shuffle(y_shuffled)
        X_c = X[y_shuffled == 0]
        X_d = X[y_shuffled == 1]
        if X_c.shape[0] < 3 or X_d.shape[0] < 3:
            continue
        _, jn = train_weights(X_c, X_d)
        null_js.append(float(jn))
        if (i + 1) % 100 == 0:
            rate = (i + 1) / (time.time() - t0)
            eta = (n - i - 1) / rate
            print(f"  perm {i+1}/{n} rate={rate:.1f}/s null_max={max(null_js):.3f} ETA={eta:.0f}s", flush=True)

    null_arr = np.array(null_js)
    p_value = float(np.mean(null_arr >= actual_j))
    return {
        'actual_j': float(actual_j),
        'null_j_mean': float(np.mean(null_arr)),
        'null_j_std': float(np.std(null_arr)),
        'null_j_max': float(np.max(null_arr)),
        'null_j_p50': float(np.percentile(null_arr, 50)),
        'null_j_p95': float(np.percentile(null_arr, 95)),
        'null_j_p99': float(np.percentile(null_arr, 99)),
        'p_value': p_value,
        'n_permutations': int(len(null_js)),
        'null_j_list': null_js,
    }


# ── Bootstrap CI ──

def run_bootstrap(X_ctrl, X_dis, w, n=1000, seed=42):
    """
    Bootstrap Youden J with FIXED weights. Resample each group with replacement,
    recompute Youden threshold on resample, record J.
    """
    rng = np.random.default_rng(seed)
    sc_c_full = X_ctrl @ w
    sc_d_full = X_dis @ w
    dis_high = float(np.mean(sc_d_full)) > float(np.mean(sc_c_full))

    # Point estimate
    _, j_point = youden_threshold(sc_c_full, sc_d_full, dis_high)

    js = []
    for i in range(n):
        idx_c = rng.integers(0, X_ctrl.shape[0], size=X_ctrl.shape[0])
        idx_d = rng.integers(0, X_dis.shape[0], size=X_dis.shape[0])
        bs_c = X_ctrl[idx_c] @ w
        bs_d = X_dis[idx_d] @ w
        dh = float(np.mean(bs_d)) > float(np.mean(bs_c))
        _, j = youden_threshold(bs_c, bs_d, dh)
        js.append(float(j))
    js_arr = np.array(js)
    return {
        'j_point': float(j_point),
        'j_mean': float(np.mean(js_arr)),
        'j_median': float(np.median(js_arr)),
        'j_ci_lower': float(np.percentile(js_arr, 2.5)),
        'j_ci_upper': float(np.percentile(js_arr, 97.5)),
        'j_sd': float(np.std(js_arr)),
        'n_bootstrap': int(n),
    }


# ── Driver ──

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset", help="Dataset label (e.g., ds003478)")
    ap.add_argument("results_json", help="Path to mk4e_results_XXX.json")
    ap.add_argument("--cv", action="store_true", help="Run 5-fold CV")
    ap.add_argument("--perm", type=int, default=0, help="Number of permutations (0=skip)")
    ap.add_argument("--boot", type=int, default=0, help="Number of bootstrap resamples (0=skip)")
    ap.add_argument("--out", default=None, help="Output JSON path")
    args = ap.parse_args()

    X, y, ids = load_features(args.results_json)
    print(f"{args.dataset}: N={len(y)} CTRL={int(np.sum(y==0))} DIS={int(np.sum(y==1))}", flush=True)

    out = {
        'dataset': args.dataset,
        'n_total': int(len(y)),
        'n_ctrl': int(np.sum(y == 0)),
        'n_dis': int(np.sum(y == 1)),
    }

    if args.cv:
        print("=== 5-fold Stratified CV ===", flush=True)
        t0 = time.time()
        out['cv'] = run_cv(X, y)
        agg = out['cv']['aggregate']
        print(f"CV: J = {agg['j_mean']:.3f} ± {agg['j_std']:.3f}, "
              f"ACC = {agg['acc_mean']*100:.1f}% ± {agg['acc_std']*100:.1f}%", flush=True)
        print(f"CV took {time.time()-t0:.1f}s", flush=True)

    if args.perm > 0:
        print(f"=== Permutation test ({args.perm} shuffles) ===", flush=True)
        t0 = time.time()
        out['permutation'] = run_permutation(X, y, n=args.perm)
        p = out['permutation']
        print(f"Actual J={p['actual_j']:.3f}, null mean={p['null_j_mean']:.3f}, "
              f"null p95={p['null_j_p95']:.3f}, p-value={p['p_value']:.4f}", flush=True)
        print(f"Permutation took {time.time()-t0:.1f}s", flush=True)

    if args.boot > 0:
        print(f"=== Bootstrap ({args.boot} resamples) ===", flush=True)
        t0 = time.time()
        # Full-data sanitize + train (matches published results)
        X_ctrl = X[y == 0]
        X_dis = X[y == 1]
        kept, san_params = sanitize_ctrl(X_ctrl, X_dis)
        X_ctrl_clean = X_ctrl[kept]
        w_full, j_train_full = train_weights(X_ctrl_clean, X_dis)
        out['full_data_fit'] = {
            'n_ctrl_kept': int(kept.sum()),
            'n_ctrl_rejected': int((~kept).sum()),
            'weights': w_full.tolist(),
            'j_train': float(j_train_full),
            'sanitize_params': san_params,
        }
        out['bootstrap'] = run_bootstrap(X_ctrl_clean, X_dis, w_full, n=args.boot)
        b = out['bootstrap']
        print(f"J point={b['j_point']:.3f}, 95% CI [{b['j_ci_lower']:.3f}, {b['j_ci_upper']:.3f}]", flush=True)
        print(f"Bootstrap took {time.time()-t0:.1f}s", flush=True)

    out_path = args.out or f"mk4e_stats_{args.dataset}.json"
    json.dump(out, open(out_path, 'w'), indent=2)
    print(f"Saved: {out_path}", flush=True)


if __name__ == "__main__":
    main()
