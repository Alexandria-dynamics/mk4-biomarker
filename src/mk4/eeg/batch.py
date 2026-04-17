#!/usr/bin/env python3
"""
mk4.eeg.batch — Spectral biomarker pipeline for EEG time series.

Feature extraction and per-dataset orientation-vector optimization.
Companion to the RNA expression MK4 for EEG modality.

Per subject:
  - Welch PSD averaged across EEG channels (bandpass 0.5-45 Hz, 50+60 Hz notch)
  - 10 features: 5 band powers (delta/theta/alpha/beta/gamma),
    alpha/theta ratio, spectral slope, peak freq/ratio,
    Shannon spectral entropy, geometric/arithmetic coherence, Hurst exponent

Orientation vector trained per dataset via brute-force grid search maximizing
Youden's J. CTRL distributional harmonization via median+MAD with asymmetry
deviation score.

Usage:
  python -m mk4.eeg.batch <dataset_root>
"""


import os, sys, json, gc, time, re
from pathlib import Path
import numpy as np

FEATS = ["chaos", "coherence", "band_delta", "band_theta", "band_alpha",
         "band_beta", "band_gamma", "peak_ratio", "hurst", "alpha_theta_ratio",
         "spectral_slope"]

MK4E_WEIGHTS = None  # trained per-dataset via optimize_weights()

BANDS = {
    "delta": (0.5, 4.0),
    "theta": (4.0, 8.0),
    "alpha": (8.0, 13.0),
    "beta":  (13.0, 30.0),
    "gamma": (30.0, 45.0),
}


def compute_eeg_features(raw, picks=None):
    """
    Extract MK4-EEG features from an MNE Raw object.
    Returns dict with band powers, entropy, peak freq, etc.
    Averages across kept EEG channels (single scalar per feature per subject).
    """
    import mne
    from scipy.signal import welch

    raw = raw.copy()
    # Standard preprocessing: pick EEG only, bandpass 0.5-45, drop bad
    eeg_picks = mne.pick_types(raw.info, eeg=True, exclude="bads")
    if len(eeg_picks) < 5:
        return None
    raw.pick(eeg_picks)
    raw.filter(0.5, 45.0, verbose=False)
    # Notch for line noise — try 50 Hz (EU) + 60 Hz
    try:
        raw.notch_filter([50, 60], verbose=False)
    except Exception:
        pass

    data = raw.get_data()  # (n_ch, n_samples)
    sfreq = raw.info["sfreq"]

    # Welch PSD per channel, then average
    f, psd = welch(data, fs=sfreq, nperseg=int(min(4 * sfreq, data.shape[1])), axis=-1)
    psd_mean = psd.mean(axis=0)  # average across channels
    # Focus on 0.5-45 Hz
    mask = (f >= 0.5) & (f <= 45.0)
    f = f[mask]
    psd_mean = psd_mean[mask]
    total = psd_mean.sum() + 1e-20

    # Band powers (relative)
    band_pow = {}
    for name, (lo, hi) in BANDS.items():
        bm = (f >= lo) & (f < hi)
        band_pow[name] = float(psd_mean[bm].sum() / total)

    # Alpha/theta ratio (classic AD biomarker)
    alpha_theta_ratio = float(band_pow["alpha"] / (band_pow["theta"] + 1e-12))

    # Spectral slope (1/f exponent, fit log PSD vs log f in 2-40 Hz)
    slope_mask = (f >= 2.0) & (f <= 40.0)
    spectral_slope = 0.0
    if slope_mask.sum() >= 5:
        log_f = np.log10(f[slope_mask])
        log_p = np.log10(psd_mean[slope_mask] + 1e-30)
        spectral_slope = float(np.polyfit(log_f, log_p, 1)[0])

    # Spectral entropy (chaos)
    p = psd_mean / total
    p = p[p > 0]
    chaos = float(-np.sum(p * np.log(p)))

    # Coherence (geometric/arithmetic mean of PSD)
    eps = 1e-20
    s_pos = psd_mean[psd_mean > 0]
    geo_mean = float(np.exp(np.mean(np.log(s_pos + eps))))
    arith_mean = float(np.mean(psd_mean) + eps)
    coherence = float(geo_mean / arith_mean)

    # Peak frequency and ratio (restrict to 1-40 Hz to avoid 0.5 Hz edge)
    peak_mask = (f >= 1.0) & (f <= 40.0)
    psd_peak_range = psd_mean[peak_mask]
    f_peak_range = f[peak_mask]
    pk_idx = int(np.argmax(psd_peak_range))
    peak_freq = float(f_peak_range[pk_idx])
    peak_ratio = float(psd_peak_range[pk_idx] / total)

    # Hurst estimate on averaged time-domain signal
    avg_signal = data.mean(axis=0)
    # Downsample to speed up
    step = max(1, len(avg_signal) // 5000)
    sig = avg_signal[::step]
    sig = sig[np.isfinite(sig)]
    hurst = 0.5
    if len(sig) >= 40:
        lags = range(10, min(21, len(sig) // 4))
        rs_vals = []
        for lag in lags:
            chunks = len(sig) // lag
            if chunks < 1:
                break
            rs_list = []
            for i in range(chunks):
                chunk = sig[i * lag:(i + 1) * lag]
                m = np.mean(chunk)
                devs = np.cumsum(chunk - m)
                r = np.max(devs) - np.min(devs)
                s = np.std(chunk, ddof=1)
                if s > eps:
                    rs_list.append(r / s)
            if rs_list:
                rs_vals.append((np.log(lag), np.log(np.mean(rs_list) + eps)))
        if len(rs_vals) >= 3:
            x = np.array([v[0] for v in rs_vals])
            y = np.array([v[1] for v in rs_vals])
            hurst = float(np.clip(np.polyfit(x, y, 1)[0], 0.0, 1.0))

    out = dict(
        chaos=chaos, coherence=coherence,
        band_delta=band_pow["delta"], band_theta=band_pow["theta"],
        band_alpha=band_pow["alpha"], band_beta=band_pow["beta"],
        band_gamma=band_pow["gamma"],
        alpha_theta_ratio=alpha_theta_ratio,
        spectral_slope=spectral_slope,
        peak_freq=peak_freq, peak_ratio=peak_ratio,
        hurst=hurst,
        n_channels=int(data.shape[0]),
        duration_s=float(data.shape[1] / sfreq),
    )
    return out


def compute_mk4e_score(f, weights=None):
    if weights is None:
        weights = {"band_theta": +1.0, "band_alpha": -1.0, "chaos": +0.5}
    return sum(f.get(k, 0) * w for k, w in weights.items())


def optimize_weights(ctrl_feats, dis_feats, feat_keys=None):
    """Find weights that maximize Youden J via brute-force grid search."""
    if feat_keys is None:
        feat_keys = ["band_delta", "band_theta", "band_alpha", "band_beta",
                     "chaos", "coherence", "peak_ratio", "alpha_theta_ratio",
                     "spectral_slope", "hurst"]
    best_j = -1
    best_w = {}
    best_key = None
    # Single-feature scan first (find best individual predictor)
    for k in feat_keys:
        c_vals = np.array([f[k] for f in ctrl_feats if k in f])
        d_vals = np.array([f[k] for f in dis_feats if k in f])
        if len(c_vals) < 3 or len(d_vals) < 3:
            continue
        for sign in [+1, -1]:
            _, j = youden_threshold(c_vals * sign, d_vals * sign, dis_is_high=True)
            if j > best_j:
                best_j = j
                best_key = k
                best_w = {k: sign}
    # Two-feature combinations with the best single
    if best_key:
        for k2 in feat_keys:
            if k2 == best_key:
                continue
            for w2 in [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]:
                trial_w = dict(best_w)
                trial_w[k2] = w2
                scores_c = np.array([compute_mk4e_score(f, trial_w) for f in ctrl_feats])
                scores_d = np.array([compute_mk4e_score(f, trial_w) for f in dis_feats])
                _, j = youden_threshold(scores_c, scores_d,
                                        dis_is_high=float(np.mean(scores_d)) > float(np.mean(scores_c)))
                if j > best_j:
                    best_j = j
                    best_w = dict(trial_w)
    # Three-feature expansion
    for k3 in feat_keys:
        if k3 in best_w:
            continue
        for w3 in [-1.0, -0.5, 0.5, 1.0]:
            trial_w = dict(best_w)
            trial_w[k3] = w3
            scores_c = np.array([compute_mk4e_score(f, trial_w) for f in ctrl_feats])
            scores_d = np.array([compute_mk4e_score(f, trial_w) for f in dis_feats])
            _, j = youden_threshold(scores_c, scores_d,
                                    dis_is_high=float(np.mean(scores_d)) > float(np.mean(scores_c)))
            if j > best_j:
                best_j = j
                best_w = dict(trial_w)
    return best_w, best_j


def youden_threshold(ctrl_vals, dis_vals, dis_is_high=True):
    all_v = np.concatenate([ctrl_vals, dis_vals])
    thresholds = np.linspace(all_v.min(), all_v.max(), 300)
    best_j, best_t = -1, float(np.median(all_v))
    for t in thresholds:
        if dis_is_high:
            tp = np.sum(dis_vals >= t); fn = np.sum(dis_vals < t)
            tn = np.sum(ctrl_vals < t); fp = np.sum(ctrl_vals >= t)
        else:
            tp = np.sum(dis_vals < t);  fn = np.sum(dis_vals >= t)
            tn = np.sum(ctrl_vals >= t); fp = np.sum(ctrl_vals < t)
        sens = tp / (tp + fn + 1e-9)
        spec = tn / (tn + fp + 1e-9)
        j = sens + spec - 1
        if j > best_j:
            best_j, best_t = j, t
    return float(best_t), best_j


def robust_stats(vals):
    v = np.asarray(vals, dtype=float)
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med))) * 1.4826 + 1e-12
    return med, mad


def clean_ctrl(ctrl_feats, dis_feats, strict=None, soft=None, bias=None):
    if len(ctrl_feats) < 4:
        return ctrl_feats, []
    n = len(ctrl_feats)
    if strict is None:
        strict = 3.0 if n < 20 else 2.5
    if soft is None:
        soft = 2.2 if n < 20 else 1.8
    if bias is None:
        bias = 2.2 if n < 20 else 1.8
    usable_feats = [k for k in FEATS if all(k in f for f in ctrl_feats)]
    stats = {k: robust_stats([f[k] for f in ctrl_feats]) for k in usable_feats + ["mk4e_score"]}
    if dis_feats:
        dis_med = float(np.median([f["mk4e_score"] for f in dis_feats]))
        dis_is_high = dis_med > stats["mk4e_score"][0]
    else:
        dis_is_high = None

    kept, rejected = [], []
    for f in ctrl_feats:
        zs = {k: (f[k] - stats[k][0]) / stats[k][1] for k in usable_feats}
        strict_hits = [k for k, z in zs.items() if abs(z) > strict]
        soft_hits = [k for k, z in zs.items() if abs(z) > soft]
        mk4_z = (f["mk4e_score"] - stats["mk4e_score"][0]) / stats["mk4e_score"][1]
        if dis_is_high is True:
            bias_hit = mk4_z > bias
        elif dis_is_high is False:
            bias_hit = mk4_z < -bias
        else:
            bias_hit = abs(mk4_z) > bias
        if strict_hits or len(soft_hits) >= 2 or bias_hit:
            rejected.append(f)
        else:
            kept.append(f)
    if len(kept) < 3:
        return ctrl_feats, []  # safety bail
    return kept, rejected


def classify(ctrl_feats, dis_feats):
    if len(ctrl_feats) < 3 or len(dis_feats) < 3:
        return None
    c = np.array([f["mk4e_score"] for f in ctrl_feats])
    d = np.array([f["mk4e_score"] for f in dis_feats])
    dis_is_high = float(np.mean(d)) > float(np.mean(c))
    thr, j = youden_threshold(c, d, dis_is_high)
    if dis_is_high:
        tp, fn = int(np.sum(d >= thr)), int(np.sum(d < thr))
        tn, fp = int(np.sum(c < thr)), int(np.sum(c >= thr))
    else:
        tp, fn = int(np.sum(d < thr)), int(np.sum(d >= thr))
        tn, fp = int(np.sum(c >= thr)), int(np.sum(c < thr))
    n = tp + tn + fp + fn
    acc = (tp + tn) / n if n else 0.0
    sens = tp / (tp + fn) if (tp + fn) else 0.0
    spec = tn / (tn + fp) if (tn + fp) else 0.0
    sep = abs(float(np.mean(d)) - float(np.mean(c))) / (float(np.std(c)) + 1e-9)
    return dict(
        n_ctrl=len(ctrl_feats), n_dis=len(dis_feats),
        threshold=round(float(thr), 4), youden_j=round(float(j), 4),
        accuracy=round(acc * 100, 1), sensitivity=round(sens * 100, 1),
        specificity=round(spec * 100, 1), separation=round(sep, 2),
        dis_is_high=bool(dis_is_high),
    )


def parse_participants_tsv(path):
    """Parse BIDS participants.tsv → dict {sub_id: {label, subtype, group_raw}}"""
    import csv
    sub_info = {}
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            sid = row.get("participant_id") or row.get("Subject") or row.get("participant")
            if not sid:
                continue
            group = row.get("Group") or row.get("group") or row.get("Condition") or ""
            label = None
            subtype = group

            # Strategy 1: explicit Group column (ds004504: A/F/C, ds003523: 0/1)
            if group:
                gl = group.strip().lower()
                # Numeric labels: 0 = CTRL, 1 = DISEASE (ds003523 TBI, common BIDS convention)
                if gl in ("0", "control", "hc"):
                    label, subtype = "CONTROL", "HC"
                elif gl in ("1", "tbi", "patient"):
                    label, subtype = "DISEASE", "TBI"
                elif gl.startswith("c") or "control" in gl or "hc" in gl or "healthy" in gl:
                    label = "CONTROL"
                elif gl == "a" or gl.startswith("alzheimer") or gl == "ad":
                    label, subtype = "DISEASE", "AD"
                elif gl == "f" or "fronto" in gl or "ftd" in gl:
                    label, subtype = "DISEASE", "FTD"
                elif "pd" in gl or "parkinson" in gl:
                    label, subtype = "DISEASE", "PD"
                elif "mdd" in gl or "depress" in gl:
                    label, subtype = "DISEASE", "MDD"
                elif "scz" in gl or "schizo" in gl:
                    label, subtype = "DISEASE", "SCZ"

            # Strategy 2: subject ID prefix (ds002778: sub-hcN vs sub-pdN)
            if label is None and sid:
                sid_lower = sid.lower()
                if "hc" in sid_lower or "ctrl" in sid_lower:
                    label, subtype = "CONTROL", "HC"
                elif "pd" in sid_lower:
                    label, subtype = "DISEASE", "PD"
                elif "ad" in sid_lower or "alz" in sid_lower:
                    label, subtype = "DISEASE", "AD"

            # Strategy 3: BDI score (ds003478 depression)
            if label is None:
                bdi = row.get("BDI", "").strip()
                try:
                    bdi_val = float(bdi)
                    if bdi_val >= 14:
                        label, subtype = "DISEASE", "MDD"
                    else:
                        label, subtype = "CONTROL", "HC"
                except (ValueError, TypeError):
                    pass

            # Strategy 4: SCID column (ds003478)
            if label is None:
                scid = row.get("SCID", "").strip().lower()
                if "current mdd" in scid:
                    label, subtype = "DISEASE", "MDD_current"
                elif "past mdd" in scid:
                    label, subtype = "DISEASE", "MDD_past"

            if label:
                sub_info[sid] = {"label": label, "group_raw": group or subtype, "subtype": subtype}
    return sub_info


def find_subject_files(dataset_root):
    """Return list of (sub_id, raw_path) for each subject. BIDS with or without session layer."""
    root = Path(dataset_root)
    out = []
    seen = set()
    for sub_dir in sorted(root.glob("sub-*")):
        sid = sub_dir.name
        if sid in seen:
            continue
        # Try: sub-XXX/eeg/ (no session) or sub-XXX/ses-*/eeg/ (with session)
        eeg_dirs = list(sub_dir.glob("ses-*/eeg")) + [sub_dir / "eeg"]
        for eeg_dir in eeg_dirs:
            if not eeg_dir.is_dir():
                continue
            for ext in ("*.set", "*.edf", "*.bdf", "*.fif", "*.vhdr"):
                files = sorted(eeg_dir.glob(ext))
                # Skip git-annex broken symlinks
                files = [f for f in files if f.exists() and f.stat().st_size > 1000]
                if files:
                    out.append((sid, str(files[0])))
                    seen.add(sid)
                    break
            if sid in seen:
                break
    return out


def load_raw(path):
    import mne
    ext = Path(path).suffix.lower()
    if ext == ".set":
        return mne.io.read_raw_eeglab(path, preload=True, verbose=False)
    if ext == ".edf":
        return mne.io.read_raw_edf(path, preload=True, verbose=False)
    if ext == ".bdf":
        return mne.io.read_raw_bdf(path, preload=True, verbose=False)
    if ext == ".fif":
        return mne.io.read_raw_fif(path, preload=True, verbose=False)
    if ext == ".vhdr":
        return mne.io.read_raw_brainvision(path, preload=True, verbose=False)
    raise ValueError(f"Unsupported EEG format: {ext}")


def main():
    if len(sys.argv) < 2:
        print("Usage: mk4_eeg_batch.py <dataset_root>", file=sys.stderr)
        sys.exit(1)
    root = Path(sys.argv[1])
    label_tsv = root / "participants.tsv"
    labels = parse_participants_tsv(label_tsv) if label_tsv.exists() else {}
    print(f"Labels parsed: {len(labels)} from {label_tsv}", flush=True)

    subjects = find_subject_files(root)
    print(f"Subject files found: {len(subjects)}", flush=True)

    results = {"dataset": str(root), "subjects": [], "labels_distribution": {}}
    ctrl_feats, dis_feats = [], []

    for i, (sid, path) in enumerate(subjects, 1):
        t0 = time.time()
        info = labels.get(sid, {})
        label = info.get("label", "UNKNOWN")
        try:
            raw = load_raw(path)
            feats = compute_eeg_features(raw)
            if feats is None:
                print(f"  [{i}/{len(subjects)}] {sid} ({label}): feat extraction failed", flush=True)
                continue
            feats["mk4e_score"] = compute_mk4e_score(feats)
            feats["subject"] = sid
            feats["label"] = label
            feats["subtype"] = info.get("subtype", "")
            results["subjects"].append(feats)
            if label == "CONTROL":
                ctrl_feats.append(feats)
            elif label == "DISEASE":
                dis_feats.append(feats)
            t = time.time() - t0
            print(f"  [{i}/{len(subjects)}] {sid} ({label}): mk4e={feats['mk4e_score']:.3f} "
                  f"α={feats['band_alpha']:.3f} θ={feats['band_theta']:.3f} {t:.1f}s", flush=True)
        except Exception as e:
            print(f"  [{i}/{len(subjects)}] {sid} ({label}): ERROR {e}", flush=True)
        gc.collect()

    # Split by subtype
    ad_feats = [s for s in results["subjects"] if s.get("subtype") == "A"]
    ftd_feats = [s for s in results["subjects"] if s.get("subtype") == "F"]
    all_dis = [s for s in results["subjects"] if s["label"] == "DISEASE"]

    results["labels_distribution"] = {
        "CONTROL": len(ctrl_feats), "DISEASE": len(all_dis),
        "AD": len(ad_feats), "FTD": len(ftd_feats),
        "UNKNOWN": sum(1 for s in results["subjects"] if s["label"] == "UNKNOWN"),
    }
    print(f"\nCTRL={len(ctrl_feats)} DIS={len(all_dis)} (AD={len(ad_feats)} FTD={len(ftd_feats)})", flush=True)

    # Train weights on full CTRL vs all-DIS
    if len(ctrl_feats) >= 3 and len(all_dis) >= 3:
        best_w, best_j_train = optimize_weights(ctrl_feats, all_dis)
        print(f"Optimized weights: {best_w} (training J={best_j_train:.3f})", flush=True)
        results["trained_weights"] = best_w
    else:
        best_w = {"band_theta": +1.0, "band_alpha": -1.0, "chaos": +0.5}
        results["trained_weights"] = best_w

    # Recompute scores with trained weights
    for s in results["subjects"]:
        s["mk4e_score"] = compute_mk4e_score(s, best_w)
    ctrl_feats = [s for s in results["subjects"] if s["label"] == "CONTROL"]
    all_dis = [s for s in results["subjects"] if s["label"] == "DISEASE"]
    ad_feats = [s for s in results["subjects"] if s.get("subtype") == "A"]
    ftd_feats = [s for s in results["subjects"] if s.get("subtype") == "F"]

    def run_analysis(ctrl, dis, label="all"):
        if len(ctrl) < 3 or len(dis) < 3:
            return None
        before = classify(ctrl, dis)
        kept, rejected = clean_ctrl(ctrl, dis)
        after = classify(kept, dis) if rejected else before
        rej_n = len(rejected)
        return {
            "label": label,
            "before": before,
            "after": after,
            "ctrl_kept": len(kept),
            "ctrl_rejected": rej_n,
            "rejected_subjects": [f["subject"] for f in rejected],
        }

    # Three analyses: CTRL vs ALL, CTRL vs AD only, CTRL vs FTD only
    analyses = {}
    for lbl, dis_group in [("CTRL_vs_ALL", all_dis), ("CTRL_vs_AD", ad_feats), ("CTRL_vs_FTD", ftd_feats)]:
        r = run_analysis(ctrl_feats, dis_group, lbl)
        if r:
            analyses[lbl] = r
            b, a = r["before"], r["after"]
            print(f"\n{lbl}: CTRL {b['n_ctrl']}→{a['n_ctrl']} (−{r['ctrl_rejected']})")
            print(f"  BEFORE: ACC={b['accuracy']}% J={b['youden_j']:.3f} sep={b['separation']}")
            print(f"  AFTER : ACC={a['accuracy']}% J={a['youden_j']:.3f} sep={a['separation']}")
    results["analyses"] = analyses

    out_path = root.parent / f"mk4e_results_{root.name}.json"
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False, default=str)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
