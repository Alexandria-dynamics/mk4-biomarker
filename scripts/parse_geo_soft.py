#!/usr/bin/env python3
"""
MK4 Biomarker - GEO Parser
Alexandria Dynamics

Podporuje:
1) GEO SOFT (.soft / .soft.gz)
2) GEO Series Matrix (.txt / .txt.gz)

Detekce formátu je automatická.
"""

import gzip
import json
import re
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.stats import entropy, ttest_ind
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ==============================================================
# PATHS
# ==============================================================

PROJECT_ROOT = Path(__file__).parent.parent
INPUT_DIR = PROJECT_ROOT / "data" / "input"
RESULTS_DIR = PROJECT_ROOT / "results"

INPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ==============================================================
# FILE TYPE DETECTION
# ==============================================================

def open_file(filepath):
    """Open file, handling .gz compression."""
    filepath = Path(filepath)
    if filepath.suffix == ".gz" or str(filepath).endswith(".gz"):
        return gzip.open(filepath, "rt", errors="ignore")
    return open(filepath, "r", errors="ignore")


def detect_file_type(filepath):
    """
    Detect if file is SOFT or MATRIX format.
    Returns: "soft", "matrix", or None
    """
    filepath = Path(filepath)

    try:
        with open_file(filepath) as f:
            # Read first 100 lines for detection
            lines = []
            for i, line in enumerate(f):
                if i >= 100:
                    break
                lines.append(line)

        content = "".join(lines)

        # SOFT format markers
        if "^SAMPLE" in content or "^SERIES" in content or "^PLATFORM" in content:
            return "soft"

        # Series Matrix format markers
        if "!Series_title" in content or "!series_matrix_table_begin" in content:
            return "matrix"

        # Try to detect tabular expression data
        # Look for header line with sample IDs and data rows with gene IDs
        for i, line in enumerate(lines):
            if line.startswith("!") or line.startswith("#") or line.startswith("^"):
                continue

            # First non-comment line - check if it's a header
            parts = line.strip().split("\t")
            if len(parts) > 2:
                # Check next few lines for gene IDs
                for j in range(i + 1, min(i + 5, len(lines))):
                    next_parts = lines[j].strip().split("\t")
                    if len(next_parts) > 2:
                        # Looks like gene expression matrix
                        first_col = next_parts[0].upper()
                        if (first_col.startswith("ENSG") or
                            first_col.startswith("ILMN") or
                            first_col.isdigit() or
                            "_at" in first_col):
                            return "matrix"
            break

        return None

    except Exception as e:
        print(f"Error detecting file type: {e}")
        return None

# ==============================================================
# LABEL INFERENCE
# ==============================================================

CONTROL_KEYWORDS = [
    "control", "ctrl", "normal", "healthy", "non-diabetic", "nondiabetic",
    "wild type", "wildtype", "wt", "sham", "untreated", "vehicle",
    "baseline", "reference", "neg", "benign", "nontumor", "non-tumor",
    "adjacent normal", "nad", "matched normal"
]

DISEASE_KEYWORDS = [
    "disease", "patient", "diabetic", "diabetes", "t2d", "t1d",
    "tumor", "tumour", "cancer", "carcinoma", "adenoma", "malignant",
    "multiple sclerosis", "ms", "parkinson", "pd", "alzheimer", "ad",
    "treated", "case", "affected", "positive", "metastatic",
    "stage", "grade", "relapse", "progressive"
]


def infer_label_from_text(text: str):
    """Infer CONTROL/DISEASE label from metadata text."""
    if not text:
        return None

    t = text.lower()

    # Check control keywords first (more specific)
    for kw in CONTROL_KEYWORDS:
        if kw in t:
            return "CONTROL"

    # Then check disease keywords
    for kw in DISEASE_KEYWORDS:
        if kw in t:
            return "DISEASE"

    return None

# ==============================================================
# SOFT PARSER
# ==============================================================

def parse_soft_file(filepath):
    """
    Parse GEO SOFT file format.
    Handles both full SOFT and family SOFT formats.
    """
    filepath = Path(filepath)

    # Extract dataset name
    name = filepath.name
    for suffix in [".soft.gz", ".soft", "_family.soft.gz", "_family.soft", ".gz"]:
        name = name.replace(suffix, "")
    dataset_name = name

    print("=" * 70)
    print(f"Parsing SOFT: {filepath.name}")
    print("=" * 70)

    samples = {}
    current_sample = None
    in_table = False
    current_data = {}
    meta_buffer = []

    try:
        with open_file(filepath) as f:
            for line in f:
                line = line.rstrip("\n\r")

                # New sample block
                if line.startswith("^SAMPLE"):
                    # Save previous sample if exists
                    if current_sample and current_data:
                        samples[current_sample]["expression"] = current_data
                        if samples[current_sample]["label"] is None:
                            lbl = infer_label_from_text(" ".join(meta_buffer))
                            samples[current_sample]["label"] = lbl

                    # Start new sample
                    parts = line.split("=", 1)
                    current_sample = parts[1].strip() if len(parts) > 1 else f"sample_{len(samples)}"
                    samples[current_sample] = {"label": None, "expression": {}}
                    current_data = {}
                    meta_buffer = []
                    in_table = False

                # Collect metadata for label inference
                elif current_sample and line.startswith("!Sample_"):
                    meta_buffer.append(line)

                    # Try to extract label from specific fields
                    line_lower = line.lower()
                    if any(x in line_lower for x in ["characteristic", "source", "title", "description"]):
                        lbl = infer_label_from_text(line)
                        if lbl and samples[current_sample]["label"] is None:
                            samples[current_sample]["label"] = lbl

                # Table markers
                elif line == "!sample_table_begin":
                    in_table = True
                    current_data = {}

                elif line == "!sample_table_end":
                    in_table = False
                    if current_sample and current_data:
                        samples[current_sample]["expression"] = current_data
                        if samples[current_sample]["label"] is None:
                            lbl = infer_label_from_text(" ".join(meta_buffer))
                            samples[current_sample]["label"] = lbl

                # Expression data rows
                elif in_table and "\t" in line:
                    if line.startswith("ID_REF") or line.startswith("IDENTIFIER"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        try:
                            gene_id = parts[0].strip()
                            value = float(parts[1])
                            if gene_id and not np.isnan(value):
                                current_data[gene_id] = value
                        except (ValueError, IndexError):
                            pass

        # Handle last sample
        if current_sample and current_data and current_sample in samples:
            samples[current_sample]["expression"] = current_data
            if samples[current_sample]["label"] is None:
                lbl = infer_label_from_text(" ".join(meta_buffer))
                samples[current_sample]["label"] = lbl

    except Exception as e:
        print(f"Error parsing SOFT: {e}")
        raise

    # Filter to valid samples only
    valid_samples = {
        k: v for k, v in samples.items()
        if v.get("label") is not None and len(v.get("expression", {})) > 0
    }

    n_ctrl = sum(1 for v in valid_samples.values() if v["label"] == "CONTROL")
    n_dis = sum(1 for v in valid_samples.values() if v["label"] == "DISEASE")

    print(f"Found {len(valid_samples)} valid samples (Control: {n_ctrl}, Disease: {n_dis})")

    if len(valid_samples) == 0:
        raise ValueError(f"No valid samples found in {filepath.name}")

    return {
        "dataset_name": dataset_name,
        "samples": valid_samples,
        "source_file": filepath.name
    }

# ==============================================================
# MATRIX PARSER (Series Matrix)
# ==============================================================

def parse_matrix_file(filepath):
    """
    Parse GEO Series Matrix file format.
    These have metadata at top, then expression matrix.
    """
    filepath = Path(filepath)

    # Extract dataset name
    name = filepath.name
    for suffix in [".txt.gz", ".txt", "_series_matrix.txt.gz", "_series_matrix.txt", ".gz"]:
        name = name.replace(suffix, "")
    dataset_name = name

    print("=" * 70)
    print(f"Parsing MATRIX: {filepath.name}")
    print("=" * 70)

    # Read metadata section
    sample_titles = {}
    sample_chars = {}
    data_start_line = 0

    try:
        with open_file(filepath) as f:
            for i, line in enumerate(f):
                line = line.rstrip("\n\r")

                if line.startswith("!Sample_title"):
                    parts = line.split("\t")
                    for j, title in enumerate(parts[1:], 1):
                        sample_titles[j] = title.strip('"')

                elif line.startswith("!Sample_characteristics") or line.startswith("!Sample_source"):
                    parts = line.split("\t")
                    for j, char in enumerate(parts[1:], 1):
                        if j not in sample_chars:
                            sample_chars[j] = []
                        sample_chars[j].append(char.strip('"'))

                elif line.startswith("!series_matrix_table_begin"):
                    data_start_line = i + 1
                    break

                elif line.startswith("ID_REF") or (not line.startswith("!") and "\t" in line and not line.startswith("#")):
                    # Found header without explicit marker
                    data_start_line = i
                    break

        # Read expression data
        with open_file(filepath) as f:
            # Skip to data section
            for _ in range(data_start_line):
                next(f)

            # Read as DataFrame
            df = pd.read_csv(f, sep="\t", index_col=0, na_values=["null", "NA", "NaN", ""])

        # Remove any trailing marker rows
        df = df[~df.index.astype(str).str.startswith("!")]

        # Convert to samples dict
        samples = {}
        for col_idx, col_name in enumerate(df.columns, 1):
            # Get sample ID
            sample_id = str(col_name).strip('"')

            # Infer label from metadata
            label = None
            meta_text = sample_titles.get(col_idx, "")
            if col_idx in sample_chars:
                meta_text += " " + " ".join(sample_chars.get(col_idx, []))

            label = infer_label_from_text(meta_text)

            # If still no label, try column name
            if label is None:
                label = infer_label_from_text(col_name)

            # Get expression values
            expr_series = pd.to_numeric(df[col_name], errors="coerce").dropna()
            expression = dict(zip(df.index[expr_series.index], expr_series.values))

            if len(expression) > 0:
                samples[sample_id] = {
                    "label": label,
                    "expression": expression
                }

    except Exception as e:
        print(f"Error parsing MATRIX: {e}")
        raise

    # Filter to valid samples
    valid_samples = {
        k: v for k, v in samples.items()
        if v.get("label") is not None and len(v.get("expression", {})) > 0
    }

    n_ctrl = sum(1 for v in valid_samples.values() if v["label"] == "CONTROL")
    n_dis = sum(1 for v in valid_samples.values() if v["label"] == "DISEASE")

    print(f"Found {len(valid_samples)} valid samples (Control: {n_ctrl}, Disease: {n_dis})")

    if len(valid_samples) == 0:
        raise ValueError(f"No valid samples found in {filepath.name}")

    return {
        "dataset_name": dataset_name,
        "samples": valid_samples,
        "source_file": filepath.name
    }

# ==============================================================
# UNIFIED PARSER
# ==============================================================

def parse_file(filepath):
    """
    Unified entry point - detects format and parses accordingly.
    Returns dataset dict or raises exception.
    """
    filepath = Path(filepath)

    ftype = detect_file_type(filepath)

    if ftype == "soft":
        return parse_soft_file(filepath)
    elif ftype == "matrix":
        return parse_matrix_file(filepath)
    else:
        raise ValueError(f"Unknown or unsupported file format: {filepath.name}")

# ==============================================================
# CHAOS ANALYSIS
# ==============================================================

def compute_sample_entropy(expression_dict):
    """Compute Shannon entropy for one sample's expression profile."""
    values = np.array(list(expression_dict.values()), dtype=float)
    values = values[~np.isnan(values)]
    values = values[values > 0]

    if len(values) == 0:
        return np.nan

    # Normalize to probability distribution
    values = values / np.sum(values)
    return float(entropy(values, base=2))


def compute_chaos_metrics(dataset):
    """
    Compute chaos metrics for all samples in dataset.
    Returns dict with entropies and statistics.
    """
    entropies = {"CONTROL": [], "DISEASE": []}
    sample_details = {}

    for sample_id, sample_data in dataset["samples"].items():
        label = sample_data.get("label")
        if label not in ("CONTROL", "DISEASE"):
            continue

        ent = compute_sample_entropy(sample_data["expression"])

        sample_details[sample_id] = {
            "label": label,
            "entropy": ent,
            "n_genes": len(sample_data["expression"])
        }

        if not np.isnan(ent):
            entropies[label].append(ent)

    # Calculate statistics
    c = entropies["CONTROL"]
    d = entropies["DISEASE"]

    stats = {
        "n_control": len(c),
        "n_disease": len(d)
    }

    if len(c) >= 2 and len(d) >= 2:
        t_stat, p_value = ttest_ind(c, d, equal_var=False)
        stats.update({
            "control_mean": float(np.mean(c)),
            "control_std": float(np.std(c, ddof=1)),
            "disease_mean": float(np.mean(d)),
            "disease_std": float(np.std(d, ddof=1)),
            "difference": float(np.mean(d) - np.mean(c)),
            "t_statistic": float(t_stat),
            "p_value": float(p_value),
            "significant": bool(p_value < 0.05)
        })
    else:
        stats["error"] = "Insufficient samples for statistical test"

    return {
        "entropies": entropies,
        "statistics": stats,
        "sample_details": sample_details
    }

# ==============================================================
# OUTPUT FUNCTIONS
# ==============================================================

def create_output_directory(base_name):
    """Create timestamped output directory."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = RESULTS_DIR / f"{timestamp}_{base_name}"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def create_dataset_subdirectory(run_dir, dataset_name):
    """Create subdirectory for individual dataset within run directory."""
    dataset_dir = run_dir / dataset_name
    dataset_dir.mkdir(parents=True, exist_ok=True)
    (dataset_dir / "figures").mkdir(exist_ok=True)
    return dataset_dir


def save_results(output_dir, dataset, chaos_results):
    """Save analysis results as JSON."""
    output_dir = Path(output_dir)

    results = {
        "dataset_name": dataset.get("dataset_name"),
        "source_file": dataset.get("source_file"),
        "generated_at": datetime.now().isoformat(),
        "chaos_results": chaos_results
    }

    (output_dir / "results.json").write_text(
        json.dumps(results, indent=2, default=str),
        encoding="utf-8"
    )


def generate_plot(output_dir, chaos_results, dataset_name=""):
    """Generate boxplot visualization."""
    output_dir = Path(output_dir)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    c = chaos_results["entropies"]["CONTROL"]
    d = chaos_results["entropies"]["DISEASE"]
    stats = chaos_results["statistics"]

    fig, ax = plt.subplots(figsize=(8, 6))

    bp = ax.boxplot([c, d], labels=["Control", "Disease"], patch_artist=True)
    bp["boxes"][0].set_facecolor("#2ecc71")
    bp["boxes"][1].set_facecolor("#e74c3c")
    for box in bp["boxes"]:
        box.set_alpha(0.7)

    ax.set_ylabel("Shannon Entropy", fontsize=12)
    ax.set_title(f"Transcriptional Chaos: {dataset_name}", fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)

    # Add p-value annotation
    if "p_value" in stats:
        p = stats["p_value"]
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.text(0.5, 0.95, f"p = {p:.4f} {sig}", transform=ax.transAxes, ha="center", fontsize=11)

    plt.tight_layout()
    plt.savefig(fig_dir / "chaos_boxplot.png", dpi=150, bbox_inches="tight")
    plt.close()


def generate_html_report(output_dir, dataset, chaos_results):
    """Generate HTML report for single dataset."""
    import base64
    from io import BytesIO

    output_dir = Path(output_dir)
    stats = chaos_results["statistics"]

    # Generate plot as base64
    c = chaos_results["entropies"]["CONTROL"]
    d = chaos_results["entropies"]["DISEASE"]

    fig, ax = plt.subplots(figsize=(10, 6))
    if c and d:
        bp = ax.boxplot([c, d], labels=["Control", "Disease"], patch_artist=True)
        bp["boxes"][0].set_facecolor("#2ecc71")
        bp["boxes"][1].set_facecolor("#e74c3c")
        for box in bp["boxes"]:
            box.set_alpha(0.7)
    ax.set_ylabel("Shannon Entropy")
    ax.set_title(f"Chaos Analysis: {dataset['dataset_name']}")
    ax.grid(axis="y", alpha=0.3)

    if "p_value" in stats:
        sig = "***" if stats["p_value"] < 0.001 else "**" if stats["p_value"] < 0.01 else "*" if stats["p_value"] < 0.05 else "ns"
        ax.text(0.5, 0.95, f"p = {stats['p_value']:.4f} {sig}", transform=ax.transAxes, ha="center")

    buf = BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close()
    buf.seek(0)
    plot_b64 = base64.b64encode(buf.read()).decode("utf-8")

    # Build HTML
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>MK4 Report - {dataset['dataset_name']}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #f5f7fa; }}
        .container {{ max-width: 1000px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #667eea; }}
        .stats {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin: 20px 0; }}
        .stat-card {{ padding: 20px; border-radius: 8px; background: #f8f9fa; border-left: 4px solid #667eea; }}
        .stat-value {{ font-size: 24px; font-weight: bold; color: #333; }}
        .stat-label {{ color: #666; font-size: 14px; }}
        img {{ max-width: 100%; margin: 20px 0; }}
        pre {{ background: #272822; color: #f8f8f2; padding: 15px; border-radius: 8px; overflow-x: auto; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>MK4 Biomarker Analysis</h1>
        <p>Dataset: <strong>{dataset['dataset_name']}</strong> | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

        <div class="stats">
            <div class="stat-card">
                <div class="stat-value">{stats.get('n_control', 0)}</div>
                <div class="stat-label">Control Samples</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{stats.get('n_disease', 0)}</div>
                <div class="stat-label">Disease Samples</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{stats.get('p_value', 'N/A'):.4f if isinstance(stats.get('p_value'), float) else 'N/A'}</div>
                <div class="stat-label">p-value</div>
            </div>
        </div>

        <h2>Visualization</h2>
        <img src="data:image/png;base64,{plot_b64}" alt="Chaos Distribution">

        <h2>Raw Data</h2>
        <details>
            <summary>View JSON</summary>
            <pre>{json.dumps(chaos_results, indent=2, default=str)}</pre>
        </details>
    </div>
</body>
</html>"""

    (output_dir / "index.html").write_text(html, encoding="utf-8")


# ==============================================================
# MAIN
# ==============================================================

def run(filepath):
    """Run analysis on single file."""
    dataset = parse_file(filepath)
    chaos_results = compute_chaos_metrics(dataset)

    output_dir = create_output_directory(dataset["dataset_name"])
    save_results(output_dir, dataset, chaos_results)
    generate_plot(output_dir, chaos_results, dataset["dataset_name"])
    generate_html_report(output_dir, dataset, chaos_results)

    print(f"\nResults saved to: {output_dir}")
    return dataset, chaos_results, output_dir


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python parse_geo_soft.py <file>")
        sys.exit(1)

    run(sys.argv[1])
