#!/usr/bin/env python3
"""
MK4 Biomarker - Universal GEO Auto-Parser
Alexandria Dynamics

Zero-configuration parser that automatically:
- Detects file format (SOFT/Matrix)
- Finds label fields (any format)
- Matches samples across files
- Computes chaos metrics
- Provides interpretation

Usage:
    python auto_parse_geo.py <file.gz>
    python auto_parse_geo.py <file1.gz> <file2.gz>  # SOFT + Matrix pair
"""

import gzip
import json
import re
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import numpy as np
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
# KEYWORD LISTS
# ==============================================================

CONTROL_KEYWORDS = [
    # Direct terms
    "control", "ctrl", "normal", "healthy", "hc",
    # Untreated
    "untreated", "vehicle", "baseline", "reference",
    # Wild type
    "wild type", "wildtype", "wt", "wild-type",
    # Negative
    "non-diabetic", "nondiabetic", "non diabetic",
    "non-tumor", "nontumor", "non tumor",
    "non-cancer", "noncancer", "non cancer",
    "benign", "negative", "neg",
    # Adjacent/matched
    "adjacent normal", "matched normal", "nad", "nat",
    # Sham
    "sham", "mock", "placebo",
]

DISEASE_KEYWORDS = [
    # Cancer general
    "cancer", "tumor", "tumour", "carcinoma", "adenoma", "malignant",
    "neoplasm", "metastatic", "metastasis", "primary tumor",
    # Cancer types
    "breast cancer", "lung cancer", "colorectal", "crc",
    "pancreatic", "pdac", "glioblastoma", "gbm", "melanoma",
    "prostate", "ovarian", "hepatocellular", "hcc", "leukemia",
    "lymphoma", "myeloma", "sarcoma", "bladder", "kidney", "renal",
    # Diabetes
    "diabetic", "diabetes", "t2d", "t1d", "type 2 diabetes", "type 1 diabetes",
    "t2dm", "t1dm", "dm2", "dm1", "hyperglycemic",
    # Neurological
    "parkinson", "pd", "alzheimer", "ad", "multiple sclerosis", "ms",
    "huntington", "als", "dementia",
    # Other diseases
    "disease", "patient", "case", "affected", "positive",
    "treated", "resistant", "progressive", "relapse", "recurrent",
    "stage", "grade", "severe", "advanced", "chronic",
    # Inflammatory
    "inflammatory", "inflamed", "crohn", "colitis", "arthritis",
    "lupus", "psoriasis", "asthma",
]

# Fields that commonly contain labels
LABEL_FIELD_PATTERNS = [
    r"disease\s*state",
    r"disease\s*status",
    r"sample\s*type",
    r"sample\s*group",
    r"tissue\s*type",
    r"tissue\s*status",
    r"cancer\s*type",
    r"tumor\s*type",
    r"tumou?r\s*status",
    r"diagnosis",
    r"condition",
    r"phenotype",
    r"clinical\s*status",
    r"health\s*status",
    r"patient\s*status",
    r"subject\s*type",
    r"subject\s*group",
    r"histology",
    r"pathology",
    r"sample\s*source",
    r"source\s*name",
    r"title",
    r"description",
    r"characteristic",
]

# ==============================================================
# FILE HANDLING
# ==============================================================

def open_file(filepath):
    """Open file, handling .gz compression."""
    filepath = Path(filepath)
    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, "rt", errors="ignore")
    return open(filepath, "r", errors="ignore")


def read_lines(filepath, max_lines=None):
    """Read lines from file."""
    lines = []
    with open_file(filepath) as f:
        for i, line in enumerate(f):
            if max_lines and i >= max_lines:
                break
            lines.append(line.rstrip("\n\r"))
    return lines

# ==============================================================
# FORMAT DETECTION
# ==============================================================

def detect_file_format(filepath):
    """
    Auto-detect file format.
    Returns: 'SOFT', 'MATRIX', or 'UNKNOWN'
    """
    lines = read_lines(filepath, max_lines=200)
    content = "\n".join(lines)

    # SOFT markers (GEO SOFT format)
    soft_markers = ["^DATABASE", "^SERIES", "^PLATFORM", "^SAMPLE"]
    if any(marker in content for marker in soft_markers):
        return "SOFT"

    # Series Matrix markers
    if "!Series_title" in content or "!series_matrix_table_begin" in content:
        return "MATRIX"

    # Try to detect plain matrix (TSV with gene rows)
    for i, line in enumerate(lines):
        if line.startswith("!") or line.startswith("#") or line.startswith("^"):
            continue

        # First data line - check if it's tabular
        parts = line.split("\t")
        if len(parts) > 3:
            # Check subsequent lines for gene IDs
            for j in range(i + 1, min(i + 10, len(lines))):
                if j < len(lines):
                    row_parts = lines[j].split("\t")
                    if len(row_parts) > 3:
                        first_col = row_parts[0].upper()
                        # Gene ID patterns
                        if (first_col.startswith("ENSG") or
                            first_col.startswith("ILMN") or
                            "_at" in first_col or
                            first_col.isdigit() or
                            re.match(r"^[A-Z0-9]+$", first_col)):
                            return "MATRIX"
        break

    return "UNKNOWN"


def has_expression_tables(filepath):
    """Check if SOFT file has embedded expression tables."""
    lines = read_lines(filepath, max_lines=5000)
    content = "\n".join(lines)
    return "!sample_table_begin" in content.lower()

# ==============================================================
# LABEL DETECTION
# ==============================================================

def smart_label_detection(text):
    """
    Detect CONTROL or DISEASE from any text.
    Returns: 'CONTROL', 'DISEASE', or None
    """
    if not text:
        return None

    t = text.lower()

    # Check CONTROL keywords first (they're more specific)
    for kw in CONTROL_KEYWORDS:
        if kw in t:
            # Make sure it's not negated
            # e.g., "non-control" should not match
            idx = t.find(kw)
            prefix = t[max(0, idx-4):idx]
            if "non" not in prefix and "not" not in prefix:
                return "CONTROL"

    # Check DISEASE keywords
    for kw in DISEASE_KEYWORDS:
        if kw in t:
            return "DISEASE"

    return None


def find_label_fields(lines):
    """
    Scan lines for potential label fields.
    Returns: list of (field_name, sample_values) tuples
    """
    potential_fields = []

    for line in lines:
        line_lower = line.lower()

        # Check if this line matches any label field pattern
        for pattern in LABEL_FIELD_PATTERNS:
            if re.search(pattern, line_lower):
                # This might be a label field
                # Check if it contains control/disease keywords
                if smart_label_detection(line) is not None:
                    potential_fields.append(line)
                break

    return potential_fields

# ==============================================================
# SOFT PARSER
# ==============================================================

def parse_soft_file(filepath):
    """
    Parse SOFT file with auto-detection of label fields.
    Handles both embedded expression tables and metadata-only files.
    """
    filepath = Path(filepath)
    dataset_name = extract_dataset_name(filepath)

    print(f"  Parsing SOFT: {filepath.name}")

    # Get file size for progress
    file_size = filepath.stat().st_size
    print(f"  File size: {file_size / 1024 / 1024:.1f} MB")

    has_tables = has_expression_tables(filepath)
    print(f"  Expression tables embedded: {'Yes' if has_tables else 'No'}")

    samples = {}
    current_sample = None
    in_table = False
    current_data = {}
    meta_buffer = []

    # Track which fields contain labels
    label_field_stats = defaultdict(lambda: {"control": 0, "disease": 0, "unknown": 0})

    # Progress tracking
    lines_read = 0
    bytes_read = 0
    last_progress = 0

    with open_file(filepath) as f:
        for line in f:
            lines_read += 1
            bytes_read += len(line)

            # Progress every 10%
            progress = int((bytes_read / file_size) * 100) // 10 * 10
            if progress > last_progress:
                last_progress = progress
                print(f"    Progress: {progress}% ({lines_read:,} lines)")

            line = line.rstrip("\n\r")

            # New sample block
            if line.startswith("^SAMPLE"):
                # Save previous sample
                if current_sample and current_sample in samples:
                    if current_data:
                        samples[current_sample]["expression"] = current_data

                    # Infer label from metadata
                    if samples[current_sample]["label"] is None:
                        label, field = infer_label_from_metadata(meta_buffer)
                        samples[current_sample]["label"] = label
                        if field:
                            if label == "CONTROL":
                                label_field_stats[field]["control"] += 1
                            elif label == "DISEASE":
                                label_field_stats[field]["disease"] += 1

                # Start new sample
                parts = line.split("=", 1)
                current_sample = parts[1].strip() if len(parts) > 1 else f"sample_{len(samples)}"
                samples[current_sample] = {"label": None, "expression": {}, "metadata": {}}
                current_data = {}
                meta_buffer = []
                in_table = False

            # Collect metadata
            elif current_sample and line.startswith("!Sample_"):
                meta_buffer.append(line)

                # Extract field name and value
                if "=" in line:
                    field, value = line.split("=", 1)
                    field = field.strip()
                    value = value.strip()

                    if field not in samples[current_sample]["metadata"]:
                        samples[current_sample]["metadata"][field] = []
                    samples[current_sample]["metadata"][field].append(value)

                    # Try to detect label on the fly
                    if samples[current_sample]["label"] is None:
                        label = smart_label_detection(value)
                        if label:
                            samples[current_sample]["label"] = label
                            if label == "CONTROL":
                                label_field_stats[field]["control"] += 1
                            else:
                                label_field_stats[field]["disease"] += 1

            # Expression table
            elif line.lower() == "!sample_table_begin":
                in_table = True
                current_data = {}

            elif line.lower() == "!sample_table_end":
                in_table = False
                if current_sample and current_data:
                    samples[current_sample]["expression"] = current_data

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
    if current_sample and current_sample in samples:
        if current_data:
            samples[current_sample]["expression"] = current_data
        if samples[current_sample]["label"] is None:
            label, field = infer_label_from_metadata(meta_buffer)
            samples[current_sample]["label"] = label

    # Report label field statistics
    if label_field_stats:
        print(f"  Label fields detected:")
        for field, counts in sorted(label_field_stats.items(),
                                     key=lambda x: sum(x[1].values()), reverse=True)[:3]:
            total = sum(counts.values())
            print(f"    - {field}: {counts['control']} ctrl, {counts['disease']} disease")

    return {
        "dataset_name": dataset_name,
        "samples": samples,
        "source_file": filepath.name,
        "format": "SOFT",
        "has_expression": has_tables
    }


def infer_label_from_metadata(meta_lines):
    """
    Infer label from list of metadata lines.
    Returns: (label, field_name) or (None, None)
    """
    # Priority order: more specific fields first
    priority_patterns = [
        r"disease\s*state",
        r"disease\s*status",
        r"sample\s*type",
        r"diagnosis",
        r"condition",
        r"phenotype",
        r"tissue",
        r"source",
        r"title",
        r"description",
        r"characteristic",
    ]

    for pattern in priority_patterns:
        for line in meta_lines:
            if re.search(pattern, line.lower()):
                label = smart_label_detection(line)
                if label:
                    # Extract field name
                    if "=" in line:
                        field = line.split("=")[0].strip()
                    else:
                        field = line.split(":")[0].strip() if ":" in line else "unknown"
                    return label, field

    # Fallback: check all lines
    for line in meta_lines:
        label = smart_label_detection(line)
        if label:
            return label, "metadata"

    return None, None

# ==============================================================
# MATRIX PARSER
# ==============================================================

def parse_matrix_file(filepath, soft_metadata=None):
    """
    Parse matrix file with optional SOFT metadata for labels.
    """
    filepath = Path(filepath)
    dataset_name = extract_dataset_name(filepath)

    # Get file size for progress
    file_size = filepath.stat().st_size
    print(f"  Parsing MATRIX: {filepath.name}")
    print(f"  File size: {file_size / 1024 / 1024:.1f} MB")
    print(f"  Reading metadata...")

    # First pass: extract metadata from header
    sample_metadata = {}  # col_index -> metadata dict
    data_start_line = 0

    with open_file(filepath) as f:
        for i, line in enumerate(f):
            line = line.rstrip("\n\r")

            # Sample titles
            if line.startswith("!Sample_title"):
                parts = line.split("\t")
                for j, title in enumerate(parts[1:], 1):
                    if j not in sample_metadata:
                        sample_metadata[j] = {"texts": []}
                    sample_metadata[j]["title"] = title.strip('"')
                    sample_metadata[j]["texts"].append(title)

            # Sample characteristics (can have multiple lines)
            elif line.startswith("!Sample_characteristics") or \
                 line.startswith("!Sample_source") or \
                 line.startswith("!Sample_description"):
                parts = line.split("\t")
                for j, char in enumerate(parts[1:], 1):
                    if j not in sample_metadata:
                        sample_metadata[j] = {"texts": []}
                    sample_metadata[j]["texts"].append(char.strip('"'))

            # Sample geo accession
            elif line.startswith("!Sample_geo_accession"):
                parts = line.split("\t")
                for j, gsm in enumerate(parts[1:], 1):
                    if j not in sample_metadata:
                        sample_metadata[j] = {"texts": []}
                    sample_metadata[j]["gsm"] = gsm.strip('"')

            # Data section markers
            elif line.startswith("!series_matrix_table_begin"):
                data_start_line = i + 1
                break

            elif line.startswith("ID_REF") or line.startswith('"ID_REF"'):
                data_start_line = i
                break

    # Read expression data
    import pandas as pd

    print(f"  Loading expression matrix (this may take a while for large files)...")

    with open_file(filepath) as f:
        # Skip to data section
        for _ in range(data_start_line):
            next(f)

        df = pd.read_csv(f, sep="\t", index_col=0, na_values=["null", "NA", "NaN", ""])

    # Remove marker rows
    df = df[~df.index.astype(str).str.startswith("!")]

    print(f"  Loaded: {len(df)} genes x {len(df.columns)} samples")

    # Build samples dict
    samples = {}

    for col_idx, col_name in enumerate(df.columns, 1):
        sample_id = str(col_name).strip('"')

        # Get metadata for this column
        meta = sample_metadata.get(col_idx, {"texts": []})

        # Use GSM if available
        if "gsm" in meta:
            sample_id = meta["gsm"]

        # Infer label from metadata
        label = None
        all_text = " ".join(meta.get("texts", []))
        label = smart_label_detection(all_text)

        # If no label from metadata, try column name
        if label is None:
            label = smart_label_detection(col_name)

        # If soft metadata provided, try to match
        if label is None and soft_metadata:
            label = match_sample_to_soft(sample_id, meta, soft_metadata)

        # Get expression values
        expr_series = pd.to_numeric(df[col_name], errors="coerce").dropna()
        expression = dict(zip(expr_series.index, expr_series.values))

        if len(expression) > 0:
            samples[sample_id] = {
                "label": label,
                "expression": expression,
                "metadata": meta
            }

    n_genes = len(df.index) if len(df) > 0 else 0

    return {
        "dataset_name": dataset_name,
        "samples": samples,
        "source_file": filepath.name,
        "format": "MATRIX",
        "n_genes": n_genes
    }


def match_sample_to_soft(sample_id, matrix_meta, soft_data):
    """
    Try to match matrix sample to SOFT metadata for label.
    """
    if not soft_data or "samples" not in soft_data:
        return None

    # Direct GSM match
    if sample_id in soft_data["samples"]:
        return soft_data["samples"][sample_id].get("label")

    # Try matching by title
    title = matrix_meta.get("title", "")
    for gsm, soft_sample in soft_data["samples"].items():
        soft_meta = soft_sample.get("metadata", {})
        soft_titles = soft_meta.get("!Sample_title", [])
        if title and any(title in t for t in soft_titles):
            return soft_sample.get("label")

    return None

# ==============================================================
# UTILITY
# ==============================================================

def extract_dataset_name(filepath):
    """Extract clean dataset name from filename."""
    name = filepath.name
    for suffix in [".soft.gz", ".soft", ".txt.gz", ".txt",
                   "_family.soft.gz", "_family.soft",
                   "_series_matrix.txt.gz", "_series_matrix.txt",
                   "_data_matrix.txt.gz", "_data_matrix.txt",
                   ".gz"]:
        name = name.replace(suffix, "")
    return name

# ==============================================================
# ANALYSIS
# ==============================================================

def compute_sample_entropy(expression_dict):
    """Compute Shannon entropy for expression profile."""
    values = np.array(list(expression_dict.values()), dtype=float)
    values = values[~np.isnan(values)]
    values = values[values > 0]

    if len(values) == 0:
        return np.nan

    values = values / np.sum(values)
    return float(entropy(values, base=2))


def analyze_dataset(dataset):
    """
    Compute chaos metrics and statistics.
    """
    samples = dataset.get("samples", {})

    control_entropy = []
    disease_entropy = []
    sample_details = {}

    for sample_id, sample_data in samples.items():
        label = sample_data.get("label")
        expression = sample_data.get("expression", {})

        if not expression:
            continue

        ent = compute_sample_entropy(expression)

        sample_details[sample_id] = {
            "label": label,
            "entropy": ent,
            "n_genes": len(expression)
        }

        if label == "CONTROL" and not np.isnan(ent):
            control_entropy.append(ent)
        elif label == "DISEASE" and not np.isnan(ent):
            disease_entropy.append(ent)

    # Statistics
    results = {
        "n_total": len(samples),
        "n_control": len(control_entropy),
        "n_disease": len(disease_entropy),
        "n_unlabeled": len(samples) - len(control_entropy) - len(disease_entropy),
        "sample_details": sample_details
    }

    if len(control_entropy) >= 2 and len(disease_entropy) >= 2:
        control_mean = np.mean(control_entropy)
        disease_mean = np.mean(disease_entropy)

        t_stat, p_value = ttest_ind(control_entropy, disease_entropy, equal_var=False)

        difference = disease_mean - control_mean
        percent_change = (difference / control_mean) * 100 if control_mean != 0 else 0
        direction = "INCREASE" if difference > 0 else "DECREASE"

        results.update({
            "control_mean": float(control_mean),
            "control_std": float(np.std(control_entropy, ddof=1)),
            "disease_mean": float(disease_mean),
            "disease_std": float(np.std(disease_entropy, ddof=1)),
            "difference": float(difference),
            "percent_change": float(percent_change),
            "direction": direction,
            "t_statistic": float(t_stat),
            "p_value": float(p_value),
            "significant": bool(p_value < 0.05),
            "control_entropy": control_entropy,
            "disease_entropy": disease_entropy,
        })
    else:
        results["error"] = "Insufficient samples for statistical analysis"
        results["control_entropy"] = control_entropy
        results["disease_entropy"] = disease_entropy

    return results


def interpret_results(results):
    """
    Generate human-readable interpretation.
    """
    lines = []

    if "error" in results:
        lines.append(f"⚠️  {results['error']}")
        lines.append(f"   Control samples: {results['n_control']}")
        lines.append(f"   Disease samples: {results['n_disease']}")
        lines.append(f"   Unlabeled: {results['n_unlabeled']}")
        return "\n".join(lines)

    # Direction and significance
    direction = results.get("direction", "UNKNOWN")
    p_value = results.get("p_value", 1.0)
    percent = results.get("percent_change", 0)

    if p_value < 0.001:
        sig_stars = "***"
        sig_text = "highly significant"
    elif p_value < 0.01:
        sig_stars = "**"
        sig_text = "very significant"
    elif p_value < 0.05:
        sig_stars = "*"
        sig_text = "significant"
    else:
        sig_stars = "ns"
        sig_text = "not significant"

    lines.append(f"  Control entropy: {results['control_mean']:.6f} ± {results['control_std']:.6f}")
    lines.append(f"  Disease entropy: {results['disease_mean']:.6f} ± {results['disease_std']:.6f}")
    lines.append(f"  Direction: {direction} ({percent:+.2f}%)")
    lines.append(f"  p-value: {p_value:.6f} {sig_stars}")
    lines.append(f"  Significant: {'YES ✓' if results['significant'] else 'NO'}")

    # Interpretation
    lines.append("")
    lines.append("💡 INTERPRETATION:")

    if results["significant"]:
        if direction == "INCREASE":
            lines.append(f"  Disease shows INCREASED transcriptional chaos")
            lines.append(f"  Pattern consistent with: metabolic disorder, cancer, inflammation")
        else:
            lines.append(f"  Disease shows DECREASED transcriptional chaos")
            lines.append(f"  Pattern consistent with: neurodegeneration, protein-level pathology")
    else:
        lines.append(f"  No significant difference in transcriptional chaos")
        lines.append(f"  This could indicate: insufficient samples, or disease acts at protein level")

    return "\n".join(lines)

# ==============================================================
# OUTPUT
# ==============================================================

def generate_output(output_dir, dataset, results):
    """Generate output files (JSON + plots + HTML with per-sample details)."""
    import base64
    from io import BytesIO

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)

    # Save JSON
    output_data = {
        "dataset_name": dataset.get("dataset_name"),
        "source_file": dataset.get("source_file"),
        "format": dataset.get("format"),
        "generated_at": datetime.now().isoformat(),
        "results": {k: v for k, v in results.items()
                   if k not in ["sample_details", "control_entropy", "disease_entropy"]}
    }

    (output_dir / "results.json").write_text(
        json.dumps(output_data, indent=2, default=str),
        encoding="utf-8"
    )

    # Get sample details
    sample_details = results.get("sample_details", {})
    c = results.get("control_entropy", [])
    d = results.get("disease_entropy", [])

    # ===== PLOT 1: Main boxplot =====
    fig, ax = plt.subplots(figsize=(10, 6))
    if c and d:
        bp = ax.boxplot([c, d], labels=["Control", "Disease"], patch_artist=True)
        bp["boxes"][0].set_facecolor("#2ecc71")
        bp["boxes"][1].set_facecolor("#e74c3c")
        for box in bp["boxes"]:
            box.set_alpha(0.7)
    ax.set_ylabel("Shannon Entropy", fontsize=12)
    ax.set_title(f"Chaos Analysis: {dataset['dataset_name']}", fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    if "p_value" in results:
        p = results["p_value"]
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        direction = results.get("direction", "")
        ax.text(0.5, 0.95, f"{direction} | p = {p:.4f} {sig}", transform=ax.transAxes, ha="center", fontsize=11)
    plt.savefig(output_dir / "figures" / "boxplot.png", dpi=150, bbox_inches="tight")
    buf = BytesIO()
    plt.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close()
    buf.seek(0)
    boxplot_b64 = base64.b64encode(buf.read()).decode("utf-8")

    # ===== PLOT 2: Per-sample bar chart =====
    if sample_details:
        # Sort samples by label then by entropy
        sorted_samples = sorted(sample_details.items(),
                                key=lambda x: (x[1].get('label', ''), x[1].get('entropy', 0)))

        sample_ids = [s[0] for s in sorted_samples]
        entropies = [s[1].get('entropy', 0) for s in sorted_samples]
        colors = ['#2ecc71' if s[1].get('label') == 'CONTROL' else '#e74c3c' for s in sorted_samples]

        fig, ax = plt.subplots(figsize=(max(10, len(sample_ids) * 0.5), 6))
        bars = ax.bar(range(len(sample_ids)), entropies, color=colors, alpha=0.7)
        ax.set_xticks(range(len(sample_ids)))
        ax.set_xticklabels(sample_ids, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel("Shannon Entropy", fontsize=12)
        ax.set_title("Entropy per Sample", fontsize=14, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor='#2ecc71', alpha=0.7, label='Control'),
                          Patch(facecolor='#e74c3c', alpha=0.7, label='Disease')]
        ax.legend(handles=legend_elements, loc='upper right')

        plt.tight_layout()
        plt.savefig(output_dir / "figures" / "samples_bar.png", dpi=150, bbox_inches="tight")
        buf = BytesIO()
        plt.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close()
        buf.seek(0)
        samples_bar_b64 = base64.b64encode(buf.read()).decode("utf-8")
    else:
        samples_bar_b64 = ""

    # ===== Build sample cards HTML =====
    sample_cards_html = ""
    if sample_details:
        # Calculate min/max for visual bar
        all_ent = [s.get('entropy', 0) for s in sample_details.values() if s.get('entropy')]
        min_ent = min(all_ent) if all_ent else 0
        max_ent = max(all_ent) if all_ent else 1
        range_ent = max_ent - min_ent if max_ent != min_ent else 1

        for sample_id, details in sorted(sample_details.items(), key=lambda x: x[1].get('label', '')):
            label = details.get('label', 'Unknown')
            entropy_val = details.get('entropy', 0)
            n_genes = details.get('n_genes', 0)

            # Color and style based on label
            if label == 'CONTROL':
                border_color = '#2ecc71'
                label_bg = '#d4edda'
                label_color = '#155724'
            elif label == 'DISEASE':
                border_color = '#e74c3c'
                label_bg = '#f8d7da'
                label_color = '#721c24'
            else:
                border_color = '#6c757d'
                label_bg = '#e9ecef'
                label_color = '#495057'

            # Visual bar width (0-100%)
            bar_width = ((entropy_val - min_ent) / range_ent) * 100 if entropy_val else 0

            sample_cards_html += f"""
            <div class="sample-card" style="border-left: 4px solid {border_color};">
                <div class="sample-header">
                    <span class="sample-id">{sample_id}</span>
                    <span class="sample-label" style="background: {label_bg}; color: {label_color};">{label}</span>
                </div>
                <div class="sample-entropy">
                    <span class="entropy-value">{entropy_val:.6f}</span>
                    <span class="entropy-label">Shannon Entropy</span>
                </div>
                <div class="entropy-bar-container">
                    <div class="entropy-bar" style="width: {bar_width}%; background: {border_color};"></div>
                </div>
                <div class="sample-meta">Genes: {n_genes:,}</div>
            </div>
            """

    # ===== Build data table HTML =====
    data_table_html = """
    <table class="data-table">
        <thead>
            <tr>
                <th>Sample ID</th>
                <th>Group</th>
                <th>Entropy</th>
                <th>Genes</th>
                <th>Visual</th>
            </tr>
        </thead>
        <tbody>
    """

    if sample_details:
        for sample_id, details in sorted(sample_details.items(), key=lambda x: (x[1].get('label', ''), -x[1].get('entropy', 0))):
            label = details.get('label', 'Unknown')
            entropy_val = details.get('entropy', 0)
            n_genes = details.get('n_genes', 0)

            label_class = 'control' if label == 'CONTROL' else 'disease' if label == 'DISEASE' else ''
            bar_width = ((entropy_val - min_ent) / range_ent) * 100 if entropy_val else 0
            bar_color = '#2ecc71' if label == 'CONTROL' else '#e74c3c'

            data_table_html += f"""
            <tr>
                <td><strong>{sample_id}</strong></td>
                <td><span class="label-badge {label_class}">{label}</span></td>
                <td>{entropy_val:.6f}</td>
                <td>{n_genes:,}</td>
                <td>
                    <div class="mini-bar-container">
                        <div class="mini-bar" style="width: {bar_width}%; background: {bar_color};"></div>
                    </div>
                </td>
            </tr>
            """

    data_table_html += """
        </tbody>
    </table>
    """

    # ===== Generate HTML =====
    sig_class = "sig-yes" if results.get("significant") else "sig-no"
    sig_text = "Significant ✓" if results.get("significant") else "Not Significant"

    p_val = results.get('p_value')
    p_value_str = f"{p_val:.6f}" if isinstance(p_val, float) else "N/A"
    ctrl_mean = results.get('control_mean', 0)
    ctrl_std = results.get('control_std', 0)
    dis_mean = results.get('disease_mean', 0)
    dis_std = results.get('disease_std', 0)
    pct_change = results.get('percent_change', 0)

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>MK4 Report - {dataset['dataset_name']}</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{ font-family: -apple-system, Arial, sans-serif; margin: 0; background: #f5f7fa; }}
        .header {{ background: linear-gradient(135deg, #667eea, #764ba2); color: white; padding: 30px; text-align: center; }}
        .container {{ max-width: 1200px; margin: 20px auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ margin: 0; font-size: 2em; }}
        h2 {{ color: #667eea; border-bottom: 2px solid #e0e0e0; padding-bottom: 10px; margin-top: 40px; }}

        .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }}
        .stat-card {{ padding: 20px; border-radius: 8px; background: #f8f9fa; border-left: 4px solid #667eea; }}
        .stat-card.control {{ border-color: #2ecc71; }}
        .stat-card.disease {{ border-color: #e74c3c; }}
        .stat-card.result {{ border-color: #f39c12; }}
        .stat-value {{ font-size: 1.8em; font-weight: bold; color: #333; }}
        .stat-label {{ color: #666; font-size: 0.9em; margin-top: 5px; }}

        img {{ max-width: 100%; margin: 20px 0; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }}

        .sig-yes {{ background: #d4edda; color: #155724; padding: 5px 15px; border-radius: 20px; display: inline-block; }}
        .sig-no {{ background: #f8d7da; color: #721c24; padding: 5px 15px; border-radius: 20px; display: inline-block; }}

        .interpretation {{ background: #e8f4f8; padding: 20px; border-radius: 8px; margin: 20px 0; }}

        /* Sample cards */
        .samples-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 15px; margin: 20px 0; }}
        .sample-card {{ background: #fff; border: 1px solid #e0e0e0; border-radius: 8px; padding: 15px; }}
        .sample-header {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px; }}
        .sample-id {{ font-weight: bold; font-size: 1.1em; color: #333; }}
        .sample-label {{ padding: 3px 10px; border-radius: 12px; font-size: 0.8em; font-weight: bold; }}
        .sample-entropy {{ margin: 10px 0; }}
        .entropy-value {{ font-size: 1.5em; font-weight: bold; color: #667eea; }}
        .entropy-label {{ display: block; font-size: 0.8em; color: #888; }}
        .entropy-bar-container {{ height: 8px; background: #e9ecef; border-radius: 4px; overflow: hidden; }}
        .entropy-bar {{ height: 100%; border-radius: 4px; transition: width 0.3s; }}
        .sample-meta {{ font-size: 0.85em; color: #888; margin-top: 8px; }}

        /* Data table */
        .data-table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        .data-table th {{ background: #667eea; color: white; padding: 12px; text-align: left; }}
        .data-table td {{ padding: 10px 12px; border-bottom: 1px solid #e0e0e0; }}
        .data-table tr:hover {{ background: #f8f9fa; }}
        .label-badge {{ padding: 3px 10px; border-radius: 12px; font-size: 0.85em; font-weight: bold; }}
        .label-badge.control {{ background: #d4edda; color: #155724; }}
        .label-badge.disease {{ background: #f8d7da; color: #721c24; }}
        .mini-bar-container {{ width: 100px; height: 10px; background: #e9ecef; border-radius: 5px; overflow: hidden; }}
        .mini-bar {{ height: 100%; border-radius: 5px; }}

        .footer {{ text-align: center; color: #888; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>MK4 Biomarker Analysis</h1>
        <p>Dataset: {dataset['dataset_name']} | Format: {dataset.get('format', 'Unknown')} | Samples: {len(sample_details)}</p>
    </div>

    <div class="container">
        <h2>📊 Summary Statistics</h2>
        <div class="stats">
            <div class="stat-card control">
                <div class="stat-value">{results.get('n_control', 0)}</div>
                <div class="stat-label">Control Samples</div>
                <div style="margin-top: 10px; font-size: 0.9em;">
                    Mean: {ctrl_mean:.6f}<br>Std: ±{ctrl_std:.6f}
                </div>
            </div>
            <div class="stat-card disease">
                <div class="stat-value">{results.get('n_disease', 0)}</div>
                <div class="stat-label">Disease Samples</div>
                <div style="margin-top: 10px; font-size: 0.9em;">
                    Mean: {dis_mean:.6f}<br>Std: ±{dis_std:.6f}
                </div>
            </div>
            <div class="stat-card result">
                <div class="stat-value">{p_value_str}</div>
                <div class="stat-label">p-value (t-test)</div>
                <div style="margin-top: 10px;">
                    <span class="{sig_class}">{sig_text}</span>
                </div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{results.get('direction', 'N/A')}</div>
                <div class="stat-label">Direction</div>
                <div style="margin-top: 10px; font-size: 0.9em;">
                    Change: {pct_change:+.2f}%
                </div>
            </div>
        </div>

        <div class="interpretation">
            <h3>💡 Interpretation</h3>
            <p>
                {"Disease group shows <strong>INCREASED</strong> transcriptional chaos compared to controls. This pattern is consistent with metabolic disorders, cancer, or inflammatory conditions where cellular regulation is disrupted."
                 if results.get('direction') == 'INCREASE' and results.get('significant') else
                 "Disease group shows <strong>DECREASED</strong> transcriptional chaos compared to controls. This pattern is consistent with neurodegenerative conditions or protein-level pathologies."
                 if results.get('direction') == 'DECREASE' and results.get('significant') else
                 "No statistically significant difference in transcriptional chaos between groups. This could indicate insufficient sample size, high variability, or that the disease mechanism operates at the protein level rather than transcriptional level."}
            </p>
        </div>

        <h2>📈 Group Comparison</h2>
        <img src="data:image/png;base64,{boxplot_b64}" alt="Boxplot Comparison">

        <h2>📊 Per-Sample Entropy</h2>
        {"<img src='data:image/png;base64," + samples_bar_b64 + "' alt='Per-sample entropy'>" if samples_bar_b64 else "<p>No sample data available</p>"}

        <h2>🧬 Individual Sample Details</h2>
        <div class="samples-grid">
            {sample_cards_html}
        </div>

        <h2>📋 Data Table</h2>
        {data_table_html}

        <div class="footer">
            <p><strong>MK4 Biomarker Analysis Engine</strong></p>
            <p>Alexandria Dynamics | Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>"""

    (output_dir / "index.html").write_text(html, encoding="utf-8")

    return output_dir

# ==============================================================
# MAIN AUTO-ANALYZE
# ==============================================================

def auto_analyze(filepaths):
    """
    Main entry point - auto-detect and analyze any GEO file(s).

    Args:
        filepaths: Single file path or list of paths

    Returns:
        dict with all results
    """
    if isinstance(filepaths, (str, Path)):
        filepaths = [Path(filepaths)]
    else:
        filepaths = [Path(f) for f in filepaths]

    print("\n" + "=" * 70)
    print("  MK4 BIOMARKER - Universal Auto-Parser")
    print("  Alexandria Dynamics")
    print("=" * 70)

    # Detect formats
    soft_file = None
    matrix_file = None

    for fp in filepaths:
        print(f"\n📂 Checking: {fp.name}")
        fmt = detect_file_format(fp)
        print(f"   Format: {fmt}")

        if fmt == "SOFT":
            soft_file = fp
        elif fmt == "MATRIX":
            matrix_file = fp
        else:
            print(f"   ⚠️ Unknown format, skipping")

    # Parse based on what we have
    dataset = None

    if soft_file and has_expression_tables(soft_file):
        # SOFT with embedded expression - use directly
        print(f"\n✅ Using SOFT file with embedded expression tables")
        dataset = parse_soft_file(soft_file)

    elif soft_file and matrix_file:
        # SOFT metadata + Matrix data - combine them
        print(f"\n✅ Combining SOFT metadata with Matrix data")
        soft_data = parse_soft_file(soft_file)
        dataset = parse_matrix_file(matrix_file, soft_metadata=soft_data)
        dataset["combined_from"] = [soft_file.name, matrix_file.name]

    elif matrix_file:
        # Matrix only - try to infer labels from column names
        print(f"\n✅ Using Matrix file (inferring labels)")
        dataset = parse_matrix_file(matrix_file)

    elif soft_file:
        # SOFT without expression - can only get metadata
        print(f"\n⚠️ SOFT file without expression tables")
        print(f"   Need matching Matrix file for expression data")
        dataset = parse_soft_file(soft_file)

    else:
        print(f"\n❌ No valid files found")
        return {"error": "No valid GEO files found"}

    # Count samples
    samples = dataset.get("samples", {})
    n_with_expr = sum(1 for s in samples.values() if s.get("expression"))
    n_control = sum(1 for s in samples.values() if s.get("label") == "CONTROL")
    n_disease = sum(1 for s in samples.values() if s.get("label") == "DISEASE")

    print(f"\n📊 SAMPLE SUMMARY:")
    print(f"   Total samples: {len(samples)}")
    print(f"   With expression: {n_with_expr}")
    print(f"   Labeled CONTROL: {n_control}")
    print(f"   Labeled DISEASE: {n_disease}")
    print(f"   Unlabeled: {len(samples) - n_control - n_disease}")

    # Analyze
    if n_with_expr == 0:
        print(f"\n❌ No expression data found")
        return {"error": "No expression data found", "dataset": dataset}

    if n_control == 0 or n_disease == 0:
        print(f"\n⚠️ Missing control or disease samples")
        print(f"   Cannot perform statistical comparison")

    print(f"\n🔬 ANALYZING...")
    results = analyze_dataset(dataset)

    print(f"\n📊 ANALYSIS RESULTS:")
    print(interpret_results(results))

    # Generate output
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = RESULTS_DIR / f"{timestamp}_{dataset['dataset_name']}"
    generate_output(output_dir, dataset, results)

    print(f"\n✅ Output saved to: {output_dir.relative_to(PROJECT_ROOT)}")
    print(f"   Open index.html in browser to view report")
    print("=" * 70 + "\n")

    return {
        "dataset": dataset,
        "results": results,
        "output_dir": str(output_dir)
    }


# ==============================================================
# CLI
# ==============================================================

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python auto_parse_geo.py <file.gz> [file2.gz ...]")
        print("\nExamples:")
        print("  python auto_parse_geo.py GSE25724_family.soft.gz")
        print("  python auto_parse_geo.py GSE68086_family.soft.gz GSE68086_data_matrix.txt.gz")
        sys.exit(1)

    files = [Path(f) for f in sys.argv[1:]]

    # Check files exist
    for f in files:
        if not f.exists():
            print(f"Error: File not found: {f}")
            sys.exit(1)

    auto_analyze(files)
