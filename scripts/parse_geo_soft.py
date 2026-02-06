#!/usr/bin/env python3
"""
MK4 Biomarker - GEO Parser
Alexandria Dynamics

Podporuje:
1) GEO SOFT (.soft / .soft.gz)
2) GEO MATRIX (.txt / .txt.gz)

Detekce formátu je automatická.
"""

import gzip
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.stats import entropy, ttest_ind

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

def detect_file_type(filepath):
    filepath = Path(filepath)

    if filepath.suffix == ".gz":
        f = gzip.open(filepath, "rt")
    else:
        f = open(filepath, "r")

    try:
        lines = [f.readline() for _ in range(50)]

        if any(line.startswith("^SAMPLE") or line.startswith("^SERIES") for line in lines):
            return "soft"

        if "\t" in lines[0] and (
            lines[1].startswith("ENSG")
            or lines[1].startswith("ILMN")
            or lines[1].startswith("ID")
        ):
            return "matrix"

        return None
    finally:
        f.close()

# ==============================================================
# LABEL INFERENCE (FIX)
# ==============================================================

def infer_label_from_meta(meta_text: str):
    t = meta_text.lower()

    control_kw = [
        "control", "normal", "healthy", "non-diabetic",
        "nondiabetic", "wild type", "wt", "sham",
        "untreated", "parental"
    ]

    disease_kw = [
        "disease", "patient", "diabetic", "diabetes",
        "tumor", "cancer", "multiple sclerosis",
        "ms", "parkinson", "alzheimer",
        "treated", "resistant"
    ]

    if any(k in t for k in control_kw):
        return "CONTROL"
    if any(k in t for k in disease_kw):
        return "DISEASE"
    return None

# ==============================================================
# SOFT PARSER (OPRAVENÝ)
# ==============================================================

def parse_soft_file(filepath):
    filepath = Path(filepath)
    dataset_name = filepath.name.replace("_family.soft.gz", "").replace("_family.soft", "")

    print("=" * 70)
    print(f"📂 Parsing SOFT: {filepath.name}")
    print("=" * 70)

    samples = {}
    current_sample = None
    in_table = False
    current_data = {}
    meta_buf = []

    if filepath.suffix == ".gz":
        f = gzip.open(filepath, "rt")
    else:
        f = open(filepath, "r")

    try:
        for line in f:
            line = line.strip()

            # ---- New sample ----
            if line.startswith("^SAMPLE"):
                current_sample = line.split("=")[1].strip()
                samples[current_sample] = {
                    "label": None,
                    "expression": {}
                }
                meta_buf = []

            # ---- Collect metadata ----
            elif current_sample and (
                line.startswith("!Sample_title")
                or line.startswith("!Sample_source_name_ch1")
                or line.startswith("!Sample_characteristics_ch1")
                or line.startswith("!Sample_description")
            ):
                meta_buf.append(line)

            # ---- Legacy special case ----
            elif current_sample and "disease state:" in line.lower():
                meta_buf.append(line)

            # ---- Table begin ----
            elif line == "!sample_table_begin":
                in_table = True
                current_data = {}

            # ---- Table end ----
            elif line == "!sample_table_end":
                in_table = False

                if current_sample and current_data:
                    samples[current_sample]["expression"] = current_data

                if samples[current_sample]["label"] is None:
                    lbl = infer_label_from_meta(" | ".join(meta_buf))
                    samples[current_sample]["label"] = lbl

            # ---- Table rows ----
            elif in_table and "\t" in line and not line.startswith("ID_REF"):
                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        current_data[parts[0]] = float(parts[1])
                    except ValueError:
                        pass

    finally:
        f.close()

    valid_samples = {
        k: v for k, v in samples.items()
        if v["label"] is not None and len(v["expression"]) > 0
    }

    print(f"✅ Found {len(valid_samples)} valid samples")
    print(f"   Controls: {sum(1 for v in valid_samples.values() if v['label']=='CONTROL')}")
    print(f"   Disease:  {sum(1 for v in valid_samples.values() if v['label']=='DISEASE')}")

    return {
        "dataset_name": dataset_name,
        "samples": valid_samples
    }

# ==============================================================
# MATRIX PARSER
# ==============================================================

def parse_matrix_file(filepath):
    filepath = Path(filepath)
    dataset_name = filepath.name.replace(".txt.gz", "").replace(".txt", "")

    print("=" * 70)
    print(f"📂 Parsing MATRIX: {filepath.name}")
    print("=" * 70)

    if filepath.suffix == ".gz":
        df = pd.read_csv(filepath, sep="\t", compression="gzip")
    else:
        df = pd.read_csv(filepath, sep="\t")

    gene_col = df.columns[0]
    sample_cols = df.columns[1:]

    samples = {}

    for col in sample_cols:
        label = "CONTROL"
        col_l = col.lower()

        if any(x in col_l for x in ["disease", "tumor", "cancer", "patient", "treated"]):
            label = "DISEASE"

        expr = pd.to_numeric(df[col], errors="coerce")
        expr = expr.dropna()

        samples[col] = {
            "label": label,
            "expression": dict(zip(df[gene_col], expr))
        }

    print(f"✅ Found {len(samples)} samples")

    return {
        "dataset_name": dataset_name,
        "samples": samples
    }

# ==============================================================
# CHAOS ANALYSIS
# ==============================================================

def compute_sample_entropy(expression_dict):
    values = np.array(list(expression_dict.values()), dtype=float)
    values = values[values > 0]
    if len(values) == 0:
        return np.nan
    values = values / np.sum(values)
    return entropy(values)

def analyze_chaos(dataset):
    entropies = {"CONTROL": [], "DISEASE": []}

    for s in dataset["samples"].values():
        e = compute_sample_entropy(s["expression"])
        if not np.isnan(e):
            entropies[s["label"]].append(e)

    if len(entropies["CONTROL"]) > 1 and len(entropies["DISEASE"]) > 1:
        t, p = ttest_ind(entropies["CONTROL"], entropies["DISEASE"], equal_var=False)
    else:
        p = None

    return entropies, p

# ==============================================================
# MAIN
# ==============================================================

def run(filepath):
    ftype = detect_file_type(filepath)

    if ftype == "soft":
        dataset = parse_soft_file(filepath)
    elif ftype == "matrix":
        dataset = parse_matrix_file(filepath)
    else:
        raise ValueError("Unknown file format")

    entropies, pval = analyze_chaos(dataset)

    print("\n📊 CHAOS RESULTS")
    print(f"Controls: {len(entropies['CONTROL'])}")
    print(f"Disease:  {len(entropies['DISEASE'])}")
    print(f"p-value:  {pval}")

    return dataset, entropies, pval


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python parse_geo_soft.py <file>")
        sys.exit(1)

    run(sys.argv[1])
# ==============================================================
# GUI COMPAT LAYER
# ==============================================================

import json
import matplotlib.pyplot as plt

def parse_file(filepath):
    """
    Unified entrypoint expected by gui.py
    Returns dataset dict
    """
    ftype = detect_file_type(filepath)
    if ftype == "soft":
        return parse_soft_file(filepath)
    if ftype == "matrix":
        return parse_matrix_file(filepath)
    return None


def compute_chaos_metrics(dataset):
    entropies = {"CONTROL": [], "DISEASE": []}

    for s in dataset["samples"].values():
        e = compute_sample_entropy(s["expression"])
        if not np.isnan(e):
            entropies[s["label"]].append(float(e))

    stats = {}
    c = entropies["CONTROL"]
    d = entropies["DISEASE"]

    if len(c) >= 2 and len(d) >= 2:
        t, p = ttest_ind(c, d, equal_var=False)
        stats = {
            "control_mean": float(np.mean(c)),
            "control_std": float(np.std(c, ddof=1)),
            "disease_mean": float(np.mean(d)),
            "disease_std": float(np.std(d, ddof=1)),
            "p_value": float(p),
            "n_control": len(c),
            "n_disease": len(d),
        }
    else:
        stats = {"error": "insufficient samples"}

    return {"entropies": entropies, "statistics": stats}


def create_output_directory(dataset_name):
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    out = RESULTS_DIR / f"run_{ts}_{dataset_name}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(exist_ok=True)
    return out


def save_results(output_dir, dataset, chaos_results):
    (output_dir / "chaos_results.json").write_text(
        json.dumps(chaos_results, indent=2), encoding="utf-8"
    )
    write_index_html(output_dir)


def generate_plot(output_dir, chaos_results):
    c = chaos_results["entropies"]["CONTROL"]
    d = chaos_results["entropies"]["DISEASE"]

    plt.figure()
    plt.boxplot([c, d], labels=["CONTROL", "DISEASE"])
    plt.title("Chaos entropy")
    plt.ylabel("Entropy")
    plt.savefig(output_dir / "figures" / "chaos_boxplot.png", dpi=150)
    plt.close()

    write_index_html(output_dir)


def write_index_html(output_dir):
    figs = sorted((output_dir / "figures").glob("*.png"))
    imgs = "\n".join(
        f'<img src="figures/{f.name}" style="max-width:100%;margin:10px;">'
        for f in figs
    )

    html = f"""<!doctype html>
<html>
<body>
<h2>MK4 Analysis Output</h2>
{imgs}
</body>
</html>
"""
    (output_dir / "index.html").write_text(html, encoding="utf-8")
