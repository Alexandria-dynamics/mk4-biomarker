#!/usr/bin/env python3
"""
GEO SOFT Parser for MK4 Biomarker Pipeline
============================================
Parses GEO SOFT format files, extracts expression matrices
and auto-detects control/disease labels.

Usage:
    python scripts/parse_geo_soft.py

Input:  *.soft.gz files (GEO format)
Output: data/parsed/{GSE_ID}_expression.csv
        data/parsed/{GSE_ID}_labels.csv
"""

import gzip
import numpy as np
import pandas as pd
from pathlib import Path


# --- Configuration --------------------------------------------------------

INPUT_DIR = Path('.')  # Where SOFT files are -- adjust path!

OUTPUT_DIR = Path('data/parsed')

SOFT_FILES = [
    'GSE49036_family.soft.gz',   # Parkinson's (negative control)
    'GSE25724_family.soft.gz',   # T2D balanced (best result)
    'GSE76894_family.soft.gz',   # T2D imbalanced
    'GSE38642_family.soft.gz',   # T2D imbalanced
    'GSE70353_family.soft.gz',   # Unknown -- parser will tell us
]

# --- Label Detection Keywords ---------------------------------------------

DISEASE_KEYWORDS = [
    # Diabetes
    't2d', 'type 2 diabetes', 'type2diabetes',
    'diabetes mellitus type 2', 'diabetes mellitus, type 2',
    'diabetic',
    # MS
    'multiple sclerosis', 'relapsing', 'progressive ms',
    # Parkinson
    "parkinson", "parkinsons", "parkinson's",
    # Cancer
    'cancer', 'tumor', 'tumour', 'carcinoma',
    'melanoma', 'sarcoma', 'leukemia', 'lymphoma',
    # Generic
    'affected', 'case', 'patient',
]

CONTROL_KEYWORDS = [
    'control', 'healthy', 'normal', 'unaffected',
    'wild-type', 'wildtype',
    'non-diabetic', 'non diabetic', 'nondiabetic',
    'no disease',
]


# --- Parser ---------------------------------------------------------------

def parse_soft(filepath):
    """
    Parse GEO SOFT file. Extract metadata + expression per sample.

    Returns
    -------
    series_id : str
    samples_meta : dict  {gsm: {title, characteristics, ...}}
    samples_expr : dict  {gsm: {gene_id: value}}
    """
    filepath = Path(filepath)
    print(f"\n{'=' * 60}")
    print(f"  Parsing: {filepath.name}")
    print(f"  Size:    {filepath.stat().st_size / 1e6:.1f} MB")
    print(f"{'=' * 60}")

    opener = gzip.open if str(filepath).endswith('.gz') else open

    series_id = None
    samples_meta = {}
    samples_expr = {}
    current_gsm = None
    in_table = False
    skip_section = False

    with opener(filepath, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')

            # New section markers
            if line.startswith('^SERIES'):
                series_id = line.split(' = ', 1)[1].strip()
                skip_section = False
                continue

            if line.startswith('^SAMPLE'):
                current_gsm = line.split(' = ', 1)[1].strip()
                samples_meta[current_gsm] = {}
                samples_expr[current_gsm] = {}
                in_table = False
                skip_section = False
                continue

            if line.startswith('^PLATFORM') or line.startswith('^ENTITY'):
                skip_section = True
                current_gsm = None
                in_table = False
                continue

            if skip_section:
                continue

            # Sample metadata
            if line.startswith('!Sample_') and current_gsm:
                if 'table_begin' in line:
                    in_table = True
                    continue
                if 'table_end' in line:
                    in_table = False
                    continue

                if ' = ' in line:
                    key, val = line.split(' = ', 1)
                    key = key.replace('!Sample_', '').strip()
                    val = val.strip('"')

                    if key in ('title', 'geo_accession', 'description',
                               'source_organism', 'status',
                               'organism_ch1', 'tissue'):
                        samples_meta[current_gsm][key] = val

                    if key.startswith('characteristics_ch'):
                        samples_meta[current_gsm].setdefault(
                            'characteristics', []
                        ).append(val)
                continue

            # Expression table data
            if in_table and current_gsm:
                if line.startswith('ID_REF'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    gene_id = parts[0].strip()
                    try:
                        value = float(parts[1].strip())
                        samples_expr[current_gsm][gene_id] = value
                    except ValueError:
                        pass

    n_with_expr = sum(1 for g in samples_expr.values() if len(g) > 0)
    print(f"  Series:       {series_id}")
    print(f"  Samples:      {len(samples_meta)}")
    print(f"  With expr:    {n_with_expr}")
    if n_with_expr > 0:
        first_gsm = next(g for g in samples_expr if len(samples_expr[g]) > 0)
        print(f"  Genes/sample: {len(samples_expr[first_gsm])}")

    return series_id, samples_meta, samples_expr


# --- Label Detection ------------------------------------------------------

def detect_labels(samples_meta):
    """
    Auto-detect control(0) / disease(1) from characteristics.

    Returns {gsm: 0 | 1 | -1}
    -1 = could not determine
    """
    labels = {}

    for gsm, meta in samples_meta.items():
        chars = meta.get('characteristics', [])
        title = meta.get('title', '')
        desc = meta.get('description', '')

        text = ' '.join(chars + [title, desc]).lower()

        hit_disease = any(kw in text for kw in DISEASE_KEYWORDS)
        hit_control = any(kw in text for kw in CONTROL_KEYWORDS)

        if hit_control and not hit_disease:
            labels[gsm] = 0
        elif hit_disease and not hit_control:
            labels[gsm] = 1
        else:
            labels[gsm] = -1

    return labels


# --- Build & Save ---------------------------------------------------------

def build_and_save(series_id, samples_meta, samples_expr, labels):
    """
    Build expression DataFrame, save CSV files.
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    valid = [gsm for gsm in samples_expr if len(samples_expr[gsm]) > 0]

    if not valid:
        print(f"  No expression data -- skipping {series_id}")
        return None

    gene_sets = [set(samples_expr[g].keys()) for g in valid]
    common = sorted(set.intersection(*gene_sets))
    union = sorted(set.union(*gene_sets))

    if len(common) >= 100:
        genes = common
        print(f"  Using COMMON genes: {len(common)}")
    else:
        genes = union
        print(f"  Using UNION genes:  {len(union)} (common only {len(common)})")

    # Expression matrix
    matrix = np.full((len(valid), len(genes)), np.nan)
    for i, gsm in enumerate(valid):
        for j, gene in enumerate(genes):
            matrix[i, j] = samples_expr[gsm].get(gene, np.nan)

    expr_df = pd.DataFrame(matrix, index=valid, columns=genes)

    # Labels
    label_rows = []
    for gsm in valid:
        meta = samples_meta.get(gsm, {})
        label_rows.append({
            'sample_id': gsm,
            'label': labels.get(gsm, -1),
            'title': meta.get('title', ''),
            'characteristics': '; '.join(meta.get('characteristics', []))
        })
    labels_df = pd.DataFrame(label_rows)

    # Save
    expr_path = OUTPUT_DIR / f"{series_id}_expression.csv"
    labels_path = OUTPUT_DIR / f"{series_id}_labels.csv"

    expr_df.to_csv(expr_path)
    labels_df.to_csv(labels_path, index=False)

    print(f"\n  {expr_path}")
    print(f"       {expr_df.shape[0]} samples x {expr_df.shape[1]} genes")
    print(f"  {labels_path}")

    # Label summary
    print(f"\n  {'─' * 56}")
    print(f"  {'Sample':<12} {'Label':<10} Title / Characteristics")
    print(f"  {'─' * 56}")

    for _, row in labels_df.iterrows():
        lbl = row['label']
        lbl_str = 'CONTROL' if lbl == 0 else 'DISEASE' if lbl == 1 else '???'
        info = row['title'] if row['title'] else row['characteristics']
        print(f"  {row['sample_id']:<12} [{lbl_str:>7}]  {info[:55]}")

    print(f"  {'─' * 56}")

    n0 = (labels_df['label'] == 0).sum()
    n1 = (labels_df['label'] == 1).sum()
    nu = (labels_df['label'] == -1).sum()

    print(f"  Control:  {n0}")
    print(f"  Disease:  {n1}")
    if nu > 0:
        print(f"  UNKNOWN: {nu} samples -- REVIEW MANUALLY!")

    balance = n1 / (n0 + n1) if (n0 + n1) > 0 else 0
    print(f"  Balance:  {balance * 100:.1f}% disease")

    if 0.3 <= balance <= 0.7:
        print(f"  Well balanced")
    else:
        print(f"  Imbalanced -- MK4 sensitivity may be affected")

    return expr_df, labels_df


# --- Main -----------------------------------------------------------------

def main():
    print("\n" + "=" * 60)
    print("   GEO SOFT PARSER -- MK4 Biomarker Pipeline")
    print("   Alexandria Dynamics")
    print("=" * 60)

    all_results = {}

    for fname in SOFT_FILES:
        fpath = INPUT_DIR / fname

        if not fpath.exists():
            print(f"\n  NOT FOUND: {fpath} -- skipping")
            continue

        series_id, meta, expr = parse_soft(fpath)
        labels = detect_labels(meta)
        result = build_and_save(series_id, meta, expr, labels)

        if result:
            all_results[series_id] = result

    # Grand summary
    print("\n\n" + "=" * 60)
    print("   PARSING COMPLETE -- GRAND SUMMARY")
    print("=" * 60)

    if not all_results:
        print("  No datasets parsed. Check INPUT_DIR and file paths.")
        print(f"  Current INPUT_DIR: {INPUT_DIR.resolve()}")
        print("=" * 60)
        return

    print(f"  {'Dataset':<12} {'Samples':<10} {'Genes':<8} "
          f"{'Ctrl':<6} {'Dis':<6} {'Unk':<5} Balance")
    print(f"  {'─' * 60}")

    for sid, (edf, ldf) in sorted(all_results.items()):
        n0 = (ldf['label'] == 0).sum()
        n1 = (ldf['label'] == 1).sum()
        nu = (ldf['label'] == -1).sum()
        bal = n1 / (n0 + n1) * 100 if (n0 + n1) > 0 else 0
        unk_flag = " !" if nu > 0 else ""
        print(f"  {sid:<12} {len(ldf):<10} {edf.shape[1]:<8} "
              f"{n0:<6} {n1:<6} {nu:<5} {bal:.1f}%{unk_flag}")

    print(f"\n  Output directory: {OUTPUT_DIR}/")
    print(f"  Files: *_expression.csv, *_labels.csv")
    print("=" * 60)

    print("\n  NEXT: Load into MK4:")
    print("    import pandas as pd")
    print("    expr = pd.read_csv('data/parsed/GSE25724_expression.csv', index_col=0)")
    print("    labels = pd.read_csv('data/parsed/GSE25724_labels.csv')")
    print("=" * 60)


if __name__ == '__main__':
    main()
