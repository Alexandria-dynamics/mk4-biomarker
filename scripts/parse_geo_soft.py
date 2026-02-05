#!/usr/bin/env python3
"""
MK4 Biomarker - GEO SOFT Parser
Alexandria Dynamics

Parses GEO SOFT files (.gz or uncompressed) and computes chaos metrics.
Results saved to timestamped directories - never overwrites old results.
"""

import gzip
import json
import os
from pathlib import Path
from datetime import datetime
import numpy as np
from scipy.stats import entropy, ttest_ind

# ═══════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════

# Paths relative to project root
PROJECT_ROOT = Path(__file__).parent.parent
INPUT_DIR = PROJECT_ROOT / "data" / "input"
RESULTS_DIR = PROJECT_ROOT / "results"

# Create directories if they don't exist
INPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════
# PARSER
# ═══════════════════════════════════════════════════════════════

def parse_soft_file(filepath):
    """
    Parse GEO SOFT file (supports .gz compression).

    Returns:
        dict: {
            'dataset_name': str,
            'samples': {
                'GSM123': {'label': 'CONTROL', 'expression': {...}},
                ...
            }
        }
    """
    filepath = Path(filepath)
    dataset_name = filepath.stem.replace('_family.soft', '').replace('.soft', '').replace('.gz', '')

    print(f"\n{'='*70}")
    print(f"Parsing: {filepath.name}")
    print(f"{'='*70}")

    samples = {}
    current_sample = None
    in_table = False
    current_data = {}

    # Open with gzip if .gz, otherwise normal open
    if filepath.suffix == '.gz':
        f = gzip.open(filepath, 'rt')
    else:
        f = open(filepath, 'r')

    try:
        for line in f:
            line = line.strip()

            # New sample
            if line.startswith('^SAMPLE'):
                current_sample = line.split('=')[1].strip()
                samples[current_sample] = {'label': None, 'expression': {}}

            # Get disease label
            elif current_sample and 'disease state:' in line.lower():
                if 'non-diabetic' in line.lower() or 'control' in line.lower() or 'normal' in line.lower():
                    samples[current_sample]['label'] = 'CONTROL'
                elif 'diabetic' in line.lower() or 'diabetes' in line.lower() or 'disease' in line.lower():
                    samples[current_sample]['label'] = 'DISEASE'

            # Table data
            elif line == '!sample_table_begin':
                in_table = True
                current_data = {}

            elif line == '!sample_table_end':
                in_table = False
                if current_sample and current_data:
                    samples[current_sample]['expression'] = current_data

            elif in_table and '\t' in line and not line.startswith('ID_REF'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    try:
                        current_data[parts[0]] = float(parts[1])
                    except ValueError:
                        pass

    finally:
        f.close()

    # Filter valid samples
    valid_samples = {
        k: v for k, v in samples.items()
        if v['label'] is not None and len(v['expression']) > 0
    }

    print(f"Found {len(valid_samples)} valid samples")
    print(f"   Controls: {sum(1 for v in valid_samples.values() if v['label']=='CONTROL')}")
    print(f"   Disease:  {sum(1 for v in valid_samples.values() if v['label']=='DISEASE')}")

    return {
        'dataset_name': dataset_name,
        'samples': valid_samples
    }

# ═══════════════════════════════════════════════════════════════
# CHAOS ANALYSIS
# ═══════════════════════════════════════════════════════════════

def compute_chaos_metrics(dataset):
    """
    Compute Shannon entropy (chaos metric) for each sample.
    """
    samples = dataset['samples']

    # Get common genes
    all_genes = [set(s['expression'].keys()) for s in samples.values()]
    common_genes = sorted(set.intersection(*all_genes))

    print(f"\nCommon genes: {len(common_genes)}")

    control_entropy = []
    disease_entropy = []
    sample_details = {}

    for gsm, data in samples.items():
        label = data['label']

        # Get expression values
        values = np.array([data['expression'][g] for g in common_genes])
        values_norm = values / values.sum()

        # Shannon entropy
        ent = entropy(values_norm, base=2)

        sample_details[gsm] = {
            'label': label,
            'entropy': float(ent),
            'num_genes': len(common_genes)
        }

        if label == 'CONTROL':
            control_entropy.append(ent)
        elif label == 'DISEASE':
            disease_entropy.append(ent)

    # Statistics
    if len(control_entropy) > 0 and len(disease_entropy) > 0:
        t_stat, p_value = ttest_ind(control_entropy, disease_entropy)

        stats = {
            'control_mean': float(np.mean(control_entropy)),
            'control_std': float(np.std(control_entropy)),
            'disease_mean': float(np.mean(disease_entropy)),
            'disease_std': float(np.std(disease_entropy)),
            'difference': float(np.mean(disease_entropy) - np.mean(control_entropy)),
            'percent_change': float((np.mean(disease_entropy) - np.mean(control_entropy)) / np.mean(control_entropy) * 100),
            't_statistic': float(t_stat),
            'p_value': float(p_value),
            'significant': bool(p_value < 0.05)
        }
    else:
        stats = {'error': 'Insufficient samples'}

    return {
        'control_entropy': [float(x) for x in control_entropy],
        'disease_entropy': [float(x) for x in disease_entropy],
        'statistics': stats,
        'sample_details': sample_details
    }

# ═══════════════════════════════════════════════════════════════
# OUTPUT
# ═══════════════════════════════════════════════════════════════

def create_output_directory(dataset_name):
    """Create timestamped output directory."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = RESULTS_DIR / f"{timestamp}_{dataset_name}"

    (output_dir / "json").mkdir(parents=True, exist_ok=True)
    (output_dir / "plots").mkdir(parents=True, exist_ok=True)

    print(f"\nOutput: {output_dir.relative_to(PROJECT_ROOT)}")
    return output_dir

def save_results(output_dir, dataset, chaos_results):
    """Save JSON results."""

    # Chaos results
    with open(output_dir / "json" / "chaos_results.json", 'w') as f:
        json.dump(chaos_results, f, indent=2)

    # Metadata
    metadata = {
        'dataset_name': dataset['dataset_name'],
        'timestamp': datetime.now().isoformat(),
        'num_samples': len(dataset['samples']),
        'num_controls': sum(1 for s in dataset['samples'].values() if s['label'] == 'CONTROL'),
        'num_disease': sum(1 for s in dataset['samples'].values() if s['label'] == 'DISEASE'),
        'analysis_type': 'Shannon_entropy_chaos_metric'
    }

    with open(output_dir / "json" / "metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

    # Sample info
    sample_info = {
        gsm: {
            'label': data['label'],
            'num_genes': len(data['expression'])
        }
        for gsm, data in dataset['samples'].items()
    }

    with open(output_dir / "json" / "sample_info.json", 'w') as f:
        json.dump(sample_info, f, indent=2)

    print(f"Saved JSON files")

def generate_plots(output_dir, chaos_results):
    """Generate plots."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        control = chaos_results['control_entropy']
        disease = chaos_results['disease_entropy']
        stats = chaos_results['statistics']

        if 'error' in stats:
            print("Skipping plots - insufficient data")
            return

        fig, ax = plt.subplots(figsize=(10, 6))

        bp = ax.boxplot([control, disease],
                        tick_labels=['Control', 'Disease'],
                        patch_artist=True)

        for patch, color in zip(bp['boxes'], ['#2ecc71', '#e74c3c']):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_ylabel('Shannon Entropy', fontsize=12)
        ax.set_title('Transcriptional Chaos: Control vs Disease', fontsize=14)
        ax.grid(axis='y', alpha=0.3)

        sig = "***" if stats['p_value'] < 0.001 else "**" if stats['p_value'] < 0.01 else "*" if stats['p_value'] < 0.05 else "ns"
        ax.text(0.5, 0.95, f"p = {stats['p_value']:.4f} {sig}",
                transform=ax.transAxes, ha='center')

        plt.tight_layout()
        plt.savefig(output_dir / "plots" / "chaos_distribution.png", dpi=300)
        plt.close()

        print(f"Saved plots")

    except ImportError:
        print("matplotlib not installed - skipping plots")

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    """Main pipeline."""

    print("\n" + "="*70)
    print("  MK4 BIOMARKER - RNA Chaos Analysis")
    print("  Alexandria Dynamics")
    print("="*70)

    # Find SOFT files
    soft_files = list(INPUT_DIR.glob("*.gz")) + list(INPUT_DIR.glob("*.soft"))

    if not soft_files:
        print(f"\nNo files in: {INPUT_DIR.relative_to(PROJECT_ROOT)}")
        print(f"\nPlace .soft.gz files in: data/input/")
        return

    print(f"\nFound {len(soft_files)} file(s):")
    for f in soft_files:
        print(f"   - {f.name}")

    # Process each file
    for filepath in soft_files:
        try:
            dataset = parse_soft_file(filepath)

            if len(dataset['samples']) == 0:
                print(f"No valid samples found")
                continue

            chaos_results = compute_chaos_metrics(dataset)
            output_dir = create_output_directory(dataset['dataset_name'])
            save_results(output_dir, dataset, chaos_results)
            generate_plots(output_dir, chaos_results)

            # Summary
            if 'error' not in chaos_results['statistics']:
                stats = chaos_results['statistics']
                print(f"\nSUMMARY")
                print(f"{'='*70}")
                print(f"  Control: {stats['control_mean']:.6f} +/- {stats['control_std']:.6f}")
                print(f"  Disease: {stats['disease_mean']:.6f} +/- {stats['disease_std']:.6f}")
                print(f"  Diff:    {stats['difference']:+.6f} ({stats['percent_change']:+.2f}%)")
                print(f"  p-value: {stats['p_value']:.6f}")
                print(f"  Signif:  {'YES' if stats['significant'] else 'NO'}")
                print(f"{'='*70}")

        except Exception as e:
            print(f"\nError: {e}")
            continue

    print(f"\nComplete! Results in: {RESULTS_DIR.relative_to(PROJECT_ROOT)}/")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
