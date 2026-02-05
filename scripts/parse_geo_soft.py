#!/usr/bin/env python3
"""
MK4 Biomarker - GEO SOFT Parser
Alexandria Dynamics

Parses GEO SOFT files (.gz or uncompressed) and computes chaos metrics.
Results saved to timestamped directories - never overwrites old results.
"""

import base64
import gzip
import json
from io import BytesIO
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
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nOutput: {output_dir.relative_to(PROJECT_ROOT)}")
    return output_dir


def generate_html_report(output_dir, dataset, chaos_results):
    """
    Generate single self-contained HTML report with embedded plot and data.

    Args:
        output_dir: Path to output directory
        dataset: Parsed dataset dict
        chaos_results: Analysis results dict

    Returns:
        Path to generated index.html
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    control = chaos_results['control_entropy']
    disease = chaos_results['disease_entropy']
    stats = chaos_results['statistics']

    # Generate plot in memory
    fig, ax = plt.subplots(figsize=(10, 6))

    if 'error' not in stats:
        bp = ax.boxplot([control, disease],
                        tick_labels=['Control', 'Disease'],
                        patch_artist=True)

        for patch, color in zip(bp['boxes'], ['#2ecc71', '#e74c3c']):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax.set_ylabel('Shannon Entropy', fontsize=12)
        ax.set_title('Transcriptional Chaos: Control vs Disease',
                      fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)

        sig = "***" if stats['p_value'] < 0.001 else \
              "**" if stats['p_value'] < 0.01 else \
              "*" if stats['p_value'] < 0.05 else "ns"

        ax.text(0.5, 0.95, f"p = {stats['p_value']:.4f} {sig}",
                transform=ax.transAxes, ha='center')

    # Save plot to base64
    buf = BytesIO()
    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    buf.seek(0)
    plot_b64 = base64.b64encode(buf.read()).decode('utf-8')

    # Build sample table rows
    sample_rows = []
    for gsm, details in chaos_results['sample_details'].items():
        label_color = '#2ecc71' if details['label'] == 'CONTROL' else '#e74c3c'
        row = (
            '<tr style="border-bottom: 1px solid #e0e0e0;">'
            f'<td style="padding: 12px;">{gsm}</td>'
            f'<td style="padding: 12px;">'
            f'<span style="color: {label_color}; font-weight: bold;">'
            f'{details["label"]}</span></td>'
            f'<td style="padding: 12px; text-align: right;">{details["entropy"]:.6f}</td>'
            f'<td style="padding: 12px; text-align: right;">{details["num_genes"]:,}</td>'
            '</tr>'
        )
        sample_rows.append(row)

    # Prepare template variables
    if 'error' not in stats:
        sig_class = 'sig-yes' if stats['significant'] else 'sig-no'
        sig_text = 'Significant' if stats['significant'] else 'Not Significant'
        stats_html = f"""
                <div class="stats-grid">
                    <div class="stat-card control">
                        <h3>Control Group</h3>
                        <div class="value">{stats['control_mean']:.6f}</div>
                        <div class="label">Shannon Entropy (n={len(control)})</div>
                        <div style="margin-top: 10px; color: #666;">
                            +/- {stats['control_std']:.6f}
                        </div>
                    </div>

                    <div class="stat-card disease">
                        <h3>Disease Group</h3>
                        <div class="value">{stats['disease_mean']:.6f}</div>
                        <div class="label">Shannon Entropy (n={len(disease)})</div>
                        <div style="margin-top: 10px; color: #666;">
                            +/- {stats['disease_std']:.6f}
                        </div>
                    </div>

                    <div class="stat-card result">
                        <h3>Statistical Test</h3>
                        <div class="value">{stats['p_value']:.6f}</div>
                        <div class="label">p-value (t-test)</div>
                        <div style="margin-top: 10px;">
                            Difference: <strong>{stats['difference']:+.6f}</strong> ({stats['percent_change']:+.2f}%)
                        </div>
                        <span class="significance {sig_class}">
                            {sig_text}
                        </span>
                    </div>
                </div>"""
    else:
        stats_html = '<p>Insufficient samples for statistical analysis.</p>'

    timestamp_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    json_data = json.dumps(chaos_results, indent=2)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MK4 Analysis Report - {dataset['dataset_name']}</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f7fa;
            padding: 20px;
        }}

        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }}

        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}

        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}

        .content {{
            padding: 40px;
        }}

        .section {{
            margin-bottom: 40px;
        }}

        .section h2 {{
            color: #667eea;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #e0e0e0;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}

        .stat-card {{
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            padding: 20px;
            background: #fafafa;
        }}

        .stat-card h3 {{
            color: #555;
            margin-bottom: 15px;
            font-size: 1.2em;
        }}

        .stat-card .value {{
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
            margin: 10px 0;
        }}

        .stat-card .label {{
            color: #888;
            font-size: 0.9em;
        }}

        .control {{ border-left: 4px solid #2ecc71; }}
        .disease {{ border-left: 4px solid #e74c3c; }}
        .result {{ border-left: 4px solid #f39c12; }}

        .plot-container {{
            text-align: center;
            margin: 30px 0;
        }}

        .plot-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}

        .json-section {{
            margin-top: 40px;
        }}

        details {{
            background: #f8f9fa;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            padding: 15px;
        }}

        summary {{
            cursor: pointer;
            font-weight: bold;
            color: #667eea;
            font-size: 1.1em;
            user-select: none;
        }}

        summary:hover {{
            color: #764ba2;
        }}

        pre {{
            background: #272822;
            color: #f8f8f2;
            padding: 20px;
            border-radius: 8px;
            overflow-x: auto;
            margin-top: 15px;
            font-size: 0.9em;
            line-height: 1.5;
        }}

        .footer {{
            background: #f8f9fa;
            padding: 30px;
            text-align: center;
            color: #888;
            border-top: 1px solid #e0e0e0;
        }}

        .significance {{
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-weight: bold;
            margin-top: 10px;
        }}

        .sig-yes {{
            background: #d4edda;
            color: #155724;
        }}

        .sig-no {{
            background: #f8d7da;
            color: #721c24;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
        }}

        @media print {{
            body {{ background: white; }}
            .container {{ box-shadow: none; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>MK4 Biomarker Analysis</h1>
            <p>Dataset: {dataset['dataset_name']} | Generated: {timestamp_str}</p>
        </div>

        <div class="content">
            <div class="section">
                <h2>Results Summary</h2>
                {stats_html}
            </div>

            <div class="section">
                <h2>Visualization</h2>
                <div class="plot-container">
                    <img src="data:image/png;base64,{plot_b64}" alt="Chaos Distribution Plot">
                </div>
            </div>

            <div class="section">
                <h2>Sample Details</h2>
                <table>
                    <thead>
                        <tr style="background: #f8f9fa;">
                            <th style="padding: 12px; text-align: left; border-bottom: 2px solid #e0e0e0;">Sample ID</th>
                            <th style="padding: 12px; text-align: left; border-bottom: 2px solid #e0e0e0;">Group</th>
                            <th style="padding: 12px; text-align: right; border-bottom: 2px solid #e0e0e0;">Entropy</th>
                            <th style="padding: 12px; text-align: right; border-bottom: 2px solid #e0e0e0;">Genes</th>
                        </tr>
                    </thead>
                    <tbody>
                        {''.join(sample_rows)}
                    </tbody>
                </table>
            </div>

            <div class="section json-section">
                <details>
                    <summary>View Raw Data (JSON)</summary>
                    <pre>{json_data}</pre>
                </details>
            </div>
        </div>

        <div class="footer">
            <p><strong>MK4 Biomarker Analysis Engine</strong></p>
            <p>Alexandria Dynamics | Public Version v1.0</p>
            <p style="margin-top: 10px; font-size: 0.9em;">
                Generated with Python / scipy / numpy / matplotlib
            </p>
        </div>
    </div>
</body>
</html>"""

    output_path = output_dir / "index.html"
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"Generated: index.html")
    return output_path

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
            html_path = generate_html_report(output_dir, dataset, chaos_results)
            print(f"\nReport: {html_path.relative_to(PROJECT_ROOT)}")
            print(f"   Open in browser to view results")

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
