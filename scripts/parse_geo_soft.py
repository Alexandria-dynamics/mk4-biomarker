#!/usr/bin/env python3
"""
MK4 Biomarker - GEO Parser (OPRAVENÝ)
Alexandria Dynamics

Podporuje DVA formáty:
1. SOFT formát - klasické .soft nebo .soft.gz soubory s ^SAMPLE sekcemi
2. MATRIX formát - tab-delimited data matice (.txt nebo .txt.gz)

Automaticky detekuje formát a použije správný parser.
"""

import base64
import gzip
import json
from io import BytesIO
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
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
# FILE TYPE DETECTION
# ═══════════════════════════════════════════════════════════════

def detect_file_type(filepath):
    """
    Detekuje typ souboru (SOFT nebo MATRIX).
    
    Returns:
        'soft' - klasický SOFT formát s ^SAMPLE sekcemi
        'matrix' - tab-delimited data matice
        None - neznámý formát
    """
    filepath = Path(filepath)
    
    # Otevři soubor (gzip nebo plain)
    if filepath.suffix == '.gz':
        f = gzip.open(filepath, 'rt')
    else:
        f = open(filepath, 'r')
    
    try:
        # Přečti prvních 50 řádků
        lines = [f.readline() for _ in range(50)]
        
        # Detekce SOFT: obsahuje ^SAMPLE nebo ^SERIES
        has_soft_markers = any(
            line.startswith('^SAMPLE') or line.startswith('^SERIES')
            for line in lines
        )
        
        # Detekce MATRIX: první řádek je header s tab-separated sample names
        # druhý řádek začíná gene ID (např. ENSG)
        first_line = lines[0].strip()
        second_line = lines[1].strip() if len(lines) > 1 else ""
        
        has_tabs = '\t' in first_line
        starts_with_gene = second_line.startswith('ENSG') or second_line.startswith('ILMN')
        
        if has_soft_markers:
            return 'soft'
        elif has_tabs and starts_with_gene:
            return 'matrix'
        else:
            return None
            
    finally:
        f.close()

# ═══════════════════════════════════════════════════════════════
# SOFT PARSER (původní)
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
    print(f"📂 Parsing SOFT: {filepath.name}")
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
    
    print(f"✅ Found {len(valid_samples)} valid samples")
    print(f"   Controls: {sum(1 for v in valid_samples.values() if v['label']=='CONTROL')}")
    print(f"   Disease:  {sum(1 for v in valid_samples.values() if v['label']=='DISEASE')}")
    
    return {
        'dataset_name': dataset_name,
        'samples': valid_samples
    }

# ═══════════════════════════════════════════════════════════════
# MATRIX PARSER (nový!)
# ═══════════════════════════════════════════════════════════════

def classify_sample_name(sample_name):
    """
    Klasifikuj vzorek jako CONTROL nebo DISEASE na základě názvu.
    
    Podporované pattern:
    - HD-, Control-* → CONTROL
    - Vše ostatní → DISEASE
    """
    s = sample_name.upper()
    
    # Kontroly
    if 'HD-' in s or 'CONTROL' in s or 'NORMAL' in s or 'HEALTHY' in s:
        return 'CONTROL'
    
    # Disease keywords
    if any(kw in s for kw in ['CANCER', 'TUMOR', 'DISEASE', 'PATIENT', 
                               'BREAST', 'LUNG', 'GBM', 'CRC', 'PANCR', 'LIVER']):
        return 'DISEASE'
    
    # Default: pokud není jasné, zkus heuristiku
    # Vzorky s čísly nebo kódy obvykle = disease
    return 'DISEASE'

def parse_matrix_file(filepath):
    """
    Parse tab-delimited matrix file (gene × sample).
    
    Formát:
    - První řádek: sample names (tab-separated)
    - Zbylé řádky: gene_id \t value1 \t value2 ...
    
    Returns:
        dict: stejná struktura jako parse_soft_file
    """
    filepath = Path(filepath)
    dataset_name = filepath.stem.replace('_data_matrix', '').replace('.txt', '').replace('.gz', '')
    
    print(f"\n{'='*70}")
    print(f"📂 Parsing MATRIX: {filepath.name}")
    print(f"   (může trvat ~10-20 sekund pro velké soubory...)")
    print(f"{'='*70}")
    
    # Načti matici pomocí pandas
    if filepath.suffix == '.gz':
        df = pd.read_csv(filepath, sep='\t', index_col=0, compression='gzip')
    else:
        df = pd.read_csv(filepath, sep='\t', index_col=0)
    
    print(f"✅ Loaded: {df.shape[0]} genes × {df.shape[1]} samples")
    
    # Klasifikuj samples
    samples = {}
    
    for col in df.columns:
        label = classify_sample_name(col)
        
        # Extrahuj expression values (filtruj >0)
        expression = {}
        for gene_id, value in df[col].items():
            if value > 0:  # Ignore zero/negative
                expression[str(gene_id)] = float(value)
        
        samples[col] = {
            'label': label,
            'expression': expression
        }
    
    # Statistika
    controls = sum(1 for s in samples.values() if s['label'] == 'CONTROL')
    diseases = sum(1 for s in samples.values() if s['label'] == 'DISEASE')
    
    print(f"✅ Classified samples:")
    print(f"   Controls: {controls}")
    print(f"   Disease:  {diseases}")
    
    return {
        'dataset_name': dataset_name,
        'samples': samples
    }

# ═══════════════════════════════════════════════════════════════
# UNIFIED PARSER (automatická detekce)
# ═══════════════════════════════════════════════════════════════

def parse_file(filepath):
    """
    Automaticky detekuje formát a použije správný parser.
    
    Returns:
        dict nebo None (pokud nepodporovaný formát)
    """
    file_type = detect_file_type(filepath)
    
    if file_type == 'soft':
        return parse_soft_file(filepath)
    elif file_type == 'matrix':
        return parse_matrix_file(filepath)
    else:
        print(f"⚠️  Neznámý formát souboru: {filepath.name}")
        return None

# ═══════════════════════════════════════════════════════════════
# CHAOS ANALYSIS (beze změny)
# ═══════════════════════════════════════════════════════════════

def compute_chaos_metrics(dataset):
    """Compute Shannon entropy (chaos metric) for each sample."""
    samples = dataset['samples']
    
    # Get common genes
    all_genes = [set(s['expression'].keys()) for s in samples.values()]
    common_genes = sorted(set.intersection(*all_genes))
    
    print(f"\n🧬 Common genes: {len(common_genes)}")
    
    control_entropy = []
    disease_entropy = []
    sample_details = {}
    
    for sample_id, data in samples.items():
        label = data['label']
        
        # Get expression values
        values = np.array([data['expression'][g] for g in common_genes])
        
        # Filter positive values only
        values = values[values > 0]
        
        if len(values) == 0:
            continue
        
        # Normalize
        values_norm = values / values.sum()
        
        # Shannon entropy
        ent = entropy(values_norm, base=2)
        
        sample_details[sample_id] = {
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
# OUTPUT (stejné jako předtím)
# ═══════════════════════════════════════════════════════════════

def create_output_directory(dataset_name):
    """Create timestamped output directory."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = RESULTS_DIR / f"{timestamp}_{dataset_name}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n📁 Output: {output_dir.relative_to(PROJECT_ROOT)}")
    return output_dir

def save_results(output_dir, dataset, chaos_results):
    """Save results to JSON."""
    
    # Chaos results
    with open(output_dir / "chaos_results.json", 'w') as f:
        json.dump(chaos_results, f, indent=2)
    
    # Metadata
    metadata = {
        'dataset_name': dataset['dataset_name'],
        'timestamp': datetime.now().isoformat(),
        'num_samples': len(dataset['samples']),
        'num_controls': sum(1 for s in dataset['samples'].values() if s['label'] == 'CONTROL'),
        'num_disease': sum(1 for s in dataset['samples'].values() if s['label'] == 'DISEASE'),
    }
    
    with open(output_dir / "metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"✅ Saved JSON files")

def generate_plot(output_dir, chaos_results):
    """Generate visualization."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        control = chaos_results['control_entropy']
        disease = chaos_results['disease_entropy']
        stats = chaos_results['statistics']
        
        if 'error' in stats:
            print("⚠️  Skipping plot - insufficient data")
            return
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        bp = ax.boxplot([control, disease], 
                         labels=['Control', 'Disease'],
                         patch_artist=True)
        
        bp['boxes'][0].set_facecolor('#2ecc71')
        bp['boxes'][1].set_facecolor('#e74c3c')
        
        for box in bp['boxes']:
            box.set_alpha(0.7)
        
        ax.set_ylabel('Shannon Entropy (bits)', fontsize=12)
        ax.set_title('Transcriptional Chaos: Control vs Disease', 
                     fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # p-value
        sig_text = "***" if stats['p_value'] < 0.001 else \
                   "**" if stats['p_value'] < 0.01 else \
                   "*" if stats['p_value'] < 0.05 else "ns"
        
        ax.text(0.5, 0.98, f"p = {stats['p_value']:.2e} {sig_text}", 
                transform=ax.transAxes, ha='center', va='top', fontsize=11)
        
        plt.tight_layout()
        plt.savefig(output_dir / "chaos_plot.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✅ Saved plot")
        
    except ImportError:
        print("⚠️  matplotlib not installed - skipping plot")

# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    """Main pipeline."""
    
    print("\n" + "="*70)
    print("  MK4 BIOMARKER - RNA Chaos Analysis (OPRAVENÝ)")
    print("  Alexandria Dynamics")
    print("  Podporuje: SOFT formát i MATRIX formát")
    print("="*70)
    
    # Find all potential files
    all_files = list(INPUT_DIR.glob("*.gz")) + list(INPUT_DIR.glob("*.soft")) + list(INPUT_DIR.glob("*.txt"))
    
    if not all_files:
        print(f"\n❌ No files in: {INPUT_DIR.relative_to(PROJECT_ROOT)}")
        print(f"\n💡 Place files in: data/input/")
        print(f"   Supported: .soft.gz, .soft, .txt.gz, .txt (matrix)")
        return
    
    print(f"\n✅ Found {len(all_files)} file(s):")
    for f in all_files:
        file_type = detect_file_type(f)
        type_str = f"[{file_type.upper()}]" if file_type else "[UNKNOWN]"
        print(f"   {type_str} {f.name}")
    
    # Process each file
    processed = 0
    
    for filepath in all_files:
        try:
            # Parse (auto-detect format)
            dataset = parse_file(filepath)
            
            if dataset is None:
                continue
            
            if len(dataset['samples']) == 0:
                print(f"⚠️  No valid samples found - skipping")
                continue
            
            # Analyze
            chaos_results = compute_chaos_metrics(dataset)
            
            # Save
            output_dir = create_output_directory(dataset['dataset_name'])
            save_results(output_dir, dataset, chaos_results)
            generate_plot(output_dir, chaos_results)
            
            # Summary
            if 'error' not in chaos_results['statistics']:
                stats = chaos_results['statistics']
                print(f"\n📊 SUMMARY")
                print(f"{'='*70}")
                print(f"  Control: {stats['control_mean']:.6f} ± {stats['control_std']:.6f}")
                print(f"  Disease: {stats['disease_mean']:.6f} ± {stats['disease_std']:.6f}")
                print(f"  Diff:    {stats['difference']:+.6f} ({stats['percent_change']:+.2f}%)")
                print(f"  p-value: {stats['p_value']:.2e}")
                print(f"  Signif:  {'✓ YES' if stats['significant'] else 'NO'}")
                print(f"{'='*70}")
            
            processed += 1
            
        except Exception as e:
            print(f"\n❌ Error processing {filepath.name}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    if processed > 0:
        print(f"\n✅ Complete! Processed {processed} file(s)")
        print(f"📁 Results in: {RESULTS_DIR.relative_to(PROJECT_ROOT)}/")
    else:
        print(f"\n⚠️  No files were successfully processed")
    
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
