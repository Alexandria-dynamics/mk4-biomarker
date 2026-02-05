# Quick Start Guide

## Installation

```bash
git clone https://github.com/alexandria-dynamics/mk4-biomarker.git
cd mk4-biomarker
pip install -r requirements.txt
```

## Usage

### Step 1: Get Data
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25724/soft/GSE25724_family.soft.gz
mv GSE25724_family.soft.gz data/input/
```

### Step 2: Run
```bash
python scripts/parse_geo_soft.py
```

### Step 3: Results
```bash
cd results
ls -lt | head -5
cat 20260205_143022_GSE25724/json/chaos_results.json
```

## Example Output

```json
{
  "statistics": {
    "control_mean": 14.367534,
    "disease_mean": 14.380579,
    "difference": 0.013045,
    "percent_change": 0.09,
    "p_value": 0.0037,
    "significant": true
  }
}
```

## Troubleshooting

**No files found:** Place `.gz` files in `data/input/`
**Parser error:** Verify valid GEO SOFT format
**No plots:** Install matplotlib
