# MK4 Disease Biomarker

**Alexandria Dynamics Research Project**

Frequency-domain biomarkers for cancer and autoimmune disease detection from blood RNA-seq data.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

MK4 (Molecular Kinetics 4th-order) uses frequency decomposition of RNA expression data to detect disease signatures.

**Key Results:**
- **78% accuracy** detecting cancer from blood
- **79% accuracy** detecting multiple sclerosis
- **92% accuracy** detecting type 2 diabetes from tissue

## Quick Start

```python
from mk4 import MK4Analyzer

# Load your RNA-seq data
analyzer = MK4Analyzer()
results = analyzer.fit('your_data.csv')

print(f"Accuracy: {results.accuracy:.1%}")
```

## Results Summary

| Disease | Tissue | n | Accuracy | Sensitivity | p-value |
|---------|--------|---|----------|-------------|---------|
| Cancer  | Blood  | 542 | 78% | 83% | <0.001 |
| MS      | Blood  | 29  | 79% | 79% | <0.001 |
| T2D     | Islets | 13  | 92% | 100% | 0.005 |

## Installation

```bash
git clone https://github.com/Alexandria-dynamics/mk4-biomarker
cd mk4-biomarker
pip install -r requirements.txt
pip install -e .
```

## Documentation

Coming soon:
- Methodology
- Tutorial
- API Reference

## Publication

**Status:** In preparation

**Code by:** Architekt (Alexandria Dynamics)
**Paper by:** Tomas Vavra

## Project Status

Under Development

Core implementation in progress. See [Issues](../../issues) for roadmap.

## Contact

- **GitHub:** [@Alexandria-dynamics](https://github.com/Alexandria-dynamics)
- **Issues:** [Report bugs or request features](../../issues)

## License

MIT License - see [LICENSE](LICENSE) file

## Citation

```bibtex
@software{architekt2026mk4,
  author = {Architekt},
  title = {MK4 Disease Biomarker},
  year = {2026},
  publisher = {Alexandria Dynamics},
  url = {https://github.com/Alexandria-dynamics/mk4-biomarker}
}
```

---

**Alexandria Dynamics** | Pattern Recognition Across Scales
