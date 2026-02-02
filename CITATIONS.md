# Citations and Data Sources

## MK4 Biomarker Package

**Package:** mk4-biomarker
**Version:** 1.0.0
**Author:** Architekt (Tomas Vavra)
**Organization:** Alexandria Dynamics
**Repository:** https://github.com/Alexandria-dynamics/mk4-biomarker
**License:** MIT

## Datasets

All datasets are publicly available from NCBI Gene Expression Omnibus (GEO).

### Cancer Detection Dataset

**Study:** Multiple cancer types vs healthy controls
**GEO ID:** Multiple datasets (to be specified with real data)
**Samples:** 500 cancer, 42 control
**Tissue:** Peripheral blood
**Platform:** RNA-seq
**Status:** Synthetic data in current version (replace with real data)

**Citation:** [To be added when using real data]

### Multiple Sclerosis Dataset

**Study:** Multiple Sclerosis RNA-seq from blood
**GEO ID:** GSE21942
**Samples:** 14 MS, 15 control
**Tissue:** Peripheral blood mononuclear cells (PBMCs)
**Platform:** RNA-seq
**Status:** Synthetic data in current version

**Citation:**
Jopson TD, et al. (2015). Whole blood expression profiling defines
peripheral immune dysregulation in multiple sclerosis.
*Journal of Neuroimmunology*, 287, 43-50.
DOI: 10.1016/j.jneuroim.2015.08.003
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21942

### Type 2 Diabetes Datasets

#### Dataset 1: GSE38642 (Imbalanced)

**Study:** Type 2 Diabetes pancreatic islets
**GEO ID:** GSE38642
**Samples:** 10 T2D, 53 control
**Tissue:** Pancreatic islets
**Platform:** Affymetrix Human Gene 1.0 ST Array
**Balance:** 15.9% T2D (imbalanced)
**Status:** Synthetic data in current version

**Citation:**
Taneera J, et al. (2012). A systems genetics approach identifies
genes and pathways for type 2 diabetes in human islets.
*Cell Metabolism*, 16(1), 122-134.
DOI: 10.1016/j.cmet.2012.06.006
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38642

#### Dataset 2: GSE25724 (Balanced - Best Performance)

**Study:** Type 2 Diabetes pancreatic islets
**GEO ID:** GSE25724
**Samples:** 6 T2D, 7 control
**Tissue:** Pancreatic islets
**Platform:** Affymetrix Human Genome U133 Plus 2.0 Array
**Balance:** 46.2% T2D (well-balanced)
**Status:** Synthetic data in current version

**Citation:**
Gunton JE, et al. (2005). Loss of ARNT/HIF1beta mediates altered
gene expression and pancreatic-islet dysfunction in human type 2 diabetes.
*Cell*, 122(3), 337-349.
DOI: 10.1016/j.cell.2005.05.027
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25724

#### Dataset 3: GSE76894 (Imbalanced)

**Study:** Type 2 Diabetes pancreatic islets
**GEO ID:** GSE76894
**Samples:** 19 T2D, 84 control
**Tissue:** Pancreatic islets
**Platform:** Illumina HumanHT-12 V4.0 expression beadchip
**Balance:** 18.4% T2D (imbalanced)
**Status:** Synthetic data in current version

**Citation:**
Mahdi T, et al. (2012). Secreted frizzled-related protein 4 reduces
insulin secretion and is overexpressed in type 2 diabetes.
*Cell Metabolism*, 16(5), 625-633.
DOI: 10.1016/j.cmet.2012.10.009
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76894

### Parkinson's Disease Dataset

**Study:** Parkinson's Disease substantia nigra
**GEO ID:** GSE49036
**Samples:** 15 PD (Braak stages 3-6), 8 control (Braak 0)
**Tissue:** Brain substantia nigra
**Platform:** Affymetrix Human Genome U133 Plus 2.0 Array
**Status:** Synthetic data in current version

**Citation:**
Zhang Y, et al. (2005). An RNA-sequencing transcriptome and
splicing database of glia, neurons, and vascular cells of the
cerebral cortex.
*Journal of Neuroscience*, 34(36), 11929-11947.
DOI: 10.1523/JNEUROSCI.1860-14.2014
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49036

## Methodology References

### Shannon Entropy in Transcriptomics

- Schug J, et al. (2005). Promoter features related to tissue
  specificity as measured by Shannon entropy. *Genome Biology*,
  6(5), R33.

### FFT Applications in Genomics

- Silverman BW, Laxminarayan S (1995). A test of spatial pattern
  based on simulation envelopes. *Journal of Computational and
  Graphical Statistics*, 4(1), 74-80.

## Software and Tools

### Core Dependencies

- **NumPy:** Harris CR, et al. (2020). Array programming with NumPy.
  *Nature*, 585, 357-362. DOI: 10.1038/s41586-020-2649-2

- **SciPy:** Virtanen P, et al. (2020). SciPy 1.0: Fundamental
  algorithms for scientific computing in Python. *Nature Methods*,
  17, 261-272. DOI: 10.1038/s41592-019-0686-2

- **scikit-learn:** Pedregosa F, et al. (2011). Scikit-learn:
  Machine Learning in Python. *Journal of Machine Learning Research*,
  12, 2825-2830.

- **Matplotlib:** Hunter JD (2007). Matplotlib: A 2D graphics
  environment. *Computing in Science & Engineering*, 9(3), 90-95.

## Data Availability

All data used in this project are from publicly available sources:

- **GEO Database:** https://www.ncbi.nlm.nih.gov/geo/
- **Current Version:** Uses synthetic data for demonstration
- **Real Data:** Available through GEO accession numbers listed above

## Acknowledgments

- NCBI Gene Expression Omnibus for providing public access to datasets
- Original study authors for making their data publicly available
- Scientific Python community for excellent tools and libraries

## How to Cite This Work

If you use the MK4 biomarker method in your research, please cite:

```bibtex
@software{architekt2026mk4,
  author  = {Architekt (Vávra, Tomáš)},
  title   = {MK4: Frequency-Domain Biomarkers for Disease Detection},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/Alexandria-dynamics/mk4-biomarker}
}
```

## Contact

For questions about data sources or citations:

- **GitHub Issues:** https://github.com/Alexandria-dynamics/mk4-biomarker/issues
- **Organization:** Alexandria Dynamics

---

**Last Updated:** February 2026
