# MK4 Development Roadmap

## v1.0.0 — Current Release (February 2026)

**Status:** Released

### What We Have

- Complete core framework (engine, analysis, visualization, utils)
- Four experiment scripts (cancer, MS, diabetes, Parkinson's)
- Preliminary results across 8 datasets
- Two key discoveries:
  - Dataset balance effect
  - RNA-level vs protein-level scope definition
- Open-source, reproducible pipeline
- Professional documentation and citations

### Preliminary Performance

| Disease | Accuracy | False Positive Rate | Context |
|---------|----------|-------------------|---------|
| Cancer | 78% | ~22% | Mammography: ~25% FP |
| MS | 79% | ~21% | First blood-based test |
| T2D (balanced) | 92% | ~8% | With balanced data |
| Parkinson's | ~61% | ns | Expected: negative control |

### Known Limitations (v1.0)

- Results are preliminary, require independent validation
- Current pipeline uses synthetic data (real data integration pending)
- False positive rate needs improvement for clinical use
- Cost of RNA-seq remains a barrier for screening applications
- Small sample sizes in some datasets

---

## v1.1 — Specificity & False Detection Reduction

**Status:** Planned
**Priority:** HIGH — Top priority

### Goals

**Primary goal: Reduce false detections.**

False positives are the main barrier to clinical adoption.
Our v1.0 false positive rate (~22% for cancer) is comparable
to mammography (~25%), but we aim to significantly improve this.

### Planned Work

#### 1. Threshold Optimization
- Cross-validated threshold selection
- Disease-specific optimal thresholds
- Balance sensitivity vs specificity explicitly
- ROC-based decision point selection

#### 2. Multi-Metric Classification
- Use chaos + coherence + band power together
- Ensemble approach (not single metric)
- Machine learning classifier on MK4 features
- Reduces reliance on single threshold

#### 3. Confidence Scoring
- Add confidence scores to each prediction
- "High confidence positive" vs "uncertain"
- Allows two-tier reporting
- Clinicians can act on high-confidence results

#### 4. False Positive Analysis
- Characterize false positive profiles
- Identify what makes FP cases different
- Build FP filter
- Systematic reduction

### Target Performance (v1.1)

| Disease | v1.0 FP Rate | v1.1 Target FP Rate |
|---------|-------------|-------------------|
| Cancer | ~22% | <15% |
| MS | ~21% | <15% |
| T2D | ~8% | <5% |

---

## v1.2 — Minimal Gene Panel

**Status:** Planned
**Priority:** MEDIUM

### Goals

Identify the minimum number of genes needed for reliable
chaos calculation. This reduces cost and enables targeted
sequencing panels.

### Planned Work

#### 1. Gene Importance Analysis
- Identify which genes contribute most to chaos signature
- Feature importance ranking
- Stability analysis across datasets

#### 2. Panel Size Optimization
- Test chaos reliability with decreasing gene counts
- Find minimum panel size (target: 500-2000 genes)
- Validate panel performance vs full transcriptome

#### 3. Cost Implications
- Full RNA-seq: $1000-3000
- Targeted panel (500 genes): $50-200
- Document cost reduction pathway

### Target

- Identify minimal gene panel (<2000 genes)
- Validate chaos metric on panel data
- Document cost reduction potential

---

## v2.0 — Independent Validation

**Status:** Planned
**Priority:** HIGH (for publication)

### Goals

Validate all preliminary results on independent datasets.
This is required for any scientific publication.

### Planned Work

#### 1. Real Data Integration
- Download and process real GEO datasets
- Run full pipeline on real data
- Compare with preliminary results
- Document any discrepancies

#### 2. Independent Test Sets
- Split datasets into train/test
- Cross-validation on real data
- Report validated performance metrics
- Only claim what is validated

#### 3. Replication
- Find independent datasets for same diseases
- Test reproducibility across studies
- Assess generalizability

### Publication Criteria

Before publishing, we require:
- [ ] Real data validation for all claimed results
- [ ] Independent test set performance reported
- [ ] All code and data reproducible
- [ ] Honest reporting of limitations
- [ ] Zero false claims

---

## v3.0 — Clinical Feasibility

**Status:** Future
**Priority:** LONG-TERM

### Goals

Assess feasibility for clinical applications.

### Planned Work

- Clinical use case analysis
- Regulatory pathway assessment
- Two-stage diagnostic framework design
- Cost-effectiveness analysis
- Prospective study design

### Note

Clinical application requires extensive validation
beyond what is possible in v1.0. This is a long-term goal.

---

## Philosophy

> We build MK4 with scientific integrity as the foundation.
> Version 1.0 is a proof of concept. We will not claim more
> than what our data supports. Each version will be better,
> more validated, and more honest than the last.
>
> False detections are our enemy. Reducing them is our
> primary mission going forward.
>
> — Alexandria Dynamics

---

**Last Updated:** February 2026
**Repository:** https://github.com/Alexandria-dynamics/mk4-biomarker
