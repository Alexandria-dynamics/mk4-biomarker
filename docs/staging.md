# Disease Staging Model — MK4 v1.0

## Concept

MK4 is not designed as a binary classifier (disease / no disease). The biological reality is a continuous spectrum:

    Healthy → Destabilizing → Pre-clinical → Manifest

Our goal is to detect **where on this spectrum** a biological system currently is.

## Three-Stage Framework

### Stage 0: Baseline (Healthy)
- Normal transcriptional dynamics
- Chaos metric at population baseline
- No regulatory instability detected

### Stage 1: Early Destabilization ("On the path")
- Subtle shift in transcriptional chaos
- System beginning to lose regulatory stability
- No clinical symptoms yet
- Pre-symptomatic detection opportunity

### Stage 2: Pre-Clinical ("Near the edge")
- Significant chaos elevation
- Measurable regulatory dysfunction
- Early biomarker signals present
- High-risk state — intervention window

### Stage 3: Manifest Disease ("Active")
- Clear transcriptional dysregulation
- Clinical diagnosis criteria met
- Active disease process
- Treatment required

## Why Staging Matters

A binary test answers: "Do you have cancer?"

A staging test answers: "Where is your biological system on the path toward disease?"

The second question is clinically more valuable:

- Stage 1 detection → lifestyle intervention
- Stage 2 detection → preventive treatment
- Stage 3 detection → active treatment
- Early staging → better outcomes

## Progression Patterns

### RNA-level diseases (Cancer, MS, T2D)

Chaos INCREASES with progression:

```
    Stage 0 → Stage 1 → Stage 2 → Stage 3
    Low chaos ─────────────────→ High chaos
```

Mechanism: Transcriptional dysregulation increases as disease progresses.

### Protein-level diseases (Parkinson's, Alzheimer's)

Chaos DECREASES with progression:

```
    Stage 0 → Stage 1 → Stage 2 → Stage 3
    High chaos ──────────────→ Low chaos
```

Mechanism: Cell death reduces transcriptional activity. Different pattern = different biology.

## Available Progression Data

### Type 2 Diabetes (GSE38642)

Clinical progression marker: **HbA1c**

| HbA1c | Clinical Status | MK4 Stage |
|-------|----------------|-----------|
| <5.7% | Normal | Stage 0 |
| 5.7-6.4% | Pre-diabetes | Stage 1 |
| 6.5-7.0% | Early T2D | Stage 2 |
| >7.0% | Manifest T2D | Stage 3 |

- HbA1c values available in dataset
- Can directly map chaos to clinical progression
- Priority analysis for v1.1

### Parkinson's Disease (GSE49036)

Progression marker: **Braak staging**

| Braak Stage | Clinical Status | MK4 Stage |
|-------------|----------------|-----------|
| 0 | Control | Stage 0 |
| 1-2 | Preclinical | Stage 1 |
| 3-4 | Symptomatic | Stage 2 |
| 5-6 | Advanced | Stage 3 |

- Braak stages available in dataset
- Negative control — chaos should decrease
- Validates two-pattern model

## Planned Analysis

1. Parse HbA1c from GSE38642
2. Correlate chaos with HbA1c levels
3. Define stage thresholds from data
4. Validate on GSE76894, GSE25724
5. Same for Braak stages in GSE49036
6. Publish staging model

## Why This Changes the Paper

From:
> "MK4 detects disease"

To:
> "MK4 measures biological system stability across the disease progression spectrum"

This is:
- More clinically useful
- More scientifically honest
- More defensible
- More novel

No existing blood test does this.

## The Key Insight

> Image AI can only see what's already there.
> MK4 can potentially detect the system BEFORE the disease becomes visible.

This is not detection. This is **progression measurement**.

---

*Alexandria Dynamics — v1.0.0 — February 2026*
