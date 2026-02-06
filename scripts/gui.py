import json
import os
from datetime import datetime
import matplotlib.pyplot as plt

def parse_file(filepath):
    """
    Unified entrypoint expected by gui.py.
    Returns dataset dict or None if unknown format.
    """
    ftype = detect_file_type(filepath)
    if ftype == "soft":
        return parse_soft_file(filepath)
    if ftype == "matrix":
        return parse_matrix_file(filepath)
    return None


def compute_chaos_metrics(dataset):
    """
    Returns dict expected by gui.py:
    {
      "entropies": {"CONTROL":[...], "DISEASE":[...]},
      "statistics": {...}
    }
    """
    entropies = {"CONTROL": [], "DISEASE": []}

    for sample in dataset["samples"].values():
        lbl = sample.get("label")
        if lbl not in ("CONTROL", "DISEASE"):
            continue
        e = compute_sample_entropy(sample["expression"])
        if not np.isnan(e):
            entropies[lbl].append(float(e))

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
            "t_stat": float(t),
            "n_control": int(len(c)),
            "n_disease": int(len(d)),
        }
    else:
        stats = {
            "error": "insufficient samples for statistics",
            "n_control": int(len(c)),
            "n_disease": int(len(d)),
        }

    return {"entropies": entropies, "statistics": stats}


def create_output_directory(dataset_name):
    """
    Creates unique output folder under RESULTS_DIR and returns Path.
    """
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    out = RESULTS_DIR / f"run_{ts}_{dataset_name}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    return out


def save_results(output_dir, dataset, chaos_results):
    """
    Saves metadata + minimal JSON. Also creates index.html gallery placeholder.
    """
    output_dir = Path(output_dir)

    meta = {
        "dataset_name": dataset.get("dataset_name"),
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "n_samples": len(dataset.get("samples", {})),
        "n_control": sum(1 for s in dataset["samples"].values() if s.get("label") == "CONTROL"),
        "n_disease": sum(1 for s in dataset["samples"].values() if s.get("label") == "DISEASE"),
    }

    (output_dir / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    (output_dir / "chaos_results.json").write_text(json.dumps(chaos_results, indent=2), encoding="utf-8")

    # Create/refresh simple index.html listing figures
    write_index_html(output_dir)


def generate_plot(output_dir, chaos_results):
    """
    Generates a simple boxplot + histogram into figures/, then refreshes index.html.
    """
    output_dir = Path(output_dir)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    c = chaos_results["entropies"]["CONTROL"]
    d = chaos_results["entropies"]["DISEASE"]

    # Boxplot
    plt.figure()
    plt.boxplot([c, d], labels=["CONTROL", "DISEASE"])
    plt.title("Chaos (Entropy) by Group")
    plt.ylabel("Entropy")
    box_path = fig_dir / "01_chaos_boxplot.png"
    plt.savefig(box_path, dpi=160, bbox_inches="tight")
    plt.close()

    # Histogram
    plt.figure()
    if len(c) > 0:
        plt.hist(c, bins=20, alpha=0.6, label="CONTROL")
    if len(d) > 0:
        plt.hist(d, bins=20, alpha=0.6, label="DISEASE")
    plt.title("Chaos (Entropy) Distribution")
    plt.xlabel("Entropy")
    plt.ylabel("Count")
    plt.legend()
    hist_path = fig_dir / "02_chaos_hist.png"
    plt.savefig(hist_path, dpi=160, bbox_inches="tight")
    plt.close()

    write_index_html(output_dir)


def write_index_html(output_dir):
    """
    Creates single-page HTML gallery of all images in output_dir/figures.
    """
    output_dir = Path(output_dir)
    fig_dir = output_dir / "figures"
    imgs = []
    if fig_dir.exists():
        imgs = sorted([p for p in fig_dir.iterdir() if p.suffix.lower() in (".png", ".jpg", ".jpeg", ".webp")])

    cards = []
    for p in imgs:
        rel = f"figures/{p.name}"
        cards.append(f"""
        <figure style="border:1px solid #ddd;border-radius:10px;padding:10px;margin:0;background:#fff">
          <a href="{rel}" target="_blank" rel="noopener">
            <img src="{rel}" style="width:100%;height:auto;border-radius:8px" loading="lazy">
          </a>
          <figcaption style="font-size:12px;color:#555;margin-top:6px;word-break:break-word">{p.name}</figcaption>
        </figure>
        """)

    grid = "\n".join(cards) if cards else "<p>No figures yet.</p>"

    html_doc = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>MK4 Output</title>
</head>
<body style="font-family:Arial, sans-serif; margin:24px; background:#f5f7fa;">
  <h1 style="margin:0 0 8px 0;">MK4 Output</h1>
  <p style="color:#666;margin:0 0 18px 0;">Analytical output (non-diagnostic). Folder: {output_dir.name}</p>
  <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px;">
    {grid}
  </div>
</body>
</html>
"""
    (output_dir / "index.html").write_text(html_doc, encoding="utf-8")
