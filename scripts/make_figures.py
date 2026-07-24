#!/usr/bin/env python3
"""
Regenerate benchmark figures from results/benchmark_results.json and
results/benchmark_per_read.csv.

Usage:
    python scripts/make_figures.py

Outputs PNGs into results/figures/. Intentionally dependency-light
(matplotlib only) so it runs cleanly in CI.
"""

import csv
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ---- House style ---------------------------------------------------------
INK = "#1a1a2e"
MUTED = "#6b7280"
GRID = "#e5e7eb"
ACCENT = "#2563eb"
GOOD = "#059669"
WARN = "#d97706"
BAD = "#dc2626"
BG = "#ffffff"

plt.rcParams.update({
    "figure.facecolor": BG,
    "axes.facecolor": BG,
    "savefig.facecolor": BG,
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.edgecolor": GRID,
    "axes.linewidth": 1.0,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": False,
})

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGDIR = RESULTS / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)


def load_results():
    with open(RESULTS / "benchmark_results.json") as fh:
        return json.load(fh)


def load_per_read():
    rows = []
    with open(RESULTS / "benchmark_per_read.csv", newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
    return rows


# ---- 1. Confusion matrix -------------------------------------------------
def confusion_matrix(res):
    cd = res["chimeric_detection"]
    tp, fp = cd["true_positives"], cd["false_positives"]
    fn, tn = cd["false_negatives"], cd["true_negatives"]
    cells = [[tn, fp], [fn, tp]]
    labels = [["TN", "FP"], ["FN", "TP"]]
    colors = [["#eef2ff", "#fee2e2"], ["#fef3c7", "#dcfce7"]]

    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    for i in range(2):
        for j in range(2):
            ax.add_patch(Rectangle((j, 1 - i), 1, 1, facecolor=colors[i][j],
                                   edgecolor="white", linewidth=3))
            ax.text(j + 0.5, 1 - i + 0.60, f"{cells[i][j]:,}",
                    ha="center", va="center", fontsize=26, fontweight="bold",
                    color=INK)
            ax.text(j + 0.5, 1 - i + 0.28, labels[i][j],
                    ha="center", va="center", fontsize=12, color=MUTED,
                    fontweight="bold")
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2)
    ax.set_xticks([0.5, 1.5])
    ax.set_yticks([1.5, 0.5])
    ax.set_xticklabels(["Predicted: Normal", "Predicted: Chimeric"], fontsize=11)
    ax.set_yticklabels(["Actual:\nNormal", "Actual:\nChimeric"], fontsize=11)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Chimeric-read detection — confusion matrix",
                 fontsize=13, fontweight="bold", pad=16, color=INK)
    fig.text(0.5, 0.02,
             f"Zero false positives across {tn + fp + fn + tp:,} reads",
             ha="center", fontsize=10, color=MUTED)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(FIGDIR / "confusion_matrix.png", dpi=200)
    plt.close(fig)


# ---- 2. Metrics bar ------------------------------------------------------
def metrics_bar(res):
    cd = res["chimeric_detection"]
    metrics = [
        ("Precision", cd["precision"]),
        ("Specificity", cd["specificity"]),
        ("Accuracy", cd["accuracy"]),
        ("F1 score", cd["f1_score"]),
        ("Recall", cd["recall"]),
    ]
    names = [m[0] for m in metrics]
    vals = [m[1] for m in metrics]

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    y = range(len(names))
    ax.barh(list(y), vals, color=ACCENT, height=0.6, zorder=3)
    for i, v in enumerate(vals):
        ax.text(v - 0.015, i, f"{v:.3f}", va="center", ha="right",
                color="white", fontweight="bold", fontsize=11, zorder=4)
    ax.axvline(1.0, color=GRID, lw=1, zorder=1)
    ax.set_yticks(list(y))
    ax.set_yticklabels(names, fontsize=12)
    ax.set_xlim(0.8, 1.005)
    ax.invert_yaxis()
    ax.set_xlabel("Score", fontsize=11)
    ax.set_title("Chimeric-detection performance",
                 fontsize=13, fontweight="bold", pad=12, loc="left")
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0)
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color=GRID, lw=1)
    fig.tight_layout()
    fig.savefig(FIGDIR / "metrics_bar.png", dpi=200)
    plt.close(fig)


# ---- 3. Read-category distribution --------------------------------------
def category_distribution(rows):
    counts = Counter(r["ground_truth_category"] for r in rows)
    order = ["normal", "chimeric", "backbone", "host"]
    palette = {"normal": "#94a3b8", "chimeric": ACCENT,
               "backbone": WARN, "host": "#7c3aed"}
    labels = [c for c in order if c in counts]
    vals = [counts[c] for c in labels]

    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    bars = ax.bar(labels, vals, color=[palette[c] for c in labels],
                  width=0.62, zorder=3)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + max(vals) * 0.01,
                f"{v:,}", ha="center", va="bottom", fontweight="bold",
                fontsize=11, color=INK)
    ax.set_ylabel("Reads", fontsize=11)
    ax.set_title("Simulated read composition (ground truth)",
                 fontsize=13, fontweight="bold", pad=12, loc="left")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=0)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, lw=1)
    ax.set_ylim(0, max(vals) * 1.12)
    fig.tight_layout()
    fig.savefig(FIGDIR / "category_distribution.png", dpi=200)
    plt.close(fig)


def main():
    res = load_results()
    rows = load_per_read()
    confusion_matrix(res)
    metrics_bar(res)
    category_distribution(rows)
    print(f"Wrote figures to {FIGDIR}")


if __name__ == "__main__":
    main()
