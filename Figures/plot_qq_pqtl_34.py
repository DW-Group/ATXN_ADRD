#!/usr/bin/env python3

import os
import math
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import chi2

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "Arial"

# -----------------------------
# paths
# -----------------------------

results_dir = "results"

out_pdf = os.path.join(results_dir, "supp_fig3_pqtl_qqplots_34proteins.pdf")
out_png = os.path.join(results_dir, "supp_fig3_pqtl_qqplots_34proteins.png")

out_source_data = os.path.join(
    results_dir,
    "source_data",
    "source_data_supplementary_fig3.xlsx"
)
os.makedirs(os.path.dirname(out_source_data), exist_ok=True)

markers = [
    ("ANXA5", "ANXA5"),
    ("APOE4", "APOE4"),
    ("BD_pTau_181", "BD-pTau-181"),
    ("CCL3", "CCL3"),
    ("CCL4", "CCL4"),
    ("CHI3L1", "CHI3L1"),
    ("CHIT1", "CHIT1"),
    ("CNTN2", "CNTN2"),
    ("CRH", "CRH"),
    ("CXCL1", "CXCL1"),
    ("FCN2", "FCN2"),
    ("FGF2", "FGF2"),
    ("GDF15", "GDF15"),
    ("IL13", "IL13"),
    ("IL6R", "IL6R"),
    ("KDR", "KDR"),
    ("MSLN", "MSLN"),
    ("NPTX1", "NPTX1"),
    ("NPTX2", "NPTX2"),
    ("NPTXR", "NPTXR"),
    ("NPY", "NPY"),
    ("PDGFRB", "PDGFRB"),
    ("POSTN", "POSTN"),
    ("S100B", "S100B"),
    ("SAA1", "SAA1"),
    ("SFTPD", "SFTPD"),
    ("SLIT2", "SLIT2"),
    ("TAFA5", "TAFA5"),
    ("TEK", "TEK"),
    ("TIMP3", "TIMP3"),
    ("TREM1", "TREM1"),
    ("TREM2", "TREM2"),
    ("VEGFA", "VEGFA"),
    ("VEGFD", "VEGFD"),
]

# -----------------------------
# plotting settings
# -----------------------------

max_points_plot = 40000
tail_points_plot = 3000

# -----------------------------
# helper functions
# -----------------------------

def read_pvalues(marker_id):
    infile = os.path.join(results_dir, f"{marker_id}.mht.txt")

    if not os.path.exists(infile):
        raise FileNotFoundError(f"Input file not found: {infile}")

    df = pd.read_csv(infile, sep="\t", usecols=["Pval"])
    p = pd.to_numeric(df["Pval"], errors="coerce").to_numpy()
    p = p[np.isfinite(p)]
    p = p[(p > 0) & (p <= 1)]

    if len(p) == 0:
        raise ValueError(f"No valid P values found for {marker_id}")

    return p


def calc_lambda_gc(p):
    chisq = chi2.isf(p, df=1)
    return np.nanmedian(chisq) / chi2.ppf(0.5, df=1)


def make_qq_data(p):
    p = np.sort(p)
    n = len(p)

    i = np.arange(1, n + 1)
    expected_p = i / (n + 1.0)

    expected = -np.log10(expected_p)
    observed = -np.log10(p)

    return expected, observed


def thin_points(expected, observed, max_points=40000, tail_points=3000):
    n = len(expected)

    if n <= max_points:
        return expected, observed

    tail_n = min(tail_points, n)
    tail_idx = np.arange(0, tail_n)

    body_n = max_points - tail_n
    body_idx = np.unique(
        np.round(np.linspace(tail_n, n - 1, body_n)).astype(int)
    )

    keep_idx = np.unique(np.concatenate([tail_idx, body_idx]))
    keep_idx.sort()

    return expected[keep_idx], observed[keep_idx]

# -----------------------------
# plot
# -----------------------------

ncol = 5
nrow = 7

fig, axes = plt.subplots(nrow, ncol, figsize=(10, 14))
axes = axes.flatten()

source_records = []

for idx, (marker_id, marker_label) in enumerate(markers):

    print(f"Processing {marker_id}")

    p = read_pvalues(marker_id)
    lambda_gc = calc_lambda_gc(p)

    expected, observed = make_qq_data(p)
    expected_plot, observed_plot = thin_points(
        expected,
        observed,
        max_points=max_points_plot,
        tail_points=tail_points_plot
    )

    max_axis = math.ceil(max(np.nanmax(expected), np.nanmax(observed)))
    max_axis = max(max_axis, 1)

    ax = axes[idx]

    ax.plot(
        [0, max_axis],
        [0, max_axis],
        color="#D62728",
        linewidth=0.45
    )

    ax.scatter(
        expected_plot,
        observed_plot,
        s=1.25,
        color="black",
        alpha=0.70,
        linewidths=0,
        rasterized=True
    )

    ax.text(
        0.05 * max_axis,
        0.90 * max_axis,
        r"$\lambda_{\mathrm{GC}}$ = " + f"{lambda_gc:.2f}",
        ha="left",
        va="center",
        fontsize=7
    )

    ax.set_title(marker_label, fontsize=8.5, fontweight="bold", pad=4)

    ax.set_xlim(0, max_axis)
    ax.set_ylim(0, max_axis)
    ax.set_aspect("equal", adjustable="box")

    ax.tick_params(axis="both", labelsize=6, length=2, width=0.5)

    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    if idx % ncol == 0:
        ax.set_ylabel(r"Observed $-\log_{10}(P)$", fontsize=7)
    else:
        ax.set_ylabel("")

    if idx >= (nrow - 1) * ncol:
        ax.set_xlabel(r"Expected $-\log_{10}(P)$", fontsize=7)
    else:
        ax.set_xlabel("")

    source_records.append({
        "marker_label": marker_label,
        "lambda_gc": lambda_gc,
        "min_p": np.min(p),
        "max_expected_neg_log10_p": np.nanmax(expected),
        "max_observed_neg_log10_p": np.nanmax(observed)
    })

for j in range(len(markers), len(axes)):
    axes[j].axis("off")

fig.subplots_adjust(
    left=0.07,
    right=0.99,
    bottom=0.045,
    top=0.99,
    wspace=0.35,
    hspace=0.45
)

fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
fig.savefig(out_png, dpi=600, bbox_inches="tight")
plt.close(fig)

source_data_supp_fig3 = pd.DataFrame(source_records)

with pd.ExcelWriter(out_source_data, engine="openpyxl") as writer:
    source_data_supp_fig3.to_excel(
        writer,
        sheet_name="source_data_supp_fig3",
        index=False
    )

print(f"Saved PDF: {out_pdf}")
print(f"Saved PNG: {out_png}")
print(f"Saved source data: {out_source_data}")