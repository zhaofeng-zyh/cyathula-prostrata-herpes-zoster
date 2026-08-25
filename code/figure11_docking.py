"""
Figure 11.  AutoDock Vina binding affinities of 18 representative active
ingredients (+ 2 reference drugs) against four canonical herpes-zoster
pathway targets. Reference drugs separated by a horizontal divider; row
labels rebuilt to eliminate the "Reference drugs / Acyclovir" overlap of
the original.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

DATA = HERE.parent / "data"
set_style()

dock = pd.read_csv(DATA / "docking_scores_real.csv")  # real AutoDock Vina 1.2.7 results
dock = dock.rename(columns={"Unnamed: 0": "Compound"})
# Strip the "— ... (ref)" prefix for cleanliness
dock["Compound_clean"] = dock["Compound"].str.replace(r"^—\s*", "", regex=True)\
                                          .str.replace(r"\s*\(ref\)$", "", regex=True)
dock["is_ref"] = dock["Compound"].str.contains("ref")

# Sort: lead compounds first (Arnicolide D, Salvigenin, Decursinol angelate, Carnosol),
# then other actives by mean ΔG, then references at bottom
target_cols = ["TNF-α (2AZ5)", "IL-6 (1ALU)", "NF-κB p65 (1IKN)", "COX-2 (5IKR)"]
dock["mean_dg"] = dock[target_cols].mean(axis=1)

leads = ["Arnicolide D", "Salvigenin", "Decursinol angelate", "Carnosol"]
non_ref = dock[~dock["is_ref"]].copy()
non_ref["lead"] = non_ref["Compound_clean"].isin(leads)
non_ref = non_ref.sort_values(["lead", "mean_dg"], ascending=[False, True])
refs = dock[dock["is_ref"]].sort_values("mean_dg", ascending=True)
ordered = pd.concat([non_ref, refs], ignore_index=True)

mat = ordered[target_cols].astype(float).values
labels = ordered["Compound_clean"].tolist()
is_ref = ordered["is_ref"].values
ref_start = int(np.argmax(is_ref))  # first row index of references

fig, ax = plt.subplots(figsize=(7.0, 9.0))
hm = sns.heatmap(mat, ax=ax, cmap="RdBu", vmin=-9, vmax=1.5, center=-4,
                 annot=True, fmt=".2f", annot_kws=dict(fontsize=7.0, fontweight="bold"),
                 cbar_kws=dict(label="AutoDock Vina ΔG (kcal/mol)",
                                shrink=0.55, pad=0.02),
                 linewidths=0.6, linecolor="white",
                 xticklabels=target_cols, yticklabels=labels)
ax.set_xticklabels(ax.get_xticklabels(), rotation=20, ha="right", fontsize=7)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)
ax.set_title("Molecular docking of C. prostrata active ingredients\nagainst four herpes-zoster pathway targets",
             fontweight="bold", fontsize=10)

# Outline cells with ΔG ≤ −9 kcal/mol
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if mat[i, j] <= -8.0:
            ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                                        edgecolor="#1F8F33", linewidth=2.0))

# Divider line between active compounds and references
ax.axhline(ref_start, color="black", lw=1.2)
ax.text(-0.7, ref_start - 0.4, "Reference",
        fontsize=7, fontweight="bold", color="#666", rotation=0,
        ha="right", va="center")
ax.text(-0.7, ref_start + 1.4, "drugs ↓",
        fontsize=7, fontweight="bold", color="#666", rotation=0,
        ha="right", va="center")

# Star/marker for the four leads
for i, name in enumerate(labels):
    if name in leads:
        ax.text(mat.shape[1] + 0.18, i + 0.5, r"$\star$", color="#cc6600", fontsize=15,
                ha="left", va="center", fontweight="bold")

# Legend
green_patch = mpatches.Patch(facecolor="none", edgecolor="#1F8F33", linewidth=2.0,
                              label="ΔG ≤ −8 kcal/mol")
star = plt.Line2D([], [], marker=r"$\star$", color="#cc6600", linestyle="",
                   markersize=12, label="Chemotype-defining lead")
ax.legend(handles=[green_patch, star], loc="upper right",
           bbox_to_anchor=(1.02, 1.10), fontsize=7, frameon=False)

plt.tight_layout()
save_fig(fig, "Figure11_docking")
print(f"Docking matrix: {mat.shape[0]} compounds × {mat.shape[1]} targets")
print(f"Cells with ΔG ≤ -8: {int((mat <= -8).sum())}")
