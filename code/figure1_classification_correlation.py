"""
Figure 1.  Chemical taxonomy of the 1,382 differential metabolites and
inter-sample reproducibility between the J6-13 (stem) and Y6-13 (leaf)
groups of Cyathula prostrata.

Caption-aligned with revised manuscript: panel A is the HMDB SuperClass
distribution of the *differential* set (n = 1,382), NOT the full 4,419-
metabolite catalogue.

Layout strategy: three independent axes positioned with explicit
add_axes() rectangles to guarantee no auto-layout collision; legend
typeset with monospace counts/percentages so columns align cleanly.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig, COLOR_STEM, COLOR_LEAF

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

DATA = HERE.parent / "data"
set_style()

df = pd.read_csv(DATA / "full_metabolites.csv")
samples = ["J6-13-1", "J6-13-2", "J6-13-3", "Y6-13-1", "Y6-13-2", "Y6-13-3"]

# Panel A — HMDB SuperClass donut
sc = df["ClassI"].fillna("Unclassified").replace({"_": "Unclassified", "": "Unclassified"})
zh_to_en = {
    "脂质和类脂分子": "Lipids and lipid-like molecules",
    "苯丙素类和聚酮类": "Phenylpropanoids & polyketides",
    "有机杂环化合物": "Organoheterocyclic compounds",
    "有机酸及其衍生物": "Organic acids & derivatives",
    "有机氧化合物": "Organic oxygen compounds",
    "苯类化合物": "Benzenoids",
    "生物碱及其衍生物": "Alkaloids & derivatives",
    "Unclassified": "Unclassified",
    "木脂素类、新木脂素类和相关化合物": "Lignans & neolignans",
    "核苷酸、核苷及其类似物": "Nucleosides & nucleotides",
    "有机氮化合物": "Organic nitrogen compounds",
}
sc_en = sc.map(lambda s: zh_to_en.get(s, s))
counts = sc_en.value_counts()
threshold = 15
keep = counts[counts >= threshold].copy()
other = counts[counts < threshold].sum()
if other > 0:
    keep["Other minor classes"] = other
keep = keep.sort_values(ascending=False)
total = int(keep.sum())

# Panel B — Pearson correlation
abund = df[samples].replace(0, np.nan)
log_ab = np.log10(abund + 1)
corr = log_ab.corr(method="pearson")

# --- Manual layout ---------------------------------------------------------
fig = plt.figure(figsize=(12.0, 4.6))

# Donut on the left (square)
axA  = fig.add_axes([0.02, 0.18, 0.26, 0.68])  # [left, bottom, width, height] in fig coords
# Legend column to the right of donut
axAL = fig.add_axes([0.30, 0.10, 0.32, 0.78]); axAL.axis("off")
# Heatmap on the far right
axB  = fig.add_axes([0.66, 0.16, 0.31, 0.72])

# --- Donut ----------------------------------------------------------------
palette = sns.color_palette("Spectral", n_colors=len(keep))
axA.pie(keep.values, radius=1.0, colors=palette,
        wedgeprops=dict(width=0.42, edgecolor="white", linewidth=1.2),
        startangle=90, counterclock=False)
axA.text(0, 0.10, f"{total:,}", ha="center", va="center",
         fontsize=15, fontweight="bold")
axA.text(0, -0.08, "differential\nmetabolites",
         ha="center", va="center", fontsize=7.5, color="#444")
axA.text(0, -0.30, "(VIP ≥ 1, P ≤ 0.05)",
         ha="center", va="center", fontsize=6.5, style="italic", color="#666")
axA.set_title("HMDB SuperClass distribution of differential metabolites",
              loc="left", x=0.0, y=1.04, fontsize=9, fontweight="bold")

# --- Legend (manual, three columns: swatch | name | n + pct) ---------------
# Header row
axAL.text(0.07, 0.95, "Class", transform=axAL.transAxes,
          fontsize=7.2, fontweight="bold", va="center")
axAL.text(0.78, 0.95, "n", transform=axAL.transAxes,
          fontsize=7.2, fontweight="bold", va="center", ha="right",
          family="monospace")
axAL.text(0.99, 0.95, "%", transform=axAL.transAxes,
          fontsize=7.2, fontweight="bold", va="center", ha="right",
          family="monospace")

y_top = 0.86
y_step = 0.075
for i, (color, (name, n)) in enumerate(zip(palette, keep.items())):
    pct = 100 * n / total
    y = y_top - i * y_step
    axAL.add_patch(mpl.patches.Rectangle((0.00, y - 0.022), 0.045, 0.038,
                                         transform=axAL.transAxes,
                                         color=color, ec="white", lw=0.4,
                                         clip_on=False))
    axAL.text(0.07, y, name, transform=axAL.transAxes,
              fontsize=7.0, va="center")
    axAL.text(0.78, y, f"{n:>4,}", transform=axAL.transAxes,
              fontsize=7.0, va="center", ha="right", family="monospace")
    axAL.text(0.99, y, f"{pct:>4.1f}", transform=axAL.transAxes,
              fontsize=7.0, va="center", ha="right", family="monospace")

# --- Heatmap --------------------------------------------------------------
sns.heatmap(corr, ax=axB, annot=True, fmt=".2f",
            annot_kws=dict(fontsize=6.8, fontweight="bold"),
            cmap="RdBu_r", vmin=0.0, vmax=1.0, center=0.6,
            cbar_kws=dict(label="Pearson r", shrink=0.78, pad=0.02),
            square=True, linewidths=0.6, linecolor="white")
axB.set_title("Inter-sample Pearson correlation",
              loc="left", x=-0.02, y=1.20, fontsize=9, fontweight="bold")
axB.set_xticklabels(axB.get_xticklabels(), rotation=45, ha="right", fontsize=7)
axB.set_yticklabels(axB.get_yticklabels(), rotation=0, fontsize=7)
# Group banners below title, above heatmap
axB.text(1.5, -0.55, "Stem (J6-13)", ha="center", va="bottom",
         color=COLOR_STEM, fontsize=7.8, fontweight="bold")
axB.text(4.5, -0.55, "Leaf (Y6-13)", ha="center", va="bottom",
         color=COLOR_LEAF, fontsize=7.8, fontweight="bold")
axB.plot([0.05, 2.95], [-0.20, -0.20], color=COLOR_STEM, lw=2.0, clip_on=False)
axB.plot([3.05, 5.95], [-0.20, -0.20], color=COLOR_LEAF, lw=2.0, clip_on=False)

save_fig(fig, "Figure1_classification_correlation")

within_stem = corr.iloc[:3, :3].values[np.triu_indices(3, k=1)].mean()
within_leaf = corr.iloc[3:, 3:].values[np.triu_indices(3, k=1)].mean()
across      = corr.iloc[:3, 3:].values.mean()
print(f"Donut total: {total}; classes shown: {len(keep)}")
print(f"Mean Pearson r — within-stem: {within_stem:.3f}  within-leaf: {within_leaf:.3f}  across: {across:.3f}")
