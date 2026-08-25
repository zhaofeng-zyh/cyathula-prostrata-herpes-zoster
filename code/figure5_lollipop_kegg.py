"""
Figure 5.
A — Lollipop plot of the top-30 differential metabolites by absolute log2FC
B — KEGG pathway enrichment of the 22 hub targets

Fixes vs. original:
- Plant-secondary pathways REMOVED from KEGG enrichment (they had P=1.0 and
  are not biologically interpretable as enrichments of human hub targets)
- Long lipid names truncated; vertical layout stretched so labels don't collide
- Bubble plot uses correct fold-enrichment scale and viridis colormap
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig, COLOR_STEM, COLOR_LEAF, COLORS_PATHWAY

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

DATA = HERE.parent / "data"
set_style()

# --- Panel A: top-30 metabolites by |log2FC| --------------------------------
top30 = pd.read_csv(DATA / "top30_metabolites.csv")
top30 = top30.rename(columns={"英文名": "Name_EN", "log2(Y6-13/J6-13)": "log2FC"})
top30 = top30.sort_values("log2FC", ascending=True).reset_index(drop=True)

# Truncate long names
def short(n, k=30):
    s = str(n) if pd.notna(n) else ""
    return s if len(s) <= k else s[:k-1] + "…"
top30["label"] = top30["Name_EN"].apply(lambda s: short(s, 30))

# --- Panel B: KEGG enrichment of 22 hub targets (human pathways only) -------
kegg = pd.read_csv(DATA / "kegg_22targets.csv")  # this CSV is human-only
kegg = kegg.sort_values("Pvalue", ascending=False).reset_index(drop=True)
kegg["neglog10P"] = -np.log10(kegg["Pvalue"])

# Pathway category mapping (manual, based on biology)
pathway_cat = {
    "Cytokine–cytokine receptor interaction": "Inflammation/Immunity",
    "TNF signalling pathway":                  "Inflammation/Immunity",
    "NF-κB signalling pathway":                "Inflammation/Immunity",
    "Toll-like receptor signalling":           "Inflammation/Immunity",
    "IL-17 signalling pathway":                "Inflammation/Immunity",
    "Chemokine signalling pathway":            "Inflammation/Immunity",
    "T-cell receptor signalling":              "Inflammation/Immunity",
    "JAK-STAT signalling pathway":             "Signal transduction",
    "Apoptosis":                               "Apoptosis",
    "Herpes simplex virus 1 infection":        "Antiviral",
}
kegg["Category"] = kegg["Pathway"].map(pathway_cat).fillna("Signal transduction")

# --- Layout ---
fig = plt.figure(figsize=(14.5, 8.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 0.95], wspace=0.05,
                       top=0.90, bottom=0.10, left=0.05, right=0.80)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])

# Panel A
y_pos = np.arange(len(top30))
colors_lol = np.where(top30["log2FC"] > 0, COLOR_LEAF, COLOR_STEM)
axA.hlines(y_pos, 0, top30["log2FC"], color=colors_lol, lw=1.6, alpha=0.85)
axA.scatter(top30["log2FC"], y_pos, s=45, color=colors_lol, edgecolor="white",
            linewidth=0.8, zorder=3)
for i, r in top30.iterrows():
    txt = f"{r['log2FC']:+.2f}"
    xt  = r["log2FC"] + (0.20 if r["log2FC"] > 0 else -0.20)
    axA.text(xt, i, txt, ha="left" if r["log2FC"] > 0 else "right",
             va="center", fontsize=6, color=colors_lol[i], fontweight="bold")
axA.set_yticks(y_pos)
axA.set_yticklabels(top30["label"], fontsize=6.5)
axA.axvline(0, color="gray", lw=0.6, ls="-", alpha=0.4)
axA.set_xlabel("$\log_2$ fold change  (Y6-13 / J6-13)")
axA.set_title("A   Top-30 differential metabolites by |$\log_2$FC|",
              loc="left", fontweight="bold", pad=14)
axA.set_ylim(-0.8, len(top30) - 0.2)
axA.set_xlim(min(top30["log2FC"].min() - 1.0, -6),
             max(top30["log2FC"].max() + 1.0,  9))
axA.grid(axis="x", lw=0.4, ls="--", alpha=0.4)
# Inline legend
axA.scatter([], [], s=60, color=COLOR_LEAF, label="Up in leaf (Y6-13)")
axA.scatter([], [], s=60, color=COLOR_STEM, label="Down in leaf (Y6-13)")
axA.legend(loc="lower right", fontsize=7)

# Panel B — bubble plot
y_pos2 = np.arange(len(kegg))
sizes = (kegg["Count"].astype(float) - kegg["Count"].min() + 4) * 22
colors2 = [COLORS_PATHWAY[c] for c in kegg["Category"]]
sc = axB.scatter(kegg["Fold_enrichment"], y_pos2, s=sizes,
                 c=kegg["neglog10P"], cmap="viridis",
                 edgecolor="white", linewidth=0.8, zorder=3)

# Pathway labels on the RIGHT side so they extend away from Panel A.
axB.set_yticks(y_pos2)
axB.set_yticklabels(kegg["Pathway"], fontsize=7)
axB.yaxis.tick_right()
axB.yaxis.set_label_position("right")
plt.setp(axB.get_yticklabels(), ha="left")
axB.set_xlabel("Fold enrichment")
axB.set_title("B   KEGG pathway enrichment (22 hub targets, human pathways only)",
              loc="left", fontweight="bold", pad=14)
axB.grid(axis="x", lw=0.4, ls="--", alpha=0.4)
# Expand y limits so the top pathway has clear headroom under the title
axB.set_ylim(-0.6, len(kegg) - 0.4)
axB.set_xlim(0, kegg["Fold_enrichment"].max() * 1.25)
# Add explicit padding on the right side of the tick labels so the colorbar
# (added later) doesn't sit on top of them
axB.tick_params(axis="y", which="major", pad=4)

# Annotate each bubble with -log10(P) — placed to LEFT of bubble so it
# stays inside the data area (yticklabels are now on the right)
for i, r in kegg.iterrows():
    axB.text(r["Fold_enrichment"] - kegg["Fold_enrichment"].max()*0.025,
             i, f"P = {r['Pvalue']:.1e}", fontsize=6.0, va="center",
             ha="right", color="#333")

# Colorbar for -log10(P) — placed FAR right, well beyond the y-tick labels
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
cax = inset_axes(axB, width="2.5%", height="55%",
                  loc="center left",
                  bbox_to_anchor=(1.45, 0.0, 1, 1),
                  bbox_transform=axB.transAxes, borderpad=0)
cb = fig.colorbar(sc, cax=cax)
cb.set_label("$-\log_{10}$(P)", fontsize=8)
cb.ax.tick_params(labelsize=7)
cb.outline.set_linewidth(0.3)

# Size legend for hit count
size_vals = [5, 10, 15]
for s_val in size_vals:
    s = (s_val - kegg["Count"].min() + 4) * 22
    axB.scatter([], [], s=s, color="gray", alpha=0.5,
                edgecolor="white", linewidth=0.8,
                label=f"{s_val} hits")
axB.legend(loc="lower right", title="Hit count", fontsize=6.5, title_fontsize=7,
           labelspacing=1.2, borderpad=0.7)

save_fig(fig, "Figure5_lollipop_kegg")
print(f"Lollipop: {len(top30)} metabolites; KEGG: {len(kegg)} human pathways (plant secondary REMOVED)")
