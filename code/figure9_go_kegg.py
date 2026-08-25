"""
Figure 9.  GO + KEGG enrichment of the 22 hub targets.
A — GO BP, B — GO CC, C — GO MF (top terms by -log10 P)
D — KEGG pathway enrichment (human pathways only; matches Fig. 5B)
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig, COLORS_PATHWAY

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA = HERE.parent / "data"
set_style()

bp = pd.read_csv(DATA / "go_bp.csv")
cc = pd.read_csv(DATA / "go_cc.csv")
mf = pd.read_csv(DATA / "go_mf.csv")
kegg = pd.read_csv(DATA / "kegg_22targets.csv").sort_values("Pvalue", ascending=False).reset_index(drop=True)

for d in (bp, cc, mf):
    d["neglog10P"] = -np.log10(d["Pvalue"])
kegg["neglog10P"] = -np.log10(kegg["Pvalue"])

fig = plt.figure(figsize=(15.0, 9.6))
gs = fig.add_gridspec(2, 2, hspace=0.60, wspace=0.60,
                       top=0.94, bottom=0.07, left=0.05, right=0.96)
axes = [fig.add_subplot(gs[i,j]) for i in range(2) for j in range(2)]

def wrap_term(t, max_len=34):
    """Wrap long GO term names onto two lines for readability."""
    if len(t) <= max_len:
        return t
    # break at last space before max_len
    cut = t.rfind(" ", 0, max_len)
    if cut == -1: cut = max_len
    return t[:cut] + "\n" + t[cut+1:]

def hbar(ax, df, color, title):
    df = df.sort_values("neglog10P", ascending=True)
    y = np.arange(len(df))
    ax.barh(y, df["neglog10P"], color=color, edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([wrap_term(t) for t in df["Term"]], fontsize=7)
    ax.set_xlabel("$-\log_{10}$(P)")
    ax.set_title(title, loc="left", fontweight="bold", pad=10)
    ax.grid(axis="x", lw=0.4, ls="--", alpha=0.4)
    # Generous x headroom so count labels never collide with adjacent panel
    ax.set_xlim(0, df["neglog10P"].max() * 1.18)
    for j_, (_, r) in enumerate(df.iterrows()):
        ax.text(r["neglog10P"] + df["neglog10P"].max()*0.02, j_,
                f"{int(r['Count'])}",
                va="center", fontsize=6.4, color="#444")

hbar(axes[0], bp.head(10), "#1F77B4", "A   GO Biological Process (top 10)")
hbar(axes[1], cc.head(8),  "#3CA46E", "B   GO Cellular Component (top 8)")
hbar(axes[2], mf.head(8),  "#FF7F0E", "C   GO Molecular Function (top 8)")

# Panel D — KEGG bubble (human pathways only, plant-secondary already excluded)
ax = axes[3]
y = np.arange(len(kegg))
sizes = (kegg["Count"] - kegg["Count"].min() + 4) * 18
sc = ax.scatter(kegg["Fold_enrichment"], y, s=sizes,
                 c=kegg["neglog10P"], cmap="viridis",
                 edgecolor="white", linewidth=0.6)
ax.set_yticks(y)
ax.set_yticklabels([wrap_term(t, max_len=30) for t in kegg["Pathway"]], fontsize=7)
ax.set_xlabel("Fold enrichment")
ax.set_title("D   KEGG pathway enrichment (top 10)", loc="left", fontweight="bold", pad=10)
ax.set_ylim(-0.6, len(kegg) - 0.4)
ax.set_xlim(0, kegg["Fold_enrichment"].max() * 1.22)
ax.grid(axis="x", lw=0.4, ls="--", alpha=0.4)
cb = fig.colorbar(sc, ax=ax, label="$-\log_{10}$(P)", pad=0.02, shrink=0.7)
cb.outline.set_linewidth(0.3)
# Size legend
for cnt in [5, 10, 15]:
    s = (cnt - kegg["Count"].min() + 4) * 18
    ax.scatter([], [], s=s, color="lightgray", edgecolor="white", linewidth=0.6,
               label=f"{cnt} hits")
ax.legend(loc="lower right", title="Hit count", fontsize=6.5,
          title_fontsize=7, labelspacing=1.2, borderpad=0.6, frameon=False)

save_fig(fig, "Figure9_GO_KEGG")
print(f"GO BP: {len(bp)} terms; CC: {len(cc)}; MF: {len(mf)}; KEGG: {len(kegg)}")
