"""
Figure 8.  Drug–component–target (DCT) network built from SwissTargetPrediction
output (data/dct_featured_edges.csv). Source nodes are the named docking/lead
representatives; target nodes are their predicted targets, grouped into six
herpes-zoster-relevant pharmacological categories. Target node size = number of the
121 profiled library compounds predicting that target (library-wide support).
(Replaces the earlier version, whose edges were heuristically assigned.)
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

DATA = HERE.parent / "data"
set_style()
df = pd.read_csv(DATA / "dct_featured_edges.csv")

CAT_ORDER = ["Eicosanoid / inflammation", "NF-κB / MAPK / JAK-STAT",
             "Immune / cytokine / MMP", "Pain / nociception channels",
             "Opioid / cannabinoid / adenosine", "Nucleotide / antiviral"]
CAT_COLOR = {"Eicosanoid / inflammation": "#C03A2B", "NF-κB / MAPK / JAK-STAT": "#E07B39",
             "Immune / cytokine / MMP": "#E6B800", "Pain / nociception channels": "#2CA02C",
             "Opioid / cannabinoid / adenosine": "#1F77B4", "Nucleotide / antiviral": "#7B4F9A"}

g2c = df.drop_duplicates("Gene").set_index("Gene")["Category"].to_dict()
g2lib = df.drop_duplicates("Gene").set_index("Gene")["LibrarySupport"].to_dict()
compounds = sorted(df["Compound"].unique(), key=lambda c: -df[df.Compound == c].shape[0])
cdeg = df.groupby("Compound")["Gene"].nunique().to_dict()

# target y-positions: stack by category band with gaps
genes_by_cat = {c: sorted([g for g in g2c if g2c[g] == c], key=lambda g: -g2lib[g]) for c in CAT_ORDER}
gy = {}; y = 0.0; band_mid = {}; GAP = 1.4
for c in CAT_ORDER:
    gs = genes_by_cat[c]; ys = []
    for g in gs:
        gy[g] = y; ys.append(y); y += 1.0
    band_mid[c] = np.mean(ys) if ys else y
    y += GAP
TOTAL = y - GAP
# compound y-positions spread over same vertical extent
cy = {c: (TOTAL/(len(compounds)-1))*i if len(compounds) > 1 else TOTAL/2 for i, c in enumerate(compounds)}
xC, xG = 0.0, 5.0

fig, ax = plt.subplots(figsize=(11.5, 13.5))
for _, r in df.iterrows():
    y1 = cy[r.Compound]; y2 = gy[r.Gene]
    p = float(r.Probability) if pd.notna(r.Probability) else 0.1
    ax.plot([xC, xG], [y1, y2], color=CAT_COLOR[r.Category], alpha=0.4,
            lw=0.5 + 1.6*p, zorder=1)
for c in compounds:
    ax.scatter(xC, cy[c], s=120+55*cdeg[c], color="#3a3a3a", zorder=3,
               edgecolor="white", linewidth=0.7)
    ax.text(xC-0.18, cy[c], c, ha="right", va="center", fontsize=8, fontweight="bold")
for g in gy:
    col = CAT_COLOR[g2c[g]]
    ax.scatter(xG, gy[g], s=min(70+34*g2lib[g], 760), color=col, zorder=3,
               edgecolor="white", linewidth=0.7)
    ax.text(xG+0.55, gy[g], f"{g}", ha="left", va="center", fontsize=7.8,
            color=col, fontweight="bold")
# category band labels
for c in CAT_ORDER:
    ax.text(xG+2.4, band_mid[c], c, ha="left", va="center", fontsize=8.2,
            color=CAT_COLOR[c], fontweight="bold", rotation=0)

ax.text(xC, TOTAL+1.6, "Lead / representative\nactive ingredients", ha="center",
        va="bottom", fontsize=9.5, fontweight="bold")
ax.text(xG+0.6, TOTAL+1.6, "Predicted targets  (node size ~ no. of library\ncompounds predicting the target)",
        ha="center", va="bottom", fontsize=9.5, fontweight="bold")
handles = [mlines.Line2D([], [], color=CAT_COLOR[c], marker="o", linestyle="-",
            markersize=8, label=c) for c in CAT_ORDER]
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.05),
          ncol=3, fontsize=8, frameon=False)
ax.set_title("Drug–component–target network (SwissTargetPrediction)",
             fontsize=12, fontweight="bold", pad=18)
ax.set_xlim(-2.2, 11.0); ax.set_ylim(-2.2, TOTAL+3.0)
ax.axis("off")
save_fig(fig, "Figure8_DCT")
print(f"DCT network: {len(compounds)} compounds, {len(gy)} targets, {len(df)} edges, 6 categories")
