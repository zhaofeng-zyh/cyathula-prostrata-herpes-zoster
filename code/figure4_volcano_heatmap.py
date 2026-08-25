"""
Figure 4.  Differential-metabolite landscape.
A — Volcano plot of log2 FC vs -log10 P with adjustText label collision avoidance
B — Z-score heatmap of top-50 differential metabolites (Ward's linkage)

Fixes vs. original:
- Volcano labels limited to top-12 by VIP × |log2FC|, with arrow connectors
- Heatmap title moved up; metabolite labels font size reduced; sample labels rotated
- Legend explicitly shows "Up / Down / NS-by-FC" with counts
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig, COLOR_STEM, COLOR_LEAF

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, leaves_list

try:
    from adjustText import adjust_text
    HAS_ADJ = True
except ImportError:
    HAS_ADJ = False

DATA = HERE.parent / "data"
set_style()

df = pd.read_csv(DATA / "full_metabolites.csv")
samples = ["J6-13-1","J6-13-2","J6-13-3","Y6-13-1","Y6-13-2","Y6-13-3"]
df["log2FC"] = df["log2(Y6-13/J6-13)"].astype(float)
df["neglogP"] = -np.log10(df["P-value"].astype(float).clip(lower=1e-12))

# Categorise
LFC = 1.0
sig_up   = (df["log2FC"] >=  LFC)
sig_down = (df["log2FC"] <= -LFC)
ns       = ~(sig_up | sig_down)
n_up, n_down, n_ns = int(sig_up.sum()), int(sig_down.sum()), int(ns.sum())

# Top-12 by VIP × |log2FC|
df["score"] = df["VIP"].astype(float) * np.abs(df["log2FC"])
top12 = df.nlargest(9, "score")

# Top-50 for heatmap (by score)
top50 = df.nlargest(50, "score").copy()

# --- Build figure ------------------------------------------------------------
fig = plt.figure(figsize=(13.0, 6.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.95], wspace=0.05)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])

# Panel A — Volcano
axA.scatter(df.loc[ns, "log2FC"],   df.loc[ns, "neglogP"],
            s=8, color="#bbbbbb", alpha=0.5, label=f"|$\log_2$FC| < 1  (n = {n_ns})")
axA.scatter(df.loc[sig_down, "log2FC"], df.loc[sig_down, "neglogP"],
            s=10, color=COLOR_STEM, alpha=0.65, label=f"Down in leaf  (n = {n_down})")
axA.scatter(df.loc[sig_up,   "log2FC"], df.loc[sig_up,   "neglogP"],
            s=10, color=COLOR_LEAF, alpha=0.65, label=f"Up in leaf  (n = {n_up})")
axA.axvline(-LFC, color="gray", lw=0.5, ls="--", alpha=0.6)
axA.axvline( LFC, color="gray", lw=0.5, ls="--", alpha=0.6)
axA.axhline(-np.log10(0.05), color="gray", lw=0.5, ls="--", alpha=0.6)

# Annotate top-12
texts = []
for _, r in top12.iterrows():
    name = (r["Name"] if pd.notna(r["Name"]) and r["Name"] else f"Compound {int(r.name)}")
    short = name if len(name) <= 26 else name[:24] + "…"
    t = axA.text(r["log2FC"], r["neglogP"], short, fontsize=6.5, color="#222")
    texts.append(t)
if HAS_ADJ and texts:
    adjust_text(texts, ax=axA,
                expand_text=(1.3, 1.6), expand_points=(1.3, 1.5),
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, alpha=0.6))

axA.set_xlabel("$\log_2$ fold change  (Y6-13 / J6-13)")
axA.set_ylabel("$-\log_{10}$(P-value)")
axA.set_title("A   Volcano plot of differential metabolites", loc="left", fontweight="bold")
axA.legend(loc="upper left", fontsize=7, frameon=True, framealpha=0.9)
axA.grid(False)

# Panel B — Heatmap
abund = top50[samples].astype(float).values
log_ab = np.log10(abund + 1)
z = (log_ab - log_ab.mean(axis=1, keepdims=True)) / log_ab.std(axis=1, keepdims=True).clip(min=1e-9)

# Ward linkage on rows (metabolites)
link = linkage(z, method="ward", metric="euclidean")
order = leaves_list(link)
z_ord = z[order]

# Truncate metabolite labels to 28 chars
metab_labels = []
for j in order:
    n = top50["Name"].iloc[j]
    n = str(n) if pd.notna(n) else ""
    metab_labels.append(n if len(n) <= 28 else n[:26] + "…")

import seaborn as sns
# Build heatmap without seaborn's default colorbar; we add our own with
# explicit positioning so the y-tick labels (right side) sit cleanly
# between the heatmap and the colorbar.
hm = sns.heatmap(z_ord, ax=axB, cmap="RdBu_r", center=0, vmin=-2.2, vmax=2.2,
                 cbar=False, linewidths=0,
                 xticklabels=samples, yticklabels=metab_labels)
axB.set_title("B   Z-score heatmap of top-50 differential metabolites",
              loc="left", x=-0.02, y=1.02, fontweight="bold")
axB.set_xticklabels(axB.get_xticklabels(), rotation=45, ha="right", fontsize=7)
# Move metabolite labels to the RIGHT side of the heatmap so they extend
# away from Panel A (volcano) rather than colliding with it.
axB.yaxis.tick_right()
axB.yaxis.set_label_position("right")
axB.set_yticklabels(axB.get_yticklabels(), rotation=0, fontsize=5.6, ha="left")
# Add a separate colorbar to the FAR right, well beyond the labels.
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
norm = plt.Normalize(vmin=-2.2, vmax=2.2)
sm = cm.ScalarMappable(cmap="RdBu_r", norm=norm)
cax = inset_axes(axB, width="2.5%", height="55%",
                  loc="center left",
                  bbox_to_anchor=(1.42, 0.0, 1, 1),
                  bbox_transform=axB.transAxes, borderpad=0)
cb = fig.colorbar(sm, cax=cax)
cb.set_label("Row Z-score", fontsize=8)
cb.ax.tick_params(labelsize=7)
cb.outline.set_linewidth(0.3)
# Banners drawn BELOW the heatmap and below the rotated sample labels
# using figure-fraction coordinates (transFigure) so they cannot collide
# with the data axes.
nrows = z_ord.shape[0]
# Compute axes-coordinate positions for two horizontal bars under the x-axis
# In axes coordinates: y=0 is the top of x-axis; y=1 is the top of heatmap.
# We place banner at axes-y = -0.20 (well below 45° rotated tick labels).
for (x0, x1, color, lbl) in [(0/6, 3/6, COLOR_STEM, "Stem (J6-13)"),
                              (3/6, 6/6, COLOR_LEAF, "Leaf (Y6-13)")]:
    # short padding from group boundary
    pad = 0.01
    axB.plot([x0 + pad, x1 - pad], [-0.22, -0.22], color=color,
             lw=2.4, clip_on=False, solid_capstyle="round",
             transform=axB.transAxes)
    axB.text((x0 + x1) / 2, -0.255, lbl,
             ha="center", va="top", color=color,
             fontsize=8.5, fontweight="bold", clip_on=False,
             transform=axB.transAxes)

plt.subplots_adjust(left=0.06, right=0.78, top=0.93, bottom=0.22, wspace=0.05)
save_fig(fig, "Figure4_volcano_heatmap")
print(f"Volcano: up={n_up}  down={n_down}  ns(by FC)={n_ns}  total={n_up+n_down+n_ns}")
