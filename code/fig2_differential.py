"""
Consolidated Figure 2 — Multivariate discrimination and differential metabolites,
C. prostrata leaf (Y6-13) vs stem (J6-13).

  (A) PLS-DA score plot           <- Figure3_PLSDA.png, top-right quadrant
  (B) 200-permutation validation  <- Figure3_PLSDA.png, bottom-left quadrant
  (C) Volcano plot                <- Figure4_volcano_heatmap.png, left panel
  (D) Top-50 z-score heatmap       <- Figure4_volcano_heatmap.png, right panel

The wanted sub-panels are sliced out of the two source PNGs BELOW their baked-in
inner panel-letter/title lines, so this montage carries a single clean outer
A/B/C/D scheme (Arial bold, 15 pt) with re-added concise panel titles.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from _style import set_style, save_fig

set_style()
FIG = Path(__file__).resolve().parent.parent / "figures"

f3 = mpimg.imread(FIG / "Figure3_PLSDA.png")          # 5596 x 4474, 2x2 grid
f4 = mpimg.imread(FIG / "Figure4_volcano_heatmap.png")  # 7350 x 3712

# --- crops (rows y0:y1, cols x0:x1); title bands sliced off ------------------
score = f3[130:2200, 2805:5596]   # top-right quadrant: PLS-DA score (keep y-axis title)
perm  = f3[2405:4474, 0:2720]     # bottom-left quadrant: permutation
volc  = f4[205:3360, 0:3060]      # left panel: volcano
heat  = f4[130:3690, 3215:7350]   # right panel: heatmap (first-column tick kept, no volcano bleed)

panels = [
    (score, "A", "PLS-DA score plot"),
    (perm,  "B", "200-permutation validation"),
    (volc,  "C", "Volcano plot (leaf vs stem)"),
    (heat,  "D", "Top-50 differential-metabolite heatmap"),
]

fig = plt.figure(figsize=(7.2, 7.3))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.04], height_ratios=[1.0, 1.30],
                      wspace=0.04, hspace=0.13)
axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]

for ax, (img, lab, title) in zip(axes, panels):
    ax.imshow(img)
    ax.axis("off")
    ax.text(0.0, 1.015, lab, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="bottom", ha="left")
    ax.text(0.5, 1.02, title, transform=ax.transAxes, fontsize=8.5,
            fontweight="bold", va="bottom", ha="center")

save_fig(fig, "Figure2_differential")
plt.close(fig)
print("Figure2_differential written")
