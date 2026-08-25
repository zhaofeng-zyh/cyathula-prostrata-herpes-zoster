"""
Consolidated Figure 3 — Network pharmacology and functional enrichment.
  (A) Drug x disease target intersection (Venn, 22 hubs)  <- Figure6_venn.png (left panel)
  (B) PPI network of the 22 hub targets                    <- Figure7_PPI.png
  (C) Drug-component-target (DCT) network                  <- Figure8_DCT.png
  (D) GO and KEGG enrichment of the 22 hub targets         <- Figure9_GO_KEGG.png

Sub-images are sliced from the source PNGs below/around their baked-in inner
titles so this montage carries a single clean outer A/B/C/D scheme.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from _style import set_style, save_fig

set_style()
FIG = Path(__file__).resolve().parent.parent / "figures"

venn = mpimg.imread(FIG / "Figure6_venn.png")   # 6419 x 2962; A Venn (left), B list (right)
ppi  = mpimg.imread(FIG / "Figure7_PPI.png")    # 4418 x 4416
dct  = mpimg.imread(FIG / "Figure8_DCT.png")    # 5407 x 6798 (tall)
enr  = mpimg.imread(FIG / "Figure9_GO_KEGG.png")  # 9287 x 5422 (wide, 4 inner panels)

venn_crop = venn[300:2760, 0:2900]              # left Venn only, BELOW the inner "A ..." title

# B/C/D source PNGs already carry their own descriptive titles → add outer letter only.
# A (Venn) title was cropped off → give it a fresh title.
panels = [
    (venn_crop, "A", None),
    (ppi,       "B", None),
    (dct,       "C", None),
    (enr,       "D", None),
]

fig = plt.figure(figsize=(7.2, 8.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.18], height_ratios=[1.0, 1.10],
                      wspace=0.05, hspace=0.12)
axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]

for ax, (img, lab, title) in zip(axes, panels):
    ax.imshow(img)
    ax.axis("off")
    ax.text(-0.02, 1.01, lab, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="bottom", ha="left")
    if title:
        ax.text(0.55, 1.02, title, transform=ax.transAxes, fontsize=8,
                fontweight="bold", va="bottom", ha="center")

save_fig(fig, "Figure3_network_enrichment")
plt.close(fig)
print("Figure3_network_enrichment written")
