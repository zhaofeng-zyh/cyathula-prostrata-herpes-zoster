"""
Consolidated Figure 1 — Untargeted metabolome overview of C. prostrata
stem (J6-13) vs leaf (Y6-13).

Montage of two existing 600-dpi panels stacked vertically:
  Row 1  = Figure1_classification_correlation.png  ->  (A) chemical-taxonomy
           donut  +  (B) inter-sample Pearson correlation heatmap
  Row 2  = Figure2_PCA.png                          ->  (C) PCA ESI-  +  (D) PCA ESI+

The two source scripts (figure1_classification_correlation.py, figure2_pca.py)
were edited to drop their baked-in inner "A"/"B" panel-letter prefixes, so this
montage carries a single clean outer A/B/C/D scheme (Arial bold, 15 pt, corner).
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from _style import set_style, save_fig

set_style()
FIG = Path(__file__).resolve().parent.parent / "figures"

top = mpimg.imread(FIG / "Figure1_classification_correlation.png")   # donut | heatmap
bot = mpimg.imread(FIG / "Figure2_PCA.png")                          # PCA- | PCA+

r_top = top.shape[1] / top.shape[0]   # w/h
r_bot = bot.shape[1] / bot.shape[0]

W = 7.2
h_top = W / r_top
h_bot = W / r_bot

fig = plt.figure(figsize=(W, h_top + h_bot))
gs = fig.add_gridspec(2, 1, height_ratios=[h_top, h_bot], hspace=0.04)
axT = fig.add_subplot(gs[0]); axT.imshow(top); axT.axis("off")
axB = fig.add_subplot(gs[1]); axB.imshow(bot); axB.axis("off")

lab_kw = dict(fontsize=15, fontweight="bold", va="top", ha="left")
# Row 1 : A over donut, B over heatmap
axT.text(0.005, 1.00, "A", transform=axT.transAxes, **lab_kw)
axT.text(0.585, 1.00, "B", transform=axT.transAxes, **lab_kw)
# Row 2 : C over ESI- panel, D over ESI+ panel (kept below the shared top legend)
axB.text(0.005, 0.90, "C", transform=axB.transAxes, **lab_kw)
axB.text(0.515, 0.90, "D", transform=axB.transAxes, **lab_kw)

save_fig(fig, "Figure1_metabolome")
plt.close(fig)
print("Figure1_metabolome written")
