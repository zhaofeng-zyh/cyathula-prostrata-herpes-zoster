"""Standalone SI figures relocated from the consolidated main figures:
   FigureS6_splot   = PLS-DA S-plot (was Figure3_PLSDA panel D)
   FigureS7_lollipop = top-30 differential-metabolite lollipop (was Figure5 panel A)
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from _style import set_style, save_fig

set_style()
FIG = Path(__file__).resolve().parent.parent / "figures"

f3 = mpimg.imread(FIG / "Figure3_PLSDA.png")        # 5596x4474, 2x2 grid
splot = f3[2360:4474, 2980:5596]                    # bottom-right quadrant = S-plot (trim left bleed)
fig, ax = plt.subplots(figsize=(5.2, 4.2)); ax.imshow(splot); ax.axis("off")
# re-add the y-axis title lost by the crop
ax.text(-0.045, 0.5, "p_corr[1]  (correlation loading)", transform=ax.transAxes,
        rotation=90, ha="center", va="center", fontsize=10, fontweight="bold")
save_fig(fig, "FigureS6_splot"); plt.close(fig)

lol = mpimg.imread(FIG / "Figure5_lollipop_kegg.png")  # 9240x4380; A lollipop (left), B KEGG (right)
lollipop = lol[150:4380, 0:4360]                       # left panel
fig, ax = plt.subplots(figsize=(6.5, 5.8)); ax.imshow(lollipop); ax.axis("off")
save_fig(fig, "FigureS7_lollipop"); plt.close(fig)
print("FigureS6_splot, FigureS7_lollipop written")
