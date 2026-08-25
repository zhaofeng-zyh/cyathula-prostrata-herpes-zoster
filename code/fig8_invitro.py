"""
Consolidated Figure 8 — In-vitro validation in LPS-stimulated RAW 264.7 cells.
  (A) CCK-8 cell-viability dose screen   (Figure14_CCK8_viability.png)
  (B) TNF-alpha ELISA inhibition          (Figure15_TNFa_ELISA.png)

Both source panels are single-axis bar charts with their own descriptive
titles and no baked-in inner panel letters, so they are montaged whole and
stacked vertically with a clean outer A/B scheme (Arial bold, 15 pt, corner).
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from _style import set_style, save_fig

set_style()
FIG = Path(__file__).resolve().parent.parent / "figures"

a = mpimg.imread(FIG / "Figure14_CCK8_viability.png")   # viability
b = mpimg.imread(FIG / "Figure15_TNFa_ELISA.png")       # TNF-a ELISA

r_a = a.shape[1] / a.shape[0]
r_b = b.shape[1] / b.shape[0]

W = 7.2
h_a = W / r_a
h_b = W / r_b

fig = plt.figure(figsize=(W, h_a + h_b))
gs = fig.add_gridspec(2, 1, height_ratios=[h_a, h_b], hspace=0.05)
axA = fig.add_subplot(gs[0]); axA.imshow(a); axA.axis("off")
axB = fig.add_subplot(gs[1]); axB.imshow(b); axB.axis("off")

lab_kw = dict(fontsize=15, fontweight="bold", va="top", ha="left")


def panel_letter(ax, img, lab):
    """Place a clean outer letter in a small white left-gutter, clear of the
    axis title text that hugs the left edge of some source panels."""
    H, W = img.shape[:2]
    ax.set_xlim(-0.065 * W, W - 0.5)   # add ~6.5% white gutter on the left
    ax.text(0.004, 0.99, lab, transform=ax.transAxes, **lab_kw)


panel_letter(axA, a, "A")
panel_letter(axB, b, "B")

save_fig(fig, "Figure8_invitro")
plt.close(fig)
print("Figure8_invitro written")
