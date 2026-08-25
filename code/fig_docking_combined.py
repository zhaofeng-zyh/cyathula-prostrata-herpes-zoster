"""
Combined molecular-docking figure (Fig 5) — 3-panel layout (restored per user request):
  (A) AutoDock Vina affinity heatmap                (Figure11_docking.png)
  (B) Decursinol angelate - COX-2 50-ns MD 3D pose  (docking_poses composite, panel C)
  (C) Decursinol angelate - COX-2 2D interaction map(docking_poses composite, panel D)

Panels B/C are cropped from the 4-panel PyMOL/PLIP composite below its baked-in
inner "C"/"D" titles; clean outer A/B/C labels + short subtitles are added here.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from _style import set_style, save_fig

set_style()
FIG = Path(__file__).resolve().parent.parent / "figures"

heat = mpimg.imread(FIG / "Figure11_docking.png")
comp = mpimg.imread(FIG / "docking_poses" / "DecursinolAng_COX2_composite.png")
H, W = comp.shape[:2]
# composite = 2x2 grid; bottom-left = 3D pose (panel C), bottom-right = 2D map (panel D).
# crop below each baked inner title, and inside the central divider.
# crop each bottom panel BELOW its baked inner "C"/"D" title (~0.55H) and inside the divider
pose = comp[int(H * 0.565):, 0:int(W * 0.49)]
imap = comp[int(H * 0.585):, int(W * 0.515):]

fig = plt.figure(figsize=(7.6, 5.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1.0], height_ratios=[1.0, 1.0],
                      wspace=0.02, hspace=0.20)
axA = fig.add_subplot(gs[:, 0]); axA.imshow(heat); axA.axis("off")
axB = fig.add_subplot(gs[0, 1]); axB.imshow(pose); axB.axis("off")
axC = fig.add_subplot(gs[1, 1]); axC.imshow(imap); axC.axis("off")

axA.text(0.02, 1.00, "A", transform=axA.transAxes, fontsize=15, fontweight="bold", va="top", ha="left")
axB.text(0.0, 1.10, "B", transform=axB.transAxes, fontsize=15, fontweight="bold", va="top", ha="left")
axB.text(0.56, 1.10, "3D binding pose (50-ns MD)", transform=axB.transAxes,
         fontsize=8, ha="center", va="top", fontweight="bold")
axC.text(0.0, 1.10, "C", transform=axC.transAxes, fontsize=15, fontweight="bold", va="top", ha="left")
axC.text(0.56, 1.10, "2D interaction map", transform=axC.transAxes,
         fontsize=8, ha="center", va="top", fontweight="bold")

save_fig(fig, "Figure_docking_combined")
plt.close(fig)
print("Figure_docking_combined written (heatmap + 3D pose + 2D interaction map)")
