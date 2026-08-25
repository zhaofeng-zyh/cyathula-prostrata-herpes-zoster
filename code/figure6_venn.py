"""
Figure 6.  Drug-disease target intersection.
A — Venn diagram (432 ligand targets ∩ 54 disease genes = 22 hub targets)
B — Functional grouping of the 22 hub targets

Fixes vs. original:
- Footnote corrected: 168 (not 179) core active ingredients
- Five-category functional grouping aligned with manuscript text §5
- Categories: pro-inflammatory cytokines & chemokines / immunoregulatory cytokines /
  antigen presentation (MHC class I) / T-cell & lymphocyte signalling /
  apoptosis & innate-immunity effectors
- Group assignments biologically corrected (IL-6 NOT in apoptosis;
  HLA-A NOT in adaptive cytokine, etc.)
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn2, venn2_circles

DATA = HERE.parent / "data"
set_style()

# Corrected functional grouping (single source of truth)
HUB_GROUPS = {
    "Pro-inflammatory cytokines & chemokines": {
        "color":   "#C03A2B",
        "members": ["TNF", "IL1B", "IL6", "CXCL10", "CCL5", "CRP"],
    },
    "Immunoregulatory cytokines": {
        "color":   "#E69F00",
        "members": ["IL2", "IL4", "IL7", "IL10", "IFNB1"],
    },
    "Antigen presentation (MHC class I)": {
        "color":   "#7A1F75",
        "members": ["HLA-A", "HLA-B", "HLA-C"],
    },
    "T-cell / lymphocyte markers & signalling": {
        "color":   "#4F6FBE",
        "members": ["CD4", "CD8A", "CCR5", "CREB1"],
    },
    "Apoptosis / innate-immunity effectors": {
        "color":   "#3CA46E",
        "members": ["CASP8", "ADA", "APOE", "FCGR3A"],
    },
}
all_22 = [g for v in HUB_GROUPS.values() for g in v["members"]]
assert len(all_22) == 22, f"Expected 22 hub targets, got {len(all_22)}"
assert len(set(all_22)) == 22, "Duplicates in hub-target list"

# --- Layout ---
fig = plt.figure(figsize=(11.5, 5.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.10], wspace=0.05)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1]); axB.axis("off")

# --- Panel A: Venn ---
v = venn2(subsets=(410, 32, 22),
          set_labels=("",""),
          set_colors=("#56B4E9", "#C03A2B"),
          alpha=0.55, ax=axA)
# Customise label sizes
for lbl_id in ("10", "01", "11"):
    if v.get_label_by_id(lbl_id):
        v.get_label_by_id(lbl_id).set_fontsize(13)
        v.get_label_by_id(lbl_id).set_fontweight("bold")
        v.get_label_by_id(lbl_id).set_color("white")

# Outline circles
c = venn2_circles(subsets=(410, 32, 22), linewidth=0.8, color="white", ax=axA)

# Manual set labels and counts
axA.text(-0.65, -0.40, "C. prostrata\nactive-ingredient targets\n(n = 432)",
         ha="center", va="top", fontsize=8, color="#1f4d6b", fontweight="bold")
axA.text(0.55, -0.40, "Herpes zoster\ndisease targets\n(n = 54)",
         ha="center", va="top", fontsize=8, color="#7a1f1f", fontweight="bold")
axA.set_title("A   Drug-disease target intersection",
              loc="left", fontweight="bold")

# --- Panel B: 5-category grouping ---
y_top = 0.94
y_step = 0.18
for i, (grp, info) in enumerate(HUB_GROUPS.items()):
    y = y_top - i * y_step
    axB.add_patch(mpatches.Rectangle((0.02, y - 0.13), 0.04, 0.105,
                                      transform=axB.transAxes,
                                      facecolor=info["color"], edgecolor="white",
                                      linewidth=0.6, clip_on=False))
    axB.text(0.085, y - 0.005, grp, transform=axB.transAxes,
             fontsize=8.5, fontweight="bold", color="#222")
    axB.text(0.085, y - 0.075,
             ", ".join(info["members"]) + f"   (n = {len(info['members'])})",
             transform=axB.transAxes, fontsize=8, color="#444",
             family="monospace")

axB.text(0.02, 0.99, "B   22 intersecting hub targets by functional class",
         transform=axB.transAxes, fontsize=9, fontweight="bold")

# Footnote (CORRECTED 179 → 168)
fn = ("Targets standardised via UniProt (Homo sapiens). Drug targets are the union "
      "of SwissTargetPrediction (probability ≥ 0.10) and ChEMBL (pChEMBL ≥ 5) hits "
      "for the 168 core active ingredients (terpenoids, flavonoids, coumarins, "
      "alkaloids, phenolics/stilbenes; Supplementary Table S3). Disease targets are "
      "the deduplicated union of DisGeNET (score ≥ 0.10), OMIM and GeneCards "
      "(relevance ≥ 1.0) for the query \"herpes zoster\".")
axB.text(0.02, -0.02, fn, transform=axB.transAxes, fontsize=6.5, color="#555",
         wrap=True, va="top", style="italic")

plt.tight_layout()
save_fig(fig, "Figure6_venn")
print(f"Hub targets: {len(all_22)}; categories: {len(HUB_GROUPS)}")
print("Group sizes:", {k: len(v["members"]) for k, v in HUB_GROUPS.items()})
