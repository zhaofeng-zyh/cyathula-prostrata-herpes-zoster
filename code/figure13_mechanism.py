"""
Figure 13.  Integrated mechanism schematic.
Polished from original — same content, cleaner layout, consistent palette.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
set_style()


def box(ax, x, y, w, h, text, fc, ec="black", fontsize=8, weight="bold"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02",
                        linewidth=1.0, facecolor=fc, edgecolor=ec)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fontsize, fontweight=weight)


def arrow(ax, xy1, xy2, color="#444", style="-|>", lw=1.4):
    a = FancyArrowPatch(xy1, xy2, arrowstyle=style, color=color,
                         lw=lw, mutation_scale=14)
    ax.add_patch(a)


fig, ax = plt.subplots(figsize=(11.5, 8.0))
ax.set_xlim(0, 12); ax.set_ylim(0, 9)
ax.set_axis_off()

# Title
ax.text(6, 8.6, "Integrated mechanism of " + r"$\mathit{C.\ prostrata}$" + " against herpes-zoster-associated inflammation",
        ha="center", va="center", fontsize=11, fontweight="bold")

# Top: VZV reactivation
box(ax, 4.5, 7.0, 3.0, 0.8,
    "Varicella-zoster virus reactivation\n→ herpes zoster + neuralgia",
    fc="#F8C8C0", fontsize=8.5)
# Three disease modules
box(ax, 0.4, 5.0, 3.2, 1.4,
    "Skin keratinocyte\n• damage-associated patterns\n• interferon response (IFNB1)",
    fc="#FFE7B0", fontsize=7.5)
box(ax, 4.1, 5.0, 3.8, 1.4,
    "Sensory neuron (DRG)\n• nociceptor sensitisation\n• TRPV1 / NMDA / Ca2+",
    fc="#D7C9E4", fontsize=7.5)
box(ax, 8.4, 5.0, 3.2, 1.4,
    "Immune cell infiltration\n• Th1 / CD4+ / CD8+\n• macrophage activation",
    fc="#FFE0DA", fontsize=7.5)

# Pro-inflammatory cascade
box(ax, 1.8, 3.0, 8.4, 1.2,
    "Pro-inflammatory cascade\n"
    "NF-κB → TNF-α / IL-6 / IL-1β / CXCL10 / CCL5\n"
    "COX-2 / PGE2 · iNOS · IL-17 axis",
    fc="#F4B7B0", fontsize=7.8)

# Cyathula leaves & bioactives
box(ax, 0.2, 0.5, 3.6, 1.8,
    r"$\mathit{Cyathula\ prostrata}$ leaves" + "\n(Y6-13 vs J6-13)\n"
    "Bioactive enrichment:\n"
    "• Sesquiterpene lactones (Arnicolide D, Dehydrocostus L.)\n"
    "• Flavonoids (Salvigenin, Hibiscitrin)\n"
    "• Coumarins (Decursinol angelate)\n"
    "• Stilbenes (Pterostilbene, Gigantol)",
    fc="#CDE7CC", fontsize=6.8)

# Multi-target panel
box(ax, 4.4, 0.6, 3.2, 1.5,
    "Multi-component / multi-target synergy\n22 hub targets across\n5 enriched pathways",
    fc="#FCE19A", fontsize=7.5)

# Three downstream effects
box(ax, 8.5, 1.7, 3.2, 0.9,
    "Conserved α-herpesvirus\nhost response (not VZV-specific)",
    fc="#FDE2C8", fontsize=7.0)
box(ax, 8.5, 0.8, 3.2, 0.9,
    "Anti-inflammatory effect\n(NF-κB ↓, COX-2 ↓, TNF/IL-6 ↓)",
    fc="#F8C5BF", fontsize=7.4)
box(ax, 8.5, -0.1, 3.2, 0.9,
    "Analgesia & neuroprotection\n(Ca2+/TRPV1 modulation)",
    fc="#F1D3DC", fontsize=7.4)

# Arrows
arrow(ax, (6, 7.0), (6, 6.45))
# fan-out arrows start at the DRG box's bottom edge (y=5.0), not inside it,
# so the connector lines no longer cross the box's own text
arrow(ax, (6, 4.98), (4.05, 4.28))
arrow(ax, (6, 4.98), (6.0, 4.28))
arrow(ax, (6, 4.98), (7.95, 4.28))
arrow(ax, (6, 3.0), (6, 2.15), color="#a33", lw=1.6)
ax.text(6.4, 2.55, "Inhibits", color="#a33", fontsize=8, fontweight="bold", style="italic")
arrow(ax, (3.8, 1.4), (4.4, 1.35), color="#196f1f", style="-|>", lw=1.5)
arrow(ax, (7.6, 1.35), (8.5, 2.15), color="#444")
arrow(ax, (7.6, 1.35), (8.5, 1.30), color="#444")
arrow(ax, (7.6, 1.35), (8.5, 0.4), color="#444")

save_fig(fig, "Figure13_mechanism")
print("Figure 13 mechanism schematic written.")
