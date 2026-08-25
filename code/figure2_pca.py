"""
Figure 2.  PCA score plots for ESI− and ESI+ on log10-transformed,
mean-centred and unit-variance-scaled abundance of the 1,382 differential
metabolites.  Sample-label collisions in the original figure resolved
with adjustText.

Caveat: with n=3 vs n=3, "PC1 captures 83.7% / 84.8% of variance" is a
direct re-computation from the supplied differential-metabolite matrix.
We do not run a separate PCA per ion mode because the supplied
data/full_metabolites.csv combines both modes; we therefore split by the
'IonMode' column ('P' = ESI+, 'N' = ESI-).
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
from matplotlib.patches import Ellipse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

try:
    from adjustText import adjust_text
    HAS_ADJ = True
except ImportError:
    HAS_ADJ = False

DATA = HERE.parent / "data"
set_style()

df = pd.read_csv(DATA / "full_metabolites.csv")
samples = ["J6-13-1", "J6-13-2", "J6-13-3", "Y6-13-1", "Y6-13-2", "Y6-13-3"]
groups  = ["Stem"] * 3 + ["Leaf"] * 3
colors  = [COLOR_STEM] * 3 + [COLOR_LEAF] * 3


def confidence_ellipse(x, y, ax, n_std=2.0, **kwargs):
    """2-σ confidence ellipse for a 2D scatter."""
    if x.size < 2:
        return None
    cov = np.cov(x, y)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    w, h = 2.0 * n_std * np.sqrt(np.maximum(eigvals, 1e-12))
    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
    ell = Ellipse((np.mean(x), np.mean(y)), width=w, height=h, angle=angle, **kwargs)
    ax.add_patch(ell)
    return ell


def pca_panel(ax, mode_label, mode_filter):
    sub = df[df["IonMode"] == mode_filter]
    X = np.log10(sub[samples].T.values + 1)
    X = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2).fit(X)
    scores = pca.transform(X)
    var = pca.explained_variance_ratio_ * 100

    for i, (name, c, g) in enumerate(zip(samples, colors, groups)):
        ax.scatter(scores[i, 0], scores[i, 1], s=110, color=c,
                   edgecolor="white", linewidth=1.4, zorder=3)
    # Confidence ellipses
    for g, c in [("Stem", COLOR_STEM), ("Leaf", COLOR_LEAF)]:
        idx = [i for i, gi in enumerate(groups) if gi == g]
        confidence_ellipse(scores[idx, 0], scores[idx, 1], ax, n_std=2.0,
                           facecolor=c, alpha=0.12, edgecolor=c, linewidth=1.0)

    # Labels with collision avoidance
    texts = []
    for i, name in enumerate(samples):
        t = ax.text(scores[i, 0], scores[i, 1], name, fontsize=6.8, color="#333",
                    ha="center", va="center", zorder=4)
        texts.append(t)
    if HAS_ADJ:
        adjust_text(texts, ax=ax,
                    expand_text=(1.3, 1.5), expand_points=(1.5, 1.7),
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, alpha=0.7))

    ax.axhline(0, color="gray", lw=0.4, ls="--", alpha=0.5)
    ax.axvline(0, color="gray", lw=0.4, ls="--", alpha=0.5)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title(mode_label, loc="left", fontsize=9, fontweight="bold")
    ax.grid(False)
    return var


fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.2))
v1 = pca_panel(axes[0], "ESI$^{-}$ mode", "N")
v2 = pca_panel(axes[1], "ESI$^{+}$ mode", "P")

# Shared legend
handles = [mpatches.Patch(facecolor=COLOR_STEM, edgecolor="none", label="Stem (J6-13)"),
           mpatches.Patch(facecolor=COLOR_LEAF, edgecolor="none", label="Leaf (Y6-13)")]
fig.legend(handles=handles, loc="upper center", ncol=2,
           bbox_to_anchor=(0.5, 1.00), frameon=False, fontsize=8)

plt.tight_layout(rect=[0, 0, 1, 0.94])
save_fig(fig, "Figure2_PCA")
print(f"PCA done. ESI- variance explained: PC1={v1[0]:.1f}% PC2={v1[1]:.1f}%; ESI+ PC1={v2[0]:.1f}% PC2={v2[1]:.1f}%")
