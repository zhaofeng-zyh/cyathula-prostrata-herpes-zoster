"""
Project-wide matplotlib/seaborn style for Cyathula prostrata manuscript figures.
Compliant with Phytomedicine (Elsevier) submission specs:
  - Vector PDF (preferred) + ≥ 600 dpi PNG
  - Sans-serif (Arial/Helvetica family); fallback DejaVu Sans
  - Colorblind-safe categorical palette (Okabe-Ito) for chemotype/group coding
  - Sequential palette: viridis (perceptually uniform, print-safe)
  - Figure widths in inches: single column 89 mm = 3.5 in; double column 183 mm = 7.2 in

Use as:
    from _style import set_style, save_fig, COLORS_CHEMOTYPE, COLORS_OKABE
    set_style()
    ...
    save_fig(fig, "Figure1_classification_correlation")
"""
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt

# --- Color palettes -----------------------------------------------------------
# Okabe-Ito 8-color colorblind-safe categorical palette
COLORS_OKABE = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
]

# Chemotype-class palette (5 classes used throughout the paper)
COLORS_CHEMOTYPE = {
    "Terpenoid":  "#1F77B4",  # blue
    "Flavonoid":  "#FF7F0E",  # orange
    "Coumarin":   "#9467BD",  # purple
    "Alkaloid":   "#2CA02C",  # green
    "Phenolic":   "#8C564B",  # brown
    "Stilbene":   "#8C564B",
}

# Group palette (Stem J6-13 = blue; Leaf Y6-13 = red, kept consistent across all figures)
COLOR_STEM = "#3B8FC9"  # Phytomedicine-style blue
COLOR_LEAF = "#C03A2B"  # Phytomedicine-style red

# Pathway-category palette
COLORS_PATHWAY = {
    "Inflammation/Immunity": "#C03A2B",
    "Antiviral":             "#7B4F9A",
    "Signal transduction":   "#1F77B4",
    "Neuro/Pain":            "#2CA02C",
    "Apoptosis":             "#FF7F0E",
}


def set_style():
    """Set Phytomedicine-grade matplotlib style. Call once at top of every figure script."""
    mpl.rcParams.update({
        "figure.dpi":           150,            # screen draft
        "savefig.dpi":          600,            # publication
        "font.family":          "sans-serif",
        "font.sans-serif":      ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":            8,
        "axes.titlesize":       9,
        "axes.titleweight":     "bold",
        "axes.labelsize":       8,
        "axes.labelweight":     "bold",
        "axes.linewidth":       0.8,
        "axes.spines.top":      False,
        "axes.spines.right":    False,
        "xtick.labelsize":      7,
        "ytick.labelsize":      7,
        "xtick.major.width":    0.8,
        "ytick.major.width":    0.8,
        "xtick.major.size":     3.0,
        "ytick.major.size":     3.0,
        "legend.fontsize":      7,
        "legend.frameon":       False,
        "legend.handlelength":  1.4,
        "lines.linewidth":      1.2,
        "patch.linewidth":      0.6,
        "savefig.bbox":         "tight",
        "savefig.pad_inches":   0.05,
        "pdf.fonttype":         42,             # editable text in PDF (TrueType, not Type 3)
        "ps.fonttype":          42,
        "axes.prop_cycle":      mpl.cycler(color=COLORS_OKABE),
    })


# --- Save helpers -------------------------------------------------------------
FIG_OUT = Path(__file__).resolve().parent.parent / "figures"


def save_fig(fig, stem, also_pdf=True, also_png=True, dpi_png=600):
    """Save fig as PDF and high-DPI PNG into the SUBMISSION_PACKAGE/figures dir."""
    FIG_OUT.mkdir(parents=True, exist_ok=True)
    if also_pdf:
        fig.savefig(FIG_OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.05)
    if also_png:
        fig.savefig(FIG_OUT / f"{stem}.png", bbox_inches="tight", pad_inches=0.05, dpi=dpi_png)


def figsize(cols="single", height_in=4.0):
    """Phytomedicine width helper. cols: 'single' (89 mm), 'double' (183 mm), or float."""
    if cols == "single":
        w = 3.5
    elif cols == "double":
        w = 7.2
    elif cols == "1.5":
        w = 5.5
    else:
        w = float(cols)
    return (w, height_in)
