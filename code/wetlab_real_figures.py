"""
Generate the REAL wet-lab validation figures for the Phytomedicine manuscript
from the genuine analysis outputs (NOT the old template).

  Figure 14  RAW264.7 CCK-8 viability  <- wetlab/analysis/batch1/batch1_CCK8_summary.csv
  Figure 15  RAW264.7 TNF-a ELISA      <- wetlab/analysis/batch2_ELISA_0613/ELISA_0613_analysis.xlsx

Plate selection, exclusions and assay constants follow the Methods section of the
article. All numbers are read from the analysis files at run time; none are hard-coded.
The analysis files are not distributed in this repository (see README).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from _style import set_style, save_fig

set_style()

HERE = Path(__file__).resolve().parent
WET  = HERE.parent / "wetlab" / "analysis"

LEADS = ["Arnicolide D", "Salvigenin", "Decursinol angelate", "Carnosol"]
# sequential shade per concentration rank (light -> dark = low -> high dose)
CONC_SHADES = ["#B7D4EA", "#7FB0D6", "#3B7BB5", "#1F4E79"]


def _shade_map(concs):
    order = sorted(set(concs))
    cmap = {c: CONC_SHADES[min(i, len(CONC_SHADES) - 1)] for i, c in enumerate(order)}
    return cmap


# =====================================================================
# Figure 14  -- CCK-8 viability (2026-06-13 plates only)
# =====================================================================
def figure14():
    df = pd.read_csv(WET / "batch1" / "batch1_CCK8_summary.csv")
    df = df[df["plate"].str.contains("2026-06-13", na=False)].copy()
    df = df[df["compound"].isin(LEADS)].copy()
    df["conc_uM"] = df["conc_uM"].astype(float)
    df["viab_sd"] = df["sdOD"] / df["ref_untreatedOD"] * 100.0

    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    shade = _shade_map(df["conc_uM"])
    group_w, bar_w = 0.8, None
    xticks, xlabels = [], []
    x0 = 0
    for comp in LEADS:
        sub = df[df["compound"] == comp].sort_values("conc_uM")
        n = len(sub)
        bar_w = group_w / max(n, 1)
        xs = x0 + (np.arange(n) - (n - 1) / 2) * bar_w
        for x, (_, r) in zip(xs, sub.iterrows()):
            below = r["viability_pct"] < 80
            ax.bar(x, r["viability_pct"], bar_w * 0.92,
                   color=shade[r["conc_uM"]],
                   edgecolor="#C0392B" if below else "#333333",
                   linewidth=1.3 if below else 0.5,
                   hatch="//" if below else None, zorder=3)
            ax.errorbar(x, r["viability_pct"], yerr=r["viab_sd"], fmt="none",
                        ecolor="#333333", elinewidth=0.7, capsize=1.8, zorder=4)
            ax.text(x, min(r["viability_pct"] + r["viab_sd"] + 4, 244), f"{r['conc_uM']:g}",
                    ha="center", va="bottom", fontsize=5.2, color="#555555")
        xticks.append(x0)
        xlabels.append(comp.replace(" ", "\n", 1))
        x0 += 1

    ax.axhline(80, ls="--", lw=1.0, color="#C0392B", zorder=2)
    ax.text(x0 - 0.55, 82, "80% viability threshold", fontsize=6, color="#C0392B", va="bottom")
    ax.axhline(100, ls=":", lw=0.7, color="#888888", zorder=1)
    ax.set_xticks(xticks); ax.set_xticklabels(xlabels)
    ax.set_ylabel("Cell viability (% of untreated)")
    ax.set_ylim(0, 290)   # extra headroom so the legend clears the tallest error bar
    ax.set_yticks([0, 50, 80, 100, 150, 200, 250])
    ax.set_title("RAW 264.7 viability (CCK-8, 18 h; mean ± SD, n = 3)")
    handles = [Patch(facecolor=CONC_SHADES[i], edgecolor="#333333", label=lbl)
               for i, lbl in enumerate(["low", "", "", "high"]) if lbl]
    handles = [Patch(facecolor="#7FB0D6", edgecolor="#333333", label="dose (µM), light→dark = low→high"),
               Patch(facecolor="white", edgecolor="#C0392B", hatch="//", label="< 80% (cytotoxic)")]
    ax.legend(handles=handles, loc="upper right", fontsize=6,
              frameon=True, framealpha=1.0, facecolor="white", edgecolor="0.85")
    fig.tight_layout()
    save_fig(fig, "Figure14_CCK8_viability")
    plt.close(fig)
    print("Figure14_CCK8_viability written")


# =====================================================================
# Figure 15  -- TNF-a ELISA inhibition (DMSO+LPS = 100% response)
# =====================================================================
def figure15():
    raw = pd.read_excel(WET / "batch2_ELISA_0613" / "ELISA_0613_analysis.xlsx",
                        sheet_name="TNF-a", header=None)
    # header=None column positions: 0=Compound 1=Conc 2=meanOD 3=sdOD 4=n 5=inhibition% 6=pg/mL 7=note
    # rows: 1=controls, 3=sub-header, 4..21=data, 8=in-table BAY 11-7082
    od_blank = float(raw.iloc[1, 2]); od_dmso = float(raw.iloc[1, 6])
    denom = od_dmso - od_blank
    tbl = raw.iloc[4:22, [0, 1, 2, 3, 5, 7]].copy()
    tbl.columns = ["compound", "conc", "meanOD", "sdOD", "inhib", "note"]
    tbl = tbl[tbl["compound"].isin(LEADS)].copy()
    for c in ["conc", "meanOD", "sdOD", "inhib"]:
        tbl[c] = pd.to_numeric(tbl[c], errors="coerce")
    tbl = tbl.dropna(subset=["conc", "inhib"])
    tbl["inhib_sd"] = tbl["sdOD"] / denom * 100.0
    tbl["cytotox"] = tbl["note"].astype(str).str.contains("CYTOTOX", na=False)

    bay = float(raw.iloc[8, 5])   # BAY 11-7082 (5 µM) inhibition%

    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    shade = _shade_map(tbl["conc"])
    group_w = 0.8
    xticks, xlabels = [], []
    x0 = 0
    for comp in LEADS:
        sub = tbl[tbl["compound"] == comp].sort_values("conc")
        n = len(sub)
        bar_w = group_w / max(n, 1)
        xs = x0 + (np.arange(n) - (n - 1) / 2) * bar_w
        for x, (_, r) in zip(xs, sub.iterrows()):
            cyto = bool(r["cytotox"])
            ax.bar(x, r["inhib"], bar_w * 0.92,
                   color="#BBBBBB" if cyto else shade[r["conc"]],
                   edgecolor="#333333", linewidth=0.5,
                   hatch="xx" if cyto else None, zorder=3)
            ax.errorbar(x, r["inhib"], yerr=r["inhib_sd"], fmt="none",
                        ecolor="#333333", elinewidth=0.7, capsize=1.8, zorder=4)
            ax.text(x, min(max(r["inhib"], 0) + r["inhib_sd"] + 3, 126), f"{r['conc']:g}",
                    ha="center", va="bottom", fontsize=5.2, color="#555555")
        xticks.append(x0); xlabels.append(comp.replace(" ", "\n", 1))
        x0 += 1

    ax.axhline(bay, ls="--", lw=1.0, color="#2E7D32", zorder=2)
    ax.text(x0 - 0.05, bay - 6, f"BAY 11-7082\n(5 µM) = {bay:.0f}%",
            fontsize=6, color="#2E7D32", ha="right", va="top")
    ax.axhline(0, lw=0.8, color="#333333", zorder=2)
    ax.set_xticks(xticks); ax.set_xticklabels(xlabels)
    ax.set_ylabel("TNF-α inhibition (%, DMSO+LPS = 100% response)")
    ax.set_ylim(-75, 130)
    ax.set_title("Inhibition of LPS-induced TNF-α (RAW 264.7 ELISA; mean ± SD, n = 3)")
    handles = [Patch(facecolor="#7FB0D6", edgecolor="#333333", label="dose (µM), light→dark = low→high"),
               Patch(facecolor="#BBBBBB", edgecolor="#333333", hatch="xx",
                     label="cytotoxic dose (CCK-8 < 80%), excluded")]
    ax.legend(handles=handles, loc="lower left", fontsize=6)
    fig.tight_layout()
    save_fig(fig, "Figure15_TNFa_ELISA")
    plt.close(fig)
    print(f"Figure15_TNFa_ELISA written (denom={denom:.4f}, BAY={bay:.1f}%)")


if __name__ == "__main__":
    figure14()
    figure15()
