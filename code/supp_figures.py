"""
Supplementary Figures S1–S4.
S1 — Database annotation coverage by ion mode
S2 — Chemotype-class breakdown of the 168-compound core library
S3 — VIP × |log2FC| scatter of the 1,382 differential metabolites
S4 — Sub-class clustering heatmap of top-50 differential metabolites
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list

DATA = HERE.parent / "data"
set_style()

# ---------- S1 ----------
fig, ax = plt.subplots(figsize=(6.5, 4.5))
modes = ["ESI$^{+}$", "ESI$^{-}$"]
total = [2611, 1808]; kegg = [733, 389]; hmdb = [1146, 656]; lipid = [373, 275]
x = np.arange(len(modes)); w = 0.20
ax.bar(x - 1.5*w, total, w, color="#3B8FC9", label="Total identified", edgecolor="white")
ax.bar(x - 0.5*w, kegg,  w, color="#3CA46E", label="KEGG annotated",  edgecolor="white")
ax.bar(x + 0.5*w, hmdb,  w, color="#E69F00", label="HMDB annotated",  edgecolor="white")
ax.bar(x + 1.5*w, lipid, w, color="#7A1F75", label="LIPID-MAPS annotated", edgecolor="white")
for k, vals in enumerate([total, kegg, hmdb, lipid]):
    for i, v in enumerate(vals):
        ax.text(x[i] + (k-1.5)*w, v + 60, f"{v:,}", ha="center", va="bottom",
                fontsize=7.5, fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(modes, fontsize=10)
ax.set_ylabel("Number of metabolites")
ax.set_title("Figure S1.   Database annotation coverage by ion mode",
             loc="left", fontweight="bold")
ax.legend(loc="upper right", fontsize=8)
ax.set_ylim(0, max(total) * 1.18)
save_fig(fig, "FigureS1_db_annotation")

# ---------- S2 ----------
# Counts derived directly from the curated source library
# (Cp project results/Cyathula_Prostrata_Active_Ingredients_Library.xlsx,
# 入库类别 = "核心类"/"核心库", 168 entries total).
fig, ax = plt.subplots(figsize=(6.0, 4.5))
classes = ["Terpenoid", "Flavonoid", "Coumarin"]
counts = [84, 68, 16]
colors_s2 = ["#2CA02C", "#FF7F0E", "#9467BD"]
bars = ax.bar(classes, counts, color=colors_s2, edgecolor="white", linewidth=0.6)
for b, v in zip(bars, counts):
    pct = 100 * v / sum(counts)
    ax.text(b.get_x() + b.get_width()/2, v + 1, f"{v}  ({pct:.1f}%)",
            ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.set_ylabel("Number of compounds")
ax.set_title("Figure S2.   Chemotype-class breakdown of the 168-compound core active-ingredient library",
             loc="left", fontweight="bold", fontsize=9)
ax.set_ylim(0, max(counts) * 1.25)
plt.setp(ax.get_xticklabels(), fontsize=9.5)
fig.subplots_adjust(bottom=0.13)
save_fig(fig, "FigureS2_class_breakdown")

# ---------- S3 ----------
df = pd.read_csv(DATA / "full_metabolites.csv")
df["log2FC"] = df["log2(Y6-13/J6-13)"].astype(float)
df["VIP"] = df["VIP"].astype(float)
df["super"] = df["ClassI"].fillna("Other").replace({"_":"Other"})
zh_to_en = {
    "脂质和类脂分子": "Lipids & lipid-like",
    "苯丙素类和聚酮类": "Phenylpropanoids & polyketides",
    "有机杂环化合物": "Organoheterocyclics",
    "有机酸及其衍生物": "Organic acids & deriv.",
    "有机氧化合物": "Organic oxygen cmp.",
    "苯类化合物": "Benzenoids",
    "生物碱及其衍生物": "Alkaloids",
    "Other": "Other / Unclassified",
}
df["super_en"] = df["super"].map(lambda s: zh_to_en.get(s, "Other / Unclassified"))
top_classes = df["super_en"].value_counts().head(7).index.tolist()
df["class_for_plot"] = df["super_en"].where(df["super_en"].isin(top_classes), "Other")

fig, ax = plt.subplots(figsize=(8.5, 5.6))
# HMDB super-class annotation is not resolved for this feature set, so points are
# shown in a single colour rather than a degenerate one-category legend.
ax.scatter(np.abs(df["log2FC"]), df["VIP"], s=8, color="#4C78A8",
           alpha=0.5, edgecolor="none")
ax.axvline(1.0, color="gray", lw=0.5, ls="--", alpha=0.7)
ax.axhline(1.0, color="gray", lw=0.5, ls="--", alpha=0.7)
ax.set_xlabel("|$\log_2$ fold change|  (Y6-13 / J6-13)")
ax.set_ylabel("PLS-DA VIP score")
ax.set_title("Figure S3.   VIP × |$\log_2$FC| scatter of the 1,382 differential metabolites",
             loc="left", fontweight="bold")
save_fig(fig, "FigureS3_VIP_log2FC")

# ---------- S4 ----------
top50 = df.assign(score=df["VIP"]*np.abs(df["log2FC"])).nlargest(50, "score").copy()
samples = ["J6-13-1","J6-13-2","J6-13-3","Y6-13-1","Y6-13-2","Y6-13-3"]
abund = top50[samples].astype(float).values
log_ab = np.log10(abund + 1)
z = (log_ab - log_ab.mean(axis=1, keepdims=True)) / log_ab.std(axis=1, keepdims=True).clip(min=1e-9)

# Cluster rows
order = leaves_list(linkage(z, method="ward"))

# Build label combining metabolite name + ClassII (HMDB Class)
labels = []
for j in order:
    n = top50["Name"].iloc[j]; n = str(n) if pd.notna(n) else f"Compound {j}"
    cls = top50["ClassII"].iloc[j]; cls = str(cls) if pd.notna(cls) else ""
    label = (n if len(n) <= 22 else n[:20]+"…") + ("  ["+cls[:18]+"]" if cls and cls != "_" else "")
    labels.append(label)

fig, ax = plt.subplots(figsize=(7.5, 9.5))
sns.heatmap(z[order], ax=ax, cmap="RdBu_r", center=0, vmin=-2.2, vmax=2.2,
            cbar_kws=dict(label="Row Z-score", shrink=0.45, pad=0.02),
            xticklabels=samples, yticklabels=labels, linewidths=0)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=5.4)
ax.set_title("Figure S4.   Sub-class-annotated heatmap of top-50 differential metabolites",
             loc="left", fontweight="bold")

save_fig(fig, "FigureS4_subclass_heatmap")
print("Supplementary figures S1-S4 written.")
