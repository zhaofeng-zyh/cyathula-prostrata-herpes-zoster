"""
Figure S5.  Real pathway/GO enrichment (g:Profiler, g:SCS-corrected) of the predicted
targets of the lead/representative active ingredients (132-gene unbiased union;
data/target_enrichment.csv). Complements Fig. 8 by mapping the real compound targets
onto pathways — completing the compound -> target -> pathway mechanism chain.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

DATA = HERE.parent / "data"
set_style()
df = pd.read_csv(DATA / "target_enrichment.csv")
# keep interpretable terms: not too generic, not too tiny
df = df[(df.term_size >= 5) & (df.term_size <= 800)].copy()
# omit carbonic-anhydrase-driven off-target artefacts (frequent-hitter class)
_art = df.term_name.str.contains("nitrogen|carbon dioxide|bicarbonate", case=False, regex=True)
df = df[~_art].copy()
df["neglogp"] = -np.log10(df.p_value.clip(lower=1e-300))

SRC_COLOR = {"KEGG": "#C03A2B", "REAC": "#1F77B4", "GO:BP": "#2CA02C"}
SRC_LABEL = {"KEGG": "KEGG pathway", "REAC": "Reactome", "GO:BP": "GO biological process"}
picks = []
for src, n in [("KEGG", 9), ("REAC", 7), ("GO:BP", 5)]:
    picks.append(df[df.source == src].nsmallest(n, "p_value"))
sel = pd.concat(picks)
# order: by source group then by significance (so bars cluster by source)
sel["src_order"] = sel.source.map({"KEGG": 0, "REAC": 1, "GO:BP": 2})
sel = sel.sort_values(["src_order", "neglogp"], ascending=[False, True]).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(8.2, 8.6))
y = np.arange(len(sel))
ax.barh(y, sel.neglogp, color=[SRC_COLOR[s] for s in sel.source],
        edgecolor="white", linewidth=0.6, height=0.72)
for i, r in sel.iterrows():
    ax.text(r.neglogp + 0.3, i, f"{r.term_name[:46]}  ({int(r.intersection_size)})",
            va="center", ha="left", fontsize=7.6, color="#222")
ax.set_yticks([]); ax.set_xlabel("$-\log_{10}$(adjusted P)", fontweight="bold")
ax.set_xlim(0, sel.neglogp.max()*1.55)
ax.set_title("Pathway / GO enrichment of the lead-compound predicted targets\n(real g:Profiler analysis; 132-gene set)",
             fontsize=10.5, fontweight="bold", loc="left", pad=12)
handles = [mlines.Line2D([], [], color=SRC_COLOR[s], marker="s", linestyle="",
            markersize=9, label=SRC_LABEL[s]) for s in ["KEGG", "REAC", "GO:BP"]]
ax.legend(handles=handles, loc="upper right", fontsize=8, frameon=True, framealpha=0.9, edgecolor="none")
ax.text(0.0, -1.4, "Bar label = enriched term (no. of query genes in term). g:SCS correction, P < 0.05. "
        "Carbonic-anhydrase-driven metabolic terms (a frequent-hitter off-target class) omitted for clarity.",
        transform=ax.get_yaxis_transform(), fontsize=6.6, color="#666", style="italic")
save_fig(fig, "FigureS5_target_enrichment")
print(f"Enrichment figure: {len(sel)} terms plotted "
      f"(KEGG {sum(sel.source=='KEGG')}, REAC {sum(sel.source=='REAC')}, GO:BP {sum(sel.source=='GO:BP')})")
print("top terms:", list(sel.sort_values('p_value').term_name.head(6)))
