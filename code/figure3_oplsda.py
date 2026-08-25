"""
Figure 3.  PLS-DA discrimination of stem (J6-13) vs leaf (Y6-13)
metabolomes (note: the original code implements PLS-DA, not true OPLS-DA;
manuscript text now reflects this).

A — model fit metrics (R²X, R²Y, Q²)
B — score plot (predictive vs orthogonal-projection)
C — Eriksson 200-permutation diagnostic with regression intercepts;
    panel annotates the R²Y intercept = +0.85 vs the < 0.40 validity
    threshold and explicitly flags that with n=3 vs n=3 the test cannot
    formally exclude overfitting (see manuscript §Limitations).
D — Loadings S-plot annotated with the four largest |p[1]| × |p_corr[1]|
    species (Lys-Asp-His, Trichagmalin B, Glycobismine G, Avermectin B1b
    monosaccharide).
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig, COLOR_STEM, COLOR_LEAF

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import LeaveOneOut

DATA = HERE.parent / "data"
set_style()
rng = np.random.default_rng(2026)

df = pd.read_csv(DATA / "full_metabolites.csv")
samples = ["J6-13-1","J6-13-2","J6-13-3","Y6-13-1","Y6-13-2","Y6-13-3"]
groups  = np.array(["Stem"]*3 + ["Leaf"]*3)
y = np.array([0]*3 + [1]*3, dtype=float)

X = np.log10(df[samples].T.values + 1)
X = StandardScaler().fit_transform(X)
y_s = (y - y.mean()) / y.std()

# --- Fit and CV ---
pls = PLSRegression(n_components=2, scale=False).fit(X, y_s.reshape(-1, 1))
r2y = pls.score(X, y_s.reshape(-1, 1))
r2x = (pls.x_scores_**2).sum() / (X**2).sum()
loo = LeaveOneOut()
preds = np.zeros_like(y_s)
for tr, te in loo.split(X):
    p = PLSRegression(n_components=2, scale=False).fit(X[tr], y_s[tr].reshape(-1, 1))
    preds[te] = p.predict(X[te]).ravel()
q2 = 1.0 - ((y_s - preds)**2).sum() / ((y_s - y_s.mean())**2).sum()

# Predictive vs orthogonal scores via OPLS-style deflation
def deflate(Xc, yc):
    w = (Xc.T @ yc).ravel(); w /= np.linalg.norm(w)
    tp = Xc @ w
    p  = Xc.T @ tp / (tp @ tp)
    w_o = p - (p @ w) * w
    if np.linalg.norm(w_o) > 0:
        w_o /= np.linalg.norm(w_o)
    to = Xc @ w_o
    return tp.ravel(), to.ravel(), w, p
tp, to, w, p = deflate(X - X.mean(0), y_s - y_s.mean())

# --- Permutation test ---
n_perm = 200  # 200 random draws; note only C(6,3)=20 unique labels exist (see box)
r2y_p, q2_p, sims = [], [], []
for i in range(n_perm):
    perm = rng.permutation(y_s)
    sims.append(np.dot(perm, y_s) / np.linalg.norm(perm) / np.linalg.norm(y_s))
    p_ = PLSRegression(n_components=2, scale=False).fit(X, perm.reshape(-1, 1))
    r2y_p.append(p_.score(X, perm.reshape(-1, 1)))
    pp = np.zeros_like(perm)
    for tr, te in loo.split(X):
        m = PLSRegression(n_components=2, scale=False).fit(X[tr], perm[tr].reshape(-1, 1))
        pp[te] = m.predict(X[te]).ravel()
    q2_p.append(1 - ((perm - pp)**2).sum() / ((perm - perm.mean())**2).sum())
sims = np.array(sims); r2y_p = np.array(r2y_p); q2_p = np.array(q2_p)

# Eriksson regression: y = a + b * |sim|; the 'real' point lies at sim=1
xs = np.concatenate([np.abs(sims), [1.0]])
yr = np.concatenate([r2y_p, [r2y]])
yq = np.concatenate([q2_p,  [q2]])
A = np.vstack([np.ones_like(xs), xs]).T
ar, br = np.linalg.lstsq(A, yr, rcond=None)[0]
aq, bq = np.linalg.lstsq(A, yq, rcond=None)[0]

# Empirical p
emp_r2y = (r2y_p >= r2y).mean()
emp_q2  = (q2_p  >= q2 ).mean()

# --- S-plot loadings ---
# correlation of each metabolite's deflated abundance with the predictive component
Xc = X - X.mean(0)
tp_norm = (tp - tp.mean()) / tp.std()
p_corr = np.array([np.corrcoef(Xc[:, j], tp_norm)[0, 1] for j in range(Xc.shape[1])])

# names
names = df["Name"].fillna("").astype(str).values
# importance metric for annotation
score = np.abs(p) * np.abs(p_corr)
idx_top = np.argsort(score)[-6:]  # top 6 to annotate

# --- Build figure ---
fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.6))
axA, axB = axes[0]
axC, axD = axes[1]

# Panel A — fit metrics
metrics = ["R²X (cum)", "R²Y (cum)", "Q²"]
vals = [r2x, r2y, q2]
bar_colors = ["#4FA9D5", "#5DBD7B", "#E07A5F"]
bars = axA.bar(metrics, vals, color=bar_colors, edgecolor="white", linewidth=1.2)
for b, v in zip(bars, vals):
    axA.text(b.get_x() + b.get_width()/2, v + 0.025, f"{v:.3f}",
             ha="center", va="bottom", fontsize=9, fontweight="bold")
axA.set_ylim(0, 1.15)
axA.set_ylabel("Cumulative variance / predictive ability")
axA.set_title("A   Model fit", loc="left", fontweight="bold")
axA.grid(axis="y", lw=0.4, ls="--", alpha=0.4)

# Panel B — score plot (predictive vs orthogonal)
for i in range(6):
    c = COLOR_STEM if y[i] == 0 else COLOR_LEAF
    axB.scatter(tp[i], to[i], s=110, color=c, edgecolor="white", linewidth=1.4, zorder=3)
    axB.annotate(samples[i], (tp[i], to[i]), xytext=(4, 4),
                 textcoords="offset points", fontsize=6.8, color="#333")
axB.axhline(0, color="gray", lw=0.4, ls="--", alpha=0.5)
axB.axvline(0, color="gray", lw=0.4, ls="--", alpha=0.5)
axB.set_xlabel("Predictive component  t[1]")
axB.set_ylabel("Orthogonal projection  t(o,1)")
axB.set_title("B   PLS-DA score plot", loc="left", fontweight="bold")
# Inline group legend at top right
axB.scatter([], [], s=80, color=COLOR_STEM, label="Stem (J6-13)")
axB.scatter([], [], s=80, color=COLOR_LEAF, label="Leaf (Y6-13)")
axB.legend(loc="upper right", fontsize=7)

# Panel C — permutation test (Eriksson regression intercepts)
axC.scatter(np.abs(sims), r2y_p, s=14, color="#5DBD7B", alpha=0.55, label="R²Y (perm)")
axC.scatter(np.abs(sims), q2_p,  s=14, color="#E07A5F", alpha=0.55, label="Q² (perm)")
axC.scatter([1.0], [r2y], marker="*", s=160, color="#1f3a78",
            edgecolor="white", lw=0.8, zorder=4, label="R²Y (real)")
axC.scatter([1.0], [q2],  marker="D", s=80,  color="#7a1f1f",
            edgecolor="white", lw=0.8, zorder=4, label="Q² (real)")
xline = np.linspace(0, 1, 50)
axC.plot(xline, ar + br*xline, color="#5DBD7B", lw=1.2)
axC.plot(xline, aq + bq*xline, color="#E07A5F", lw=1.2)
axC.axhline(0.40, color="gray", lw=0.6, ls=":")
axC.text(0.05, 0.42, "Eriksson R²Y validity threshold (0.40)",
         fontsize=6.5, color="gray", style="italic")
axC.set_xlim(-0.05, 1.08); axC.set_ylim(-1.15, 1.15)
axC.set_xlabel("|Correlation between permuted and real Y|")
axC.set_ylabel("R²Y / Q²")
axC.set_title("C   200-permutation validation", loc="left", fontweight="bold")
axC.legend(loc="lower right", fontsize=6.5, ncol=2)

# Annotation box with intercept values + honest interpretation
note = (f"R²Y intercept = {ar:+.2f}  (rule: < 0.40)\n"
        f"Q²  intercept = {aq:+.2f}  (rule: < 0.05)\n"
        f"Empirical p (R²Y) = {emp_r2y:.3f}\n"
        f"Empirical p (Q²)  = {emp_q2:.3f}\n"
        f"With n=3 vs n=3 only ~20 unique label\n"
        f"permutations exist; permutation test\n"
        f"cannot formally exclude overfitting.")
# anchored in the empty mid-band (below the R²Y cloud, above the Q² line) so the
# box no longer overlaps the salmon Q² regression line at the lower left
axC.text(0.02, 0.33, note, fontsize=6.4, va="top",
         family="monospace", color="#222",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff7e6",
                   edgecolor="#999", linewidth=0.6))

# Panel D — S-plot
ax = axD
sig_mask = score > np.percentile(score, 20)
ax.scatter(p[~sig_mask], p_corr[~sig_mask], s=4, color="#cccccc", alpha=0.55, label="Other")
ax.scatter(p[sig_mask & (p > 0)], p_corr[sig_mask & (p > 0)], s=8, color=COLOR_LEAF,
           alpha=0.7, label="Up in leaf (Y6-13)")
ax.scatter(p[sig_mask & (p < 0)], p_corr[sig_mask & (p < 0)], s=8, color=COLOR_STEM,
           alpha=0.7, label="Up in stem (J6-13)")

# Annotate the most leaf-defining (top by p>0) and most stem-defining (top by p<0) species
# Use a small set, manually offset to avoid stacking
top_pos = np.argsort(score * (p > 0))[-3:][::-1]
top_neg = np.argsort(score * (p < 0))[-3:][::-1]

annot_specs = []
for k, j in enumerate(top_pos):
    annot_specs.append((j, 0.012, 0.55 - 0.20*k))   # right column
for k, j in enumerate(top_neg):
    annot_specs.append((j, -0.012, -0.55 + 0.20*k)) # left column

for j, xt, yt in annot_specs:
    n = names[j] if names[j] else f"Compound {j}"
    short = n if len(n) <= 26 else n[:24] + "…"
    ax.scatter(p[j], p_corr[j], s=22, facecolor="none",
               edgecolor="black", lw=0.8, zorder=5)
    ax.annotate(short, xy=(p[j], p_corr[j]), xytext=(xt, yt),
                fontsize=6.4, color="black",
                ha="center", va="center",
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.5,
                                connectionstyle="arc3,rad=0.1"))

ax.axhline(0, color="gray", lw=0.4, ls="--", alpha=0.5)
ax.axvline(0, color="gray", lw=0.4, ls="--", alpha=0.5)
ax.set_xlim(-0.04, 0.04); ax.set_ylim(-1.15, 1.15)
ax.set_xlabel("p[1]  (covariance loading)")
ax.set_ylabel("p_corr[1]  (correlation loading)")
ax.set_title("D   PLS-DA S-plot", loc="left", fontweight="bold")
ax.legend(loc="lower right", fontsize=6.5)

plt.tight_layout()
save_fig(fig, "Figure3_PLSDA")
print(f"R²X={r2x:.3f}  R²Y={r2y:.3f}  Q²={q2:.3f}")
print(f"R²Y intercept={ar:+.3f}  Q² intercept={aq:+.3f}  emp_p_R²Y={emp_r2y:.3f}  emp_p_Q²={emp_q2:.3f}")
print("Top S-plot annotations (idx, name):")
for j in idx_top[::-1]:
    print(f"  {j:5d}  {names[j][:50]:50s}  p={p[j]:+.3f}  p_corr={p_corr[j]:+.3f}")
