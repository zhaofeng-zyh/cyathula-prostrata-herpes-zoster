"""
Figure 12 (REAL DATA).  50-ns molecular-dynamics of the completed compound-target
complexes, computed from the actual GROMACS 2025.4 trajectories (amber99sb-ildn /
GAFF, TIP3P, 0.15 M NaCl; 2 fs; 50 ps sampling; PBC fixed via whole->nojump).

Panels:
  A  Backbone RMSD vs time
  B  Per-residue Cα RMSF
  C  Radius of gyration vs time
  D  Ligand RMSD vs time (protein-fit; pose stability in the pocket)
  E  Protein-ligand contacts (< 0.35 nm) vs time
  F  MM-GBSA binding free energy (ΔG and components)  [added once mmgbsa_real.csv exists]

Data source: SUBMISSION_PACKAGE/data/md_real/<system>__<metric>.xvg (real gmx output)
NOTE: replaces the earlier synthetic Figure12_MD (np.random) — this is real data.
Currently the 2 finished TNF-α systems; COX-2 / NF-κB systems append when their runs finish.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from _style import set_style, save_fig, COLOR_STEM, COLOR_LEAF

import numpy as np
import matplotlib.pyplot as plt

DATA = HERE.parent / "data" / "md_real"
set_style()

# --- systems (real completed runs); extend when COX-2 / NF-κB finish -----------
SYSTEMS = [
    ("DecursinolAng_TNFa", "Decursinol angelate · TNF-α", "#0072B2"),
    ("Carnosol_TNFa",      "Carnosol · TNF-α",            "#D55E00"),
    ("DecursinolAng_COX2", "Decursinol angelate · COX-2", "#009E73"),
    ("ArnicolideD_NFkB",   "Arnicolide D · NF-κB p65",   "#CC79A7"),
]

def read_xvg(path):
    xs, ys = [], []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line[0] in "@#&":
            continue
        p = line.split()
        try:
            xs.append(float(p[0])); ys.append(float(p[1]))
        except (ValueError, IndexError):
            continue
    return np.array(xs), np.array(ys)

def time_ns(x):
    """Normalize an x axis to ns (gmx writes ps unless -tu ns was passed)."""
    return x/1000.0 if x.max() > 1000 else x

def load(sysid, metric):
    return read_xvg(DATA / f"{sysid}__{metric}.xvg")

fig = plt.figure(figsize=(7.2, 8.4))
gs = fig.add_gridspec(3, 2, hspace=0.48, wspace=0.30,
                      top=0.95, bottom=0.06, left=0.10, right=0.97)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])
axE = fig.add_subplot(gs[2, 0]); axF = fig.add_subplot(gs[2, 1])

summary = {}

# A — backbone RMSD (nm -> Å)
for sid, lab, c in SYSTEMS:
    t, r = load(sid, "rmsd_backbone"); t = time_ns(t); r = r*10.0
    axA.plot(t, r, color=c, lw=0.6, alpha=0.85, label=lab)
    summary.setdefault(lab, {})["RMSD_bb_mean_A"] = float(r[t>=20].mean())
    summary[lab]["RMSD_bb_max_A"] = float(r.max())
axA.set_xlabel("Time (ns)"); axA.set_ylabel("Backbone RMSD (Å)")
axA.set_title("A  Backbone RMSD", loc="left", pad=6)
axA.set_xlim(0, 50)
# opaque background so the high-RMSD trace never shows through the legend text
axA.legend(loc="upper left", fontsize=6.5, frameon=True, framealpha=0.95,
           facecolor="white", edgecolor="0.8")

# B — Cα RMSF (nm -> Å). Plot vs a continuous Cα index (proteins differ in length; the two
# TNF-α protomer chains share residue numbers, so a continuous index avoids a spurious jump-back).
for sid, lab, c in SYSTEMS:
    res, f = load(sid, "rmsf_calpha"); f = f*10.0
    axB.plot(np.arange(f.size), f, color=c, lw=0.7, alpha=0.85, label=lab)
    summary[lab]["RMSF_mean_A"] = float(f.mean())
axB.set_xlabel("Cα index"); axB.set_ylabel("Cα RMSF (Å)")
axB.set_title("B  Per-residue Cα RMSF", loc="left", pad=6)

# C — radius of gyration (nm)
for sid, lab, c in SYSTEMS:
    t, g = load(sid, "gyrate"); t = time_ns(t)
    axC.plot(t, g, color=c, lw=0.6, alpha=0.85, label=lab)
    summary[lab]["Rg_mean_nm"] = float(g[t>=20].mean())
axC.set_xlabel("Time (ns)"); axC.set_ylabel("R$_g$ (nm)")
axC.set_title("C  Radius of gyration", loc="left", pad=6)
axC.set_xlim(0, 50)

# D — minimum protein-ligand distance (nm). Robust (minimum-image); shows the ligand never
# leaves the pocket. (Ligand RMSD via gmx rms was PBC-imaging-unreliable for the buried COX-2
# pose, so the min-distance is used as the pose-stability metric alongside the contact count.)
for sid, lab, c in SYSTEMS:
    t, d = load(sid, "mindist"); t = time_ns(t); d = d*10.0
    axD.plot(t, d, color=c, lw=0.6, alpha=0.85, label=lab)
    summary[lab]["mindist_mean_A"] = float(d[t>=20].mean())
axD.set_xlabel("Time (ns)"); axD.set_ylabel("Min. protein–ligand dist. (Å)")
axD.set_title("D  Minimum protein–ligand distance", loc="left", pad=6)
axD.set_xlim(0, 50); axD.set_ylim(0, None)

# E — protein-ligand contacts < 0.35 nm
for sid, lab, c in SYSTEMS:
    t, n = load(sid, "ncontacts"); t = time_ns(t)
    axE.plot(t, n, color=c, lw=0.5, alpha=0.8, label=lab)
    summary[lab]["contacts_mean"] = float(n[t>=20].mean())
    tm, dm = load(sid, "mindist"); dm = time_ns(tm)  # noqa
    summary[lab]["mindist_mean_nm"] = float(load(sid, "mindist")[1][ time_ns(load(sid,"mindist")[0])>=20 ].mean())
axE.set_xlabel("Time (ns)"); axE.set_ylabel("Contacts < 0.35 nm")
axE.set_title("E  Protein–ligand contacts", loc="left", pad=6)
axE.set_xlim(0, 50)

# F — MM-GBSA (if available)
mmg = HERE.parent / "data" / "mmgbsa_real.csv"
if mmg.exists():
    import csv
    rows = list(csv.DictReader(open(mmg)))
    labels = [r["System"] for r in rows]
    comps = ["VDW", "EEL", "EGB", "ESURF"]
    comp_lab = ["van der Waals", "Electrostatic", "Polar solv.", "Nonpolar solv."]
    colors_d = ["#3B8FC9", "#E07A5F", "#5DBD7B", "#9467BD"]
    x = np.arange(len(labels)); w = 0.2
    for k,(comp,cl) in enumerate(zip(comps, colors_d)):
        axF.bar(x+(k-1.5)*w, [float(r[comp]) for r in rows], w, color=cl,
                edgecolor="white", lw=0.5, label=comp_lab[k])
    axF.axhline(0, color="black", lw=0.5)
    axF.set_xticks(x); axF.set_xticklabels([l.replace(" · ","\n· ") for l in labels], fontsize=6.5)
    axF.set_ylabel("Energy (kcal mol$^{-1}$)")
    axF.set_title("F  MM-GBSA / MM-PBSA decomposition", loc="left", pad=6)
    axF.legend(loc="lower left", fontsize=6)
    ytop = axF.get_ylim()[1]
    for i,r in enumerate(rows):
        pb = f"\nΔG$_{{PB}}$={float(r['DeltaTOTAL_PB']):+.1f}±{float(r['SEM_PB']):.1f}" if r.get('DeltaTOTAL_PB') else ""
        axF.text(x[i], ytop*0.98, f"ΔG$_{{GB}}$={float(r['DeltaTOTAL']):+.1f}±{float(r['SEM']):.1f}{pb}",
                 ha="center", va="top", fontsize=6, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#888", lw=0.5))
else:
    axF.text(0.5, 0.5, "MM-GBSA panel\n(pending mmgbsa_real.csv)", ha="center", va="center",
             fontsize=8, color="gray", style="italic", transform=axF.transAxes)
    axF.set_title("F  MM-GBSA decomposition", loc="left", pad=6)
    axF.set_xticks([]); axF.set_yticks([])

save_fig(fig, "Figure12_MD")

# write summary CSV for the manuscript / Table
import csv
outcsv = HERE.parent / "data" / "md_real_summary.csv"
with open(outcsv, "w", newline="") as fh:
    wtr = csv.writer(fh)
    wtr.writerow(["System","RMSD_bb_mean_A(>=20ns)","RMSD_bb_max_A","RMSF_mean_A",
                  "Rg_mean_nm(>=20ns)","contacts_mean(>=20ns)","mindist_mean_nm(>=20ns)"])
    for lab in summary:
        s = summary[lab]
        wtr.writerow([lab, f"{s['RMSD_bb_mean_A']:.2f}", f"{s['RMSD_bb_max_A']:.2f}",
                      f"{s['RMSF_mean_A']:.2f}", f"{s['Rg_mean_nm']:.3f}",
                      f"{s['contacts_mean']:.0f}", f"{s['mindist_mean_nm']:.3f}"])
print("=== REAL MD summary (>=20 ns equilibrated window) ===")
for lab in summary:
    print(lab, summary[lab])
