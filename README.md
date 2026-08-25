# *Cyathula prostrata* — herpes-zoster chemotype study: analysis code

Analysis and figure-generation code accompanying the manuscript on the
*Cyathula prostrata* chemotype–target–pathway atlas.

---

## ⚠️ Scope of this repository (changed 2026-08-25)

**This repository now contains source code only. It no longer distributes the
processed datasets.**

Earlier versions of this repository (and the v1.0.0 archive) also contained
processed metabolomics, network, docking, molecular-dynamics and in-vitro tables.
Those data files have been removed from this repository at the authors' request
while intellectual-property matters relating to the work are being resolved.

**Consequently the scripts here will not run end-to-end as published.** They are
provided so that the analysis and plotting logic can be inspected, not so that the
figures can be regenerated from this repository alone.

### What is not here

- `data/` — processed metabolomics, network/enrichment, machine-learning, docking
  and molecular-dynamics tables
- `docking_real/` — docking score table and GROMACS run-parameter files
- `wetlab/` — processed in-vitro (CCK-8 / ELISA) results
- `figures/` — rendered figure panels and docking-pose composites
- Interpretive captions and numerical result footnotes that were previously
  hard-coded inside `code/build_supp_tables.py` and `code/wetlab_real_figures.py`.
  Where such a string was removed, the code now carries a bracketed placeholder and
  a comment marking the removal.

### What is here

- `code/` — the Python analysis and figure scripts, unchanged in structure and logic
- `CITATION.cff`, `.zenodo.json` — citation metadata

### No simulated data has been substituted

We deliberately did **not** replace the removed tables with synthetic or
proportionally-scaled surrogate values. Surrogate data that preserved the real
proportions would still convey the study's quantitative findings, and surrogate data
that did not would risk being mistaken for the real results if the pipeline were run.
Rather than publish numbers that are neither real nor obviously fake, the data are
simply absent and this notice states so.

### Obtaining the data

The processed data are archived, with the code as it stood at v1.0.0, at:

- **10.5281/zenodo.21231405** — analysis code and processed data
- **10.5281/zenodo.21231656** — molecular-dynamics trajectories (GROMACS, 4 × 50 ns)

Raw instrument data (metabolomics mzML) are being deposited separately; see the
Data Availability statement of the article. For access to material not covered by
those deposits, contact the corresponding author.

### Note on the v1.0.0 tag

The `v1.0.0` tag has been moved to this code-only commit so that references of the
form `.../tree/v1.0.0` continue to resolve. The Zenodo archive under
10.5281/zenodo.21231405 preserves the original v1.0.0 tree and is unaffected by this
change; where the two differ, **the Zenodo archive is the record of what was
originally released**.

---

## Software environment

Python 3.11 with `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`,
`networkx`, `adjustText`, `openpyxl`.

```bash
conda create -n cprostrata python=3.11 numpy pandas scipy scikit-learn \
    matplotlib seaborn networkx openpyxl -c conda-forge
conda activate cprostrata
pip install adjustText
```

## Code layout

Scripts resolve their inputs by relative path from `code/`. The pipeline has two
stages: per-analysis source panels, then composite panels montaged into the article's
final figures. `build_supp_tables.py` assembles the supplementary tables.

## Licence

Code: as stated in `CITATION.cff` / `.zenodo.json`.
