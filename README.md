# 2026_streit_opes

This repository contains Python code and Jupyter notebooks to reproduce the analyses and figures presented in:

> **Streit, J. O., Invernizzi, M., Bottaro, S., Tamiola, K., & Lindorff-Larsen, K. (2026).**
> *Transient tertiary structure in intrinsically disordered proteins revealed by multithermal enhanced sampling.*
> **bioRxiv** (URL TBA)

---

## Data availability

To run the analyses in this repository, additional data must be downloaded from the corresponding Zenodo repositories:

* **OPES multiT MD simulation data** and the **reweighted ACTR ensemble**:
  [https://zenodo.org/records/18249727](https://zenodo.org/records/18249727)

* **Unbiased simulations**, **REST2 simulations**, and **ELViM projection data**:
  [https://zenodo.org/records/18253563](https://zenodo.org/records/18253563)

Please download the data and place it in the root directory of this repository as expected by the analysis notebooks.

---

## Repository structure

The repository is organised as follows:

* **`analysis/`**
  OPES multiT and unbiased MD simulation analyses and comparisons.

* **`analysis_REST2/`**
  REST2 simulation analyses and comparisons with OPES.

* **`analysis_elvim/`**
  Analyses of conformational space exploration using ELViM projections.

* **`analysis_with_exp/`**
  Reweighting and analysis of OPES multiT ACTR simulations using experimental data.

* **`analysis_with_exp_rest2/`**
  Reweighting and analysis of REST2 ACTR simulations using experimental data.

* **`native_contacts_actr/`**
  Sets of "native" contact pairs and reference distances used to analyse ACTR simulations.

* **`pdbs_actr_state3/`**
  PDB files of the main conformational states of ACTR "state 3".

* **`plumed_files/`**
  PLUMED input files used for running and analysing OPES multiT simulations.

* **`scripts/`**
  Python scripts used to calculate experimental observables, radii of gyration, intramolecular contacts, and fractions of native contacts.

---

## Dependencies

A conda environment containing all required Python dependencies can be set up using `python_env.yml`.

Ensemble reweighting calculations require BME: [https://github.com/KULL-Centre/BME](https://github.com/KULL-Centre/BME)
