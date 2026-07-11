# Numerical Studies: Magnetic Coupling in Nanographene Systems

This repository contains computational tools and scripts to study magnetic exchange coupling ($J$) in finite nanographene systems. Two main projects are included:

1. **7-AGNR** — Magnetic coupling in finite 7-Armchair Graphene Nanoribbons.
2. **Triangulene-Phenalenyl** — Heisenberg spin-chain model and IETS simulation of phenalenyl–[3]triangulene oligomers.

---

## Repository Structure

```
.
├── 7-AGNR/
│   └── Main/
│       ├── sisl_hubbard_7agnr.py
│       ├── predict_J.py
│       ├── plot_spin_polarization_en.py
│       ├── test_geom.py
│       └── Figures/
├── Triangulene-Phenalenyl/
│   ├── Geometry checks/
│   │   ├── build_v2.py
│   │   ├── build_monomer.py
│   │   └── units.py
│   ├── Heisenberg chain/
│   │   ├── spin_chain_qutip.py
│   │   ├── verify_ed.py
│   │   ├── fss_spin_gap.py
│   │   ├── ed_result.npz
│   │   ├── model.tex / model.pdf
│   │   ├── spin_chain_N3.png
│   │   ├── spin_chain_ed.png
│   │   └── spin_gap_fss.png
│   ├── IETS_Simulation/
│   │   ├── iets_spatial_map.py
│   │   ├── iets_waterfall.py
│   │   ├── monomer3.xyz
│   │   ├── iets_theory.tex / iets_theory.pdf
│   │   ├── dIdV_Maps/           ← created on run
│   │   └── Waterfall_Plots/     ← created on run
│   └── Spin Polarization check/
│       ├── 1 monomer/
│       │   ├── SpinPlotTriangulenes.py
│       │   └── monomer1.xyz
│       └── 3 monomer/
│           ├── SpinPlotTriangulenes.py
│           └── monomer3.xyz
├── Presentations/
│   ├── presentation.tex / presentation.pdf
│   └── presentation_short.tex / presentation_short.pdf
└── figs/
    ├── J_MFH_vs_L.png
    ├── J_comparison_vs_L.png
    ├── J_dimer_vs_L.png
    ├── J_diff_vs_L.png
    ├── spin_polarization_L6.png
    ├── spin_polarization_L12.png
    ├── oligomer_closed.png
    ├── spin_density_monomer1.png
    ├── spin_density_monomer3.png
    ├── spin_chain_N3.png
    ├── spin_chain_N3_crop.png
    ├── spin_gap_fss.png
    ├── iets_didv_map_N3.png
    ├── iets_waterfall_N3.png
    └── test_geom_L4.png
```

---

## 1. 7-AGNR: Magnetic Coupling Analysis

### Overview

This project systematically compares two theoretical approaches to evaluate the magnetic exchange coupling $J$ as a function of nanoribbon length $L$ (number of DBBA monomers).

### Scripts

- **`sisl_hubbard_7agnr.py`** — Core simulation script. Uses `sisl` to construct finite 7-AGNR geometries and evaluates $J$ via two methods:
  1. **Effective Hubbard Dimer**: Extracts the effective hopping ($t_{\mathrm{eff}}$) and Coulomb repulsion ($U_{\mathrm{eff}}$) to estimate $J = 4t_{\mathrm{eff}}^2 / U_{\mathrm{eff}}$.
  2. **Mean-Field Hubbard (MFH)**: Performs a self-consistent calculation to find the total energy difference between the Ferromagnetic (FM) and Antiferromagnetic (AFM) ground states ($J = E_{\mathrm{FM}} - E_{\mathrm{AFM}}$).

- **`predict_J.py`** — CLI tool for predicting the number of DBBA monomers and extrapolated $J_{\mathrm{MFH}}$ from an experimental nanoribbon length (in nm). Fits an exponential decay to computed MFH data and extrapolates to the given length.

  **Usage:**
  ```bash
  python predict_J.py <length_nm>
  ```
  *Example:* `python predict_J.py 4.5` will determine the closest integer $L$, compute MFH data up to that $L$, fit an exponential, and report both the continuous equivalent monomer count and the predicted $J$ in meV.

- **`plot_spin_polarization_en.py`** — Utility to visualize spin polarization at the zigzag edges.

- **`test_geom.py`** — Verifies that generated nanoribbon geometries are properly closed with no dangling bonds.

### Key Results

The calculations yield plots showing the exponential decay of $J$ with the number of DBBA monomers ($L$). Reference scales included: the superconducting gap of Nb(110) ($2\Delta \approx 3.0$ meV) and the thermal energy at standard STM operating temperatures ($k_B T$ at $1.2$ K).

---

## 2. Triangulene-Phenalenyl

### Overview

This project models oligomers of phenalenyl–[3]triangulene (P–T–P) nanographene units as a 1D Heisenberg spin chain, extracts the spin gap and local spin texture via exact diagonalization (ED), and simulates the corresponding IETS/STM signatures.

**Polymer sequence:** $[\mathrm{P}$–$\mathrm{T}$–$\mathrm{P}]_N$ — each monomer contributes 4 spin-$\frac{1}{2}$ sites:
- $P$: phenalenyl radical (one $S = \frac{1}{2}$ site)
- $T_a, T_b$: the [3]triangulene $S=1$ unit (two $S = \frac{1}{2}$ sites locked by a strong FM exchange $J_T$)

**Hamiltonian convention:** $H = \sum_{\langle ij \rangle} J_{ij}\, \mathbf{S}_i \cdot \mathbf{S}_j$

**Exchange parameters (meV):**

| Bond | Type | $J$ (meV) |
|---|---|---|
| Phenalenyl – Phenalenyl ($J_{PP}$) | AFM | +38.0 |
| Triangulene – Phenalenyl ($J_{TP}$) | AFM | +40.0 |
| Intra-triangulene $T_a$–$T_b$ ($J_T$) | FM  | −500.0 |

### Subfolders

#### `Geometry checks/`

Tools to construct and visualize the phenalenyl–triangulene oligomer geometry.

- **`units.py`** — Defines the atomic positions for the phenalenyl (P) and [3]triangulene (T) monomers, including radical corner identification.
- **`build_v2.py`** — Assembles an alternating P–T–P chain with biphenyl-type C–C links and alternating apex orientation (horizontal zigzag layout).
- **`build_monomer.py`** — Builds a fully H-passivated oligomer (no radicals), writes it to `.xyz`, and saves a bond visualization (`oligomer_closed.png`).

#### `Heisenberg chain/`

Exact diagonalization of the Heisenberg model for the P–T–P chain.

- **`spin_chain_qutip.py`** — Full ED using `QuTiP`. Builds the tensor-product Hilbert space, diagonalizes the Hamiltonian, and plots:
  - The ground-state local magnetization $\langle S_z^i \rangle$ (identically zero for a singlet by SU(2) symmetry).
  - The local spin texture of the $S_z = +1$ member of the lowest triplet (reveals the spatial distribution of magnetic moment).
  
  Saves figures as `spin_chain_N{N}.png`.

- **`verify_ed.py`** — Lightweight cross-check of the ED using `scipy.sparse`. Computes the same spectrum and spin texture without QuTiP. Saves results to `ed_result.npz` for further analysis.

- **`fss_spin_gap.py`** — Finite-size scaling (FSS) of the singlet–triplet spin gap. Sweeps $N = 1\ldots5$ monomers using a sparse ARPACK solver (via QuTiP), extracts $\Delta(N) = E_1 - E_0$, and extrapolates linearly in $1/N$ to estimate the thermodynamic limit gap. Saves a two-panel figure (`spin_gap_fss.png`) showing both the FSS plot and the raw trend.

- **`ed_result.npz`** — Cached ED output (energies, spin gap, local magnetizations) for $N = 3$ monomers.

- **`model.tex / model.pdf`** — LaTeX document summarizing the model, Hamiltonian, and numerical results.

#### `IETS_Simulation/`

Simulation of the inelastic electron tunneling spectroscopy (IETS) signal expected from the P–T–P polymer in an STM experiment, using the Tersoff–Hamann approximation.

- **`iets_spatial_map.py`** — Computes the spatially resolved $dI/dV$ map of the first magnetic excitation (singlet → triplet spin flip). For each grid point $(x,y)$ above the molecule it evaluates the Tersoff–Hamann tunneling current weighted by the per-site spin-flip transition weights $W_i$ obtained from ED. Saves `dIdV_Maps/iets_didv_map_N3.png`.

- **`iets_waterfall.py`** — Simulates point-spectroscopy $dI/dV(V)$ spectra for each of the 12 spin sites (tip parked at the site centroid). The inelastic step amplitude at each site is set by the full molecule-wide Tersoff–Hamann decay sum. Spectra are stacked as a waterfall plot. Saves `Waterfall_Plots/iets_waterfall_N3.png`.

- **`monomer3.xyz`** — Optimized geometry of the $N=3$ P–T–P oligomer (input for both IETS scripts).

- **`iets_theory.tex / iets_theory.pdf`** — LaTeX document describing the IETS model and results.

> **Workflow:** Both IETS scripts read `monomer3.xyz` from the same folder and create their output subdirectories automatically on first run.

#### `Spin Polarization check/`

MFH spin density calculations on the actual molecular geometry (requires `.xyz` from Geometry checks).

- **`1 monomer/`** — Single P–T–P monomer: geometry (`monomer1.xyz`) + MFH spin density plot.
- **`3 monomer/`** — Three-monomer oligomer: geometry (`monomer3.xyz`) + MFH spin density plot.

> **Workflow:** Run `build_monomer.py` (Geometry checks) to generate the `.xyz`, then use `SpinPlotTriangulenes.py` in the corresponding monomer folder to compute and visualize the Hubbard spin density.

---

## Dependencies

### 7-AGNR

```
numpy
scipy
matplotlib
sisl       # geometry construction and tight-binding Hamiltonians
hubbard    # local module for Mean-Field Hubbard calculations
```

### Triangulene-Phenalenyl

```
numpy
scipy
matplotlib
sisl       # geometry construction (Geometry checks, Spin Polarization check)
qutip      # exact diagonalization (Heisenberg chain, IETS Simulation)
hubbard    # Mean-Field Hubbard (Spin Polarization check)
```

---

## Presentations

All presentation sources and compiled PDFs are located in `Presentations/`:

- [`presentation.pdf`](./Presentations/presentation.pdf) — Full presentation of results.
- [`presentation_short.pdf`](./Presentations/presentation_short.pdf) — Condensed version.
