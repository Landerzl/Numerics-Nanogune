# Magnetic Coupling in 7-AGNRs

This repository contains scripts and calculations to study the magnetic exchange coupling ($J$) in finite armchair graphene nanoribbons of the 7-AGNR type.

## Contents

- **`Main/sisl_hubbard_7agnr.py`**: This is the main script of the project. It builds the geometry of the nanoribbons for different lengths, defines the Tight-Binding Hamiltonian, and calculates $J$ using two different approaches for comparison:
  1. **Effective Hubbard Dimer**: A simplified model that extracts $t_{\mathrm{eff}}$ and $U_{\mathrm{eff}}$ to yield $J = 4t_{\mathrm{eff}}^2 / U_{\mathrm{eff}}$.
  2. **Mean-Field Hubbard (MFH)**: An iterative self-consistent calculation that evaluates the actual energy difference between the Ferromagnetic (FM) and Antiferromagnetic (AFM) states.
- **`Main/23-06-2026.py`**: An auxiliary script using only `numpy` and `scipy` to verify the dimer rule and visualize the probability density of the topological states.
- **`Main/plot_spin_polarization_en.py`**: A utility to visualize how spins are distributed along the zigzag edges of the nanoribbon.
- **`Main/test_geom.py`**: A quick test script to ensure that pruning the ribbon doesn't leave unwanted dangling bonds.

## Results

Running the main script generates two comparative plots directly in your directory:
1. **`J_comparison_vs_L.png`**: Compares how the magnetic coupling $J$ decays exponentially as the ribbon length ($L$) increases. It includes useful reference lines such as the superconducting gap $2\Delta$ of Nb(110) at 3.0 meV and the thermal energy scale of an STM operating at 1.2 K. This allows you to easily see whether $J$ falls within the gap.
2. **`J_diff_vs_L.png`**: A plot displaying the absolute energy difference (the deviation) between the simple Dimer model and the full MFH simulation.

## Requirements

To run the scripts successfully, your Python virtual environment needs the following installed:
- `numpy` and `scipy`
- `matplotlib` (for generating the figures)
- `sisl` (the key library used to easily define complex geometries and Hamiltonians)

All the steps are thoroughly documented within the code through comments.
