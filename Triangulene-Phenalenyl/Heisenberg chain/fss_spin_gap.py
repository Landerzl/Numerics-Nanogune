"""
Finite-Size Scaling (FSS) of the singlet-triplet spin gap of the
Triangulene-Phenalenyl polymer   P-T-P-[P-T-P]_{N-1}.

This is a standalone sweep script. It does NOT redefine the Hamiltonian
construction: build_bonds / build_spin_operators / build_hamiltonian are
imported from the sibling script spin_chain_qutip.py, which also carries the
couplings  J_PP = 38, J_TP = 40, J_T = -500 meV  (build_bonds injects them).

For each N = 1..5 we take only the two lowest eigenvalues (E0, E1) with a
sparse solver and form the gap  Delta = E1 - E0, then extrapolate Delta vs 1/N
to the thermodynamic limit.
"""

import os
import importlib.util

import numpy as np
import matplotlib
matplotlib.use("Agg")            # headless: keep the imported module's plt.show() a no-op
import matplotlib.pyplot as plt

# ------------------------------------------------------------
#  Reuse the existing builders from spin_chain_qutip.py
# ------------------------------------------------------------
_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "spin_chain_qutip", os.path.join(_here, "spin_chain_qutip.py"))
_scq = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_scq)    # runs the sibling script once (its N=3 demo)

build_bonds           = _scq.build_bonds
build_spin_operators  = _scq.build_spin_operators
build_hamiltonian     = _scq.build_hamiltonian

print("\n" + "=" * 60)
print("Finite-size scaling sweep  (Triangulene-Phenalenyl)")
print("=" * 60)

# ============================================================
#  Sweep  N = 1 .. 5
# ============================================================
N_list = np.arange(1, 6)                 # N = 1, 2, 3, 4, 5
gaps   = np.empty(len(N_list))

for k, N in enumerate(N_list):
    n             = 4 * N
    bonds_N       = build_bonds(N)                       # same J couplings
    SxN, SyN, SzN = build_spin_operators(n)
    HN            = build_hamiltonian(bonds_N, SxN, SyN, SzN)

    # Sparse solver: only the two lowest eigenvalues (E0, E1).
    # QuTiP's sparse path (ARPACK, which='SA'); equivalently one could do
    #   from scipy.sparse.linalg import eigsh
    #   E0, E1 = np.sort(eigsh(HN.data.as_scipy(), k=2, which='SA',
    #                          return_eigenvectors=False))
    E0, E1 = HN.eigenenergies(sparse=True, eigvals=2, sort='low')

    gaps[k] = E1 - E0
    print(f"N = {N}  ({n:2d} sites, dim 2^{n} = {2**n:>7d})   "
          f"E0 = {E0:11.4f}   E1 = {E1:11.4f}   Delta = {gaps[k]:8.4f} meV")

# ---- extrapolate the gap to the thermodynamic limit  1/N -> 0 ----
inv_N            = 1.0 / N_list
slope, intercept = np.polyfit(inv_N, gaps, 1)     # linear:  Delta = a*(1/N) + b
xfit             = np.linspace(0.0, inv_N.max() * 1.05, 100)
yfit             = slope * xfit + intercept
print(f"\nLinear extrapolation  Delta(1/N -> 0) = {intercept:.4f} meV "
      f"(slope = {slope:.4f})")

# ============================================================
#  Plot :  Delta vs 1/N (with extrapolation)   |   Delta vs N (raw trend)
# ============================================================
fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.5))

# left panel : Delta vs 1/N  +  linear extrapolation to 1/N = 0
axL.plot(inv_N, gaps, 'o', ms=8, color='#cc0000', label='ED data')
axL.plot(xfit, yfit, '--', color='#1f5fbf',
         label=fr'linear fit  $\to$  {intercept:.3f} meV')
axL.plot(0.0, intercept, 's', ms=10, mfc='none', mec='#1f5fbf', mew=2,
         label=r'extrapolation  $1/N \to 0$')
axL.set_xlabel(r'$1/N$')
axL.set_ylabel(r'spin gap  $\Delta = E_1 - E_0$  (meV)')
axL.set_title('Finite-size scaling of the spin gap', fontweight='bold')
axL.set_xlim(left=-0.02)
axL.legend(fontsize=9)
axL.grid(alpha=0.3)

# right panel : raw trend  Delta vs N
axR.plot(N_list, gaps, 'o-', ms=8, color='#cc0000')
axR.set_xlabel(r'$N$  (monomers)')
axR.set_ylabel(r'$\Delta$  (meV)')
axR.set_title('Raw trend', fontweight='bold')
axR.set_xticks(N_list)
axR.grid(alpha=0.3)

plt.tight_layout()
_outfile = os.path.join(_here, 'spin_gap_fss.png')
plt.savefig(_outfile, dpi=200, bbox_inches='tight')
print(f"Saved figure to {_outfile}")
