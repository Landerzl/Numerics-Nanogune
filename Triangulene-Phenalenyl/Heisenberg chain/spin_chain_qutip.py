"""
Exact diagonalization of the 1D Heisenberg Hamiltonian of a
phenalenyl-triangulene nanographene polymer, using QuTiP.

Polymer sequence  P-T-P-[P-T-P]_{N-1}  ->  4N spin-1/2 sites.
Per monomer the 4 sites are ordered  [P, Ta, Tb, P]:
    * P        : phenalenyl, one S=1/2 site
    * Ta, Tb   : the [3]triangulene S=1 unit, two S=1/2 sites locked
                 by a strong ferromagnetic exchange J_T < 0.

Convention:  H = sum_{<ij>} J_ij  S_i . S_j     with  S = sigma/2 .

Bond rules (4N-1 bonds, open boundary conditions):
    inside a monomer :  J_TP , J_T , J_TP
    between monomers :  J_PP
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import qutip as qt

# ============================================================
#  Parameters
# ============================================================
N_monomers = 3            # system size  ->  n_sites = 4 * N_monomers

J_PP =  38.0              # phenalenyl - phenalenyl        (AFM, meV)
J_TP =  40.0              # triangulene - phenalenyl       (AFM, meV)
J_T  = -500.0             # intra-triangulene Ta - Tb      (FM,  meV)


# ============================================================
#  Lattice / bond generation
# ============================================================
def build_bonds(N):
    """Return the list of (i, j, J) tuples for the 4N-1 bonds."""
    bonds = []
    for m in range(N):
        b = 4 * m                       # first site of monomer m  (a P)
        bonds.append((b,     b + 1, J_TP))   # P  - Ta
        bonds.append((b + 1, b + 2, J_T ))   # Ta - Tb   (ferromagnetic)
        bonds.append((b + 2, b + 3, J_TP))   # Tb - P
        if m < N - 1:
            bonds.append((b + 3, b + 4, J_PP))   # P - P  (inter-monomer)
    return bonds


def site_labels(N):
    labels, kinds = [], []
    for _ in range(N):
        labels += ['P', 'Ta', 'Tb', 'P']
        kinds  += ['P', 'T',  'T',  'P']
    return labels, kinds


# ============================================================
#  Spin operators on the full 2**n_sites Hilbert space
# ============================================================
def build_spin_operators(n):
    """Return lists Sx, Sy, Sz of single-site spin-1/2 operators embedded
    in the full tensor-product space (S = sigma/2)."""
    sx, sy, sz = 0.5 * qt.sigmax(), 0.5 * qt.sigmay(), 0.5 * qt.sigmaz()

    def embed(op, i):
        ops = [qt.qeye(2)] * n
        ops[i] = op
        return qt.tensor(ops)

    Sx = [embed(sx, i) for i in range(n)]
    Sy = [embed(sy, i) for i in range(n)]
    Sz = [embed(sz, i) for i in range(n)]
    return Sx, Sy, Sz


def build_hamiltonian(bonds, Sx, Sy, Sz):
    """H = sum J_ij ( Sx_i Sx_j + Sy_i Sy_j + Sz_i Sz_j )."""
    H = 0
    for i, j, J in bonds:
        H += J * (Sx[i] * Sx[j] + Sy[i] * Sy[j] + Sz[i] * Sz[j])
    return H


# ============================================================
#  Helpers
# ============================================================
def total_spin(state, S2):
    """Return S from  <S^2> = S(S+1)."""
    s2 = max(qt.expect(S2, state), 0.0)
    return 0.5 * (np.sqrt(1.0 + 4.0 * s2) - 1.0)


def polarized_member(states, energies, Sz_tot, E_target, target_sz=1.0, tol=1e-3):
    """Given a (near-)degenerate multiplet at energy E_target, return the linear
    combination that is an eigenstate of total Sz with eigenvalue ~ target_sz.
    Used to expose the local spin texture of a magnetic multiplet."""
    sub = [psi for psi, E in zip(states, energies) if abs(E - E_target) < tol]
    d = len(sub)
    M = np.array([[Sz_tot.matrix_element(sub[a], sub[b]) for b in range(d)]
                  for a in range(d)])
    w, U = np.linalg.eigh(M)
    coeff = U[:, np.argmin(np.abs(w - target_sz))]
    psi = sum(coeff[a] * sub[a] for a in range(d))
    return psi.unit()


# ============================================================
#  Build & diagonalize
# ============================================================
n_sites = 4 * N_monomers
bonds   = build_bonds(N_monomers)
labels, kinds = site_labels(N_monomers)

Sx, Sy, Sz = build_spin_operators(n_sites)
H = build_hamiltonian(bonds, Sx, Sy, Sz)

Sz_tot = sum(Sz)
S2_tot = sum(Sx) ** 2 + sum(Sy) ** 2 + sum(Sz) ** 2

n_states = min(12, 2 ** n_sites - 1)
evals, evecs = H.eigenstates(sparse=True, eigvals=n_states)

# ---- ground state ----
E0   = evals[0]
gs   = evecs[0]
Sz0  = qt.expect(Sz_tot, gs)
S_gs = total_spin(gs, S2_tot)

# ---- first excited state and spin gap ----
E1   = evals[1]
gap  = E1 - E0
S_e1 = total_spin(evecs[1], S2_tot)

print(f"N_monomers = {N_monomers}   ->   {n_sites} sites,   Hilbert dim = 2^{n_sites} = {2**n_sites}")
print(f"Ground state : E0 = {E0:.4f} meV   Sz = {Sz0:+.3f}   S = {S_gs:.3f}")
print(f"1st excited  : E1 = {E1:.4f} meV                    S = {S_e1:.3f}")
print(f"Spin gap     : Delta = {gap:.4f} meV")

# ============================================================
#  Local magnetization
# ============================================================
# (1) ground state : for a total singlet this is identically zero (SU(2)).
mz_gs = np.array([qt.expect(Sz[i], gs) for i in range(n_sites)])

# (2) spin texture : Sz=+1 member of the lowest S>0 multiplet.
S_all = [total_spin(evecs[a], S2_tot) for a in range(n_states)]
mag_idx = next((a for a in range(1, n_states) if S_all[a] > 0.9), 1)
psi_pol = polarized_member(evecs, evals, Sz_tot, evals[mag_idx], target_sz=1.0)
mz_pol  = np.array([qt.expect(Sz[i], psi_pol) for i in range(n_sites)])

# ============================================================
#  Plot
# ============================================================
x = np.arange(n_sites)
colors = ['#cc0000' if k == 'T' else '#1f5fbf' for k in kinds]

fig, axes = plt.subplots(2, 1, figsize=(1.0 * n_sites + 2, 7), sharex=True)

axes[0].bar(x, mz_gs, color=colors, edgecolor='k', linewidth=0.4)
axes[0].axhline(0, color='gray', lw=0.6)
axes[0].set_ylabel(r'$\langle S_z^i\rangle$')
axes[0].set_ylim(-0.6, 0.6)
axes[0].set_title(f'Ground state  (S = {S_gs:.2f} singlet):  '
                  r'$\langle S_z^i\rangle = 0$ by SU(2) symmetry'
                  f'    —    N = {N_monomers},  E0 = {E0:.2f} meV',
                  fontweight='bold')
axes[0].legend(handles=[Patch(facecolor='#cc0000', label='Triangulene (Ta, Tb)'),
                        Patch(facecolor='#1f5fbf', label='Phenalenyl (P)')],
               fontsize=9, loc='upper right')

axes[1].bar(x, mz_pol, color=colors, edgecolor='k', linewidth=0.4)
axes[1].axhline(0, color='gray', lw=0.6)
axes[1].set_ylabel(r'$\langle S_z^i\rangle$')
axes[1].set_xlabel('site index')
axes[1].set_title(r'Local spin texture: $S_z=+1$ member of the lowest triplet '
                  f'($\\Delta$ = {gap:.2f} meV above the singlet)',
                  fontweight='bold')
axes[1].set_xticks(x)
axes[1].set_xticklabels([f'{i}\n{labels[i]}' for i in range(n_sites)], fontsize=8)

plt.tight_layout()
plt.savefig(f'spin_chain_N{N_monomers}.png', dpi=200, bbox_inches='tight')
print(f"\nSaved figure to spin_chain_N{N_monomers}.png")
plt.show()
