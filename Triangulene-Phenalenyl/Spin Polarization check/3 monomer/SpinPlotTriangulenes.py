"""
Ground-state spin density of ONE monomer (phenalenyl-triangulene-phenalenyl,
single-bond linked, fully H-passivated) via mean-field Hubbard.

Adapted from the 7-AGNR AFM/FM example. Key physics for THIS molecule:
the pi network is bipartite with N_A = N_B = 24  ->  by Lieb's theorem the
ground state is an OPEN-SHELL SINGLET (S = 0). The triangulene carries a
local S = 1 moment that couples antiferromagnetically to the two S = 1/2
phenalenyl moments, giving zero net spin but a rich local spin texture.

So the ground state = the AFM (open-shell singlet) branch, exactly analogous
to the AFM state in the AGNR example. We converge it from a staggered
(sublattice) seed and plot m_i = n_up_i - n_dn_i.

Requires: sisl, hubbard, numpy, matplotlib   (+ monomer1.xyz in the same folder)
"""
import sisl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from collections import deque
from hubbard import sp2, HubbardHamiltonian, density

# ──────────────────── Parameters ────────────────────
U_val = 3.0          # on-site repulsion (eV)
XYZ   = 'monomer3.xyz'
# ────────────────────────────────────────────────────

# ── Geometry: keep only the sp2 carbons (the pi system) ──
full = sisl.get_sile(XYZ).read_geometry()
carbon_idx = [ia for ia in full if full.atoms[ia].Z == 6]
geom = full.sub(carbon_idx)
geom.set_nsc([1, 1, 1])
na = geom.na
print(f"Monomer pi-system: {na} carbon atoms")

# ── Bonds + sublattice 2-coloring (for the staggered seed & plotting) ──
def get_bonds(g):
    bonds = []
    for ia in g:
        for j in g.close(ia, R=1.6):
            if j > ia:
                bonds.append((ia, j))
    return bonds

def sublattices(g):
    adj = [set(g.close(ia, R=1.6)) - {ia} for ia in g]
    color = -np.ones(g.na, dtype=int); color[0] = 0
    dq = deque([0])
    while dq:
        v = dq.popleft()
        for w in adj[v]:
            if color[w] < 0:
                color[w] = 1 - color[v]; dq.append(w)
    A = np.where(color == 0)[0]
    B = np.where(color == 1)[0]
    return A, B

bonds = get_bonds(geom)
A, B = sublattices(geom)
print(f"Sublattices: N_A={len(A)}, N_B={len(B)}  (imbalance {abs(len(A)-len(B))} -> Lieb S={abs(len(A)-len(B))/2})")

# ── 3NN TB Hamiltonian (spin polarized) ──
H_tb = sp2(geom, t1=-2.7, t2=-0.2, t3=-0.18, spin='polarized')

# ══════════════════════════════════════════════════════
# GROUND STATE = open-shell singlet (S = 0)
#   -> balanced filling (default q), staggered sublattice seed
# ══════════════════════════════════════════════════════
HH = HubbardHamiltonian(H_tb, U=U_val, kT=1e-6, nkpt=[1, 1, 1])
HH.set_polarization(A, dn=B)      # sublattice A up, B down  ->  AFM singlet
dn = HH.converge(density.calc_n, tol=1e-13, print_info=True)

m = HH.n[0] - HH.n[1]
E0 = HH.Etot
print(f"\nE_singlet = {E0:.10f} eV")
print(f"Total Sz  = {np.sum(m)/2:.6f}")
print(f"|m| max   = {np.max(np.abs(m)):.6f}")
print(f"sum|m|    = {np.sum(np.abs(m)):.6f}")

# ── Spin-density plot (same style as the AGNR example) ──
xyz = geom.xyz
vmax = max(np.max(np.abs(m)), 1e-6)

fig, ax = plt.subplots(figsize=(13, 5))
segments = [[xyz[i, :2], xyz[j, :2]] for (i, j) in bonds]
ax.add_collection(LineCollection(segments, colors='#888888', linewidths=1.2, zorder=1))

sizes = 120 + 700 * np.abs(m)     # emphasize magnetic sites
sc = ax.scatter(xyz[:, 0], xyz[:, 1], c=m, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                s=sizes, edgecolors='k', linewidths=0.5, zorder=2)

ax.set_aspect('equal')
ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
ax.set_title('Ground-state spin density (open-shell singlet, S=0)\n'
             f'phenalenyl-triangulene-phenalenyl,  MFH  U={U_val} eV',
             fontsize=12, fontweight='bold')
ax.text(0.02, 0.95, f'Total $S_z$ = {np.sum(m)/2:.3f}\n$\\Sigma|m_i|$ = {np.sum(np.abs(m)):.2f}',
        transform=ax.transAxes, fontsize=10, va='top',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))
cbar = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label('$m_i = n_\\uparrow - n_\\downarrow$', fontsize=12)

plt.tight_layout()
plt.savefig('spin_density_monomer_hubbard.png', dpi=200, bbox_inches='tight')
print("\nPlot saved to 'spin_density_monomer_hubbard.png'")
plt.show()