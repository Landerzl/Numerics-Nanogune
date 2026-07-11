"""
================================================================================
 Spatial IETS  dI/dV  mapping of the first magnetic excitation (spin gap Delta)
 of the Triangulene-Phenalenyl polymer   P-T-P-[P-T-P]_{N-1}   at  N = 3
 (12 spin-1/2 sites), bridged to a phenomenological Tersoff-Hamann STM image.
--------------------------------------------------------------------------------
 HOW TO USE
   * The IETS_Simulation folder contains:
         monomer3.xyz              <-- optimized geometry (input)
         iets_spatial_map.py       <-- THIS script
         dIdV_Maps/                <-- created automatically

   * Run it with a Python environment that has numpy, matplotlib and qutip.
     On execution it creates the sub-folder  dIdV_Maps/  and writes the map
     there.  The geometry is read from the local monomer3.xyz file.

 PHYSICS OUTLINE
   1. Build the 12-site Heisenberg Hamiltonian with the couplings
        J_PP = 38,  J_TP = 40,  J_T = -500  meV.
   2. Exactly diagonalize -> singlet ground state |0>  and the (3-fold
      degenerate) triplet first-excited multiplet {|1_m>}.
   3. Per-site spin-flip transition weight (basis-independent: summed over the
      degenerate final multiplet and over the three spin components)
          W_i = sum_{m in triplet} ( |<1_m|Sx_i|0>|^2
                                     +|<1_m|Sy_i|0>|^2
                                     +|<1_m|Sz_i|0>|^2 ).
      For a single final state this reduces to the requested
          W_i = |<1|Sx_i|0>|^2 + |<1|Sy_i|0>|^2 + |<1|Sz_i|0>|^2 ,
      but the multiplet sum removes the arbitrariness of the numerical
      triplet basis returned by the solver.
   4. Phenomenological constant-height Tersoff-Hamann map
          I(x,y) ~ sum_i W_i sum_{alpha in site i}
                     exp(-kappa sqrt((x-x_a)^2 + (y-y_a)^2 + z0^2)) .
================================================================================
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

# ================================================================
#  0.  Model / imaging parameters
# ================================================================
N_MON = 3                     # number of monomers -> 4*N_MON = 12 spin sites
J_PP  =  38.0                 # phenalenyl - phenalenyl    (AFM, meV)
J_TP  =  40.0                 # triangulene - phenalenyl   (AFM, meV)
J_T   = -500.0                # intra-triangulene Ta - Tb  (FM,  meV)

KAPPA = 1.0                   # inverse vacuum decay length  (1/Angstrom)
Z0    = 4.5                   # constant tip height          (Angstrom)
GRID_STEP = 0.15              # real-space sampling of the (x,y) map (Angstrom)
MARGIN    = 4.0               # padding around the molecule in the map (Angstrom)

# ================================================================
#  1.  Heisenberg Hamiltonian for the P-T-P polymer (4N spin-1/2 sites)
#      Per monomer the 4 sites are ordered  [P, Ta, Tb, P].
# ================================================================
def build_bonds(N):
    """(i, j, J) list for the 4N-1 open-boundary bonds."""
    bonds = []
    for m in range(N):
        b = 4 * m
        bonds.append((b,     b + 1, J_TP))   # P  - Ta
        bonds.append((b + 1, b + 2, J_T ))   # Ta - Tb  (ferromagnetic)
        bonds.append((b + 2, b + 3, J_TP))   # Tb - P
        if m < N - 1:
            bonds.append((b + 3, b + 4, J_PP))   # P - P (inter-monomer)
    return bonds


def build_spin_operators(n):
    """Single-site S=sigma/2 operators embedded in the 2**n Hilbert space."""
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
    H = 0
    for i, j, J in bonds:
        H += J * (Sx[i] * Sx[j] + Sy[i] * Sy[j] + Sz[i] * Sz[j])
    return H


def total_spin(state, S2):
    s2 = max(qt.expect(S2, state), 0.0)
    return 0.5 * (np.sqrt(1.0 + 4.0 * s2) - 1.0)


# ================================================================
#  2.  Transition weights  W_i  (singlet -> triplet spin-flip)
# ================================================================
def transition_weights(n_sites):
    """Return (W, info) with W[i] the per-site spin-flip transition weight."""
    Sx, Sy, Sz = build_spin_operators(n_sites)
    H  = build_hamiltonian(build_bonds(N_MON), Sx, Sy, Sz)
    S2 = sum(Sx) ** 2 + sum(Sy) ** 2 + sum(Sz) ** 2

    # dense ED (2**12 = 4096): trivial and gives clean degenerate eigenvectors
    evals, evecs = H.eigenstates()

    E0, gs = evals[0], evecs[0]
    S_gs   = total_spin(gs, S2)

    # first-excited multiplet: all states degenerate with evals[1]
    E1  = evals[1]
    tol = 1e-6 * max(1.0, abs(E1))
    multiplet = [evecs[k] for k in range(len(evals)) if abs(evals[k] - E1) < tol]
    S_ex = total_spin(multiplet[0], S2)

    # W_i = sum over degenerate final states and spin components of |<f|S^i|0>|^2
    W = np.zeros(n_sites)
    for i in range(n_sites):
        for f in multiplet:
            for Sop in (Sx[i], Sy[i], Sz[i]):
                W[i] += abs(Sop.matrix_element(f, gs)) ** 2

    info = dict(E0=float(np.real(E0)), E1=float(np.real(E1)),
                gap=float(np.real(E1 - E0)), S_gs=S_gs, S_ex=S_ex,
                mult=len(multiplet))
    return W, info


# ================================================================
#  3.  Geometry:  monomer3.xyz  ->  12 spin sites
# ================================================================
#  Fragment sequence along +x for N=3 :  P - T - P - P - T - P - P - T - P
#  (6 phenalenyls of 13 C, 3 triangulenes of 22 C; 6*13 + 3*22 = 144 C).
#  The atoms in monomer3.xyz are listed fragment-by-fragment along the chain,
#  so a cumulative-count partition reproduces the physical fragments (verified
#  against the coordinates).  Each phenalenyl -> one spin site; each triangulene
#  -> two spin sites (Ta = lower-x half, Tb = higher-x half), matching the
#  [P, Ta, Tb, P] ordering of build_bonds.
#
#      frag :  P    T      P    P    T      P    P    T      P
#      site :  0   1,2     3    4   5,6     7    8   9,10   11
# ----------------------------------------------------------------
FRAG_KIND = ['P', 'T', 'P', 'P', 'T', 'P', 'P', 'T', 'P']
FRAG_SIZE = [13,  22,  13,  13,  22,  13,  13,  22,  13]     # carbons per fragment


def load_xyz_to_sites(filename=None):
    """Parse the Carbon atoms of monomer3.xyz and group them into the 12 spin
    sites 0..11.  Returns  site_atoms : dict{site_index -> (M,2) array of (x,y)}."""
    if filename is None:
        here = os.path.dirname(os.path.abspath(__file__))
        filename = os.path.join(here, 'monomer3.xyz')

    # --- read, skip the 2 header lines, keep only Carbons ---
    with open(filename) as fh:
        lines = fh.readlines()[2:]                 # skip count line + comment line
    carbons = []
    for ln in lines:
        tok = ln.split()
        if len(tok) >= 4 and tok[0] == 'C':
            carbons.append((float(tok[1]), float(tok[2])))
    carbons = np.array(carbons)

    if len(carbons) != sum(FRAG_SIZE):
        raise ValueError(f"Expected {sum(FRAG_SIZE)} carbons for the N=3 chain, "
                         f"found {len(carbons)} in {filename}.")

    # --- cumulative-count partition into the 9 fragments, then -> 12 sites ---
    site_atoms = {}
    cursor, site = 0, 0
    for kind, size in zip(FRAG_KIND, FRAG_SIZE):
        frag = carbons[cursor:cursor + size]
        cursor += size
        if kind == 'P':
            site_atoms[site] = frag
            site += 1
        else:                                      # triangulene -> Ta, Tb
            order = np.argsort(frag[:, 0])         # sort by x
            half  = len(order) // 2
            site_atoms[site]     = frag[order[:half]]    # Ta = lower-x half
            site_atoms[site + 1] = frag[order[half:]]    # Tb = higher-x half
            site += 2
    return site_atoms


# ================================================================
#  4.  Phenomenological constant-height Tersoff-Hamann dI/dV map
# ================================================================
def simulate_didv(site_atoms, W):
    """I(x,y) ~ sum_i W_i sum_{a in site i} exp(-kappa*sqrt(dx^2+dy^2+z0^2))."""
    all_x = np.concatenate([a[:, 0] for a in site_atoms.values()])
    all_y = np.concatenate([a[:, 1] for a in site_atoms.values()])

    xs = np.arange(all_x.min() - MARGIN, all_x.max() + MARGIN, GRID_STEP)
    ys = np.arange(all_y.min() - MARGIN, all_y.max() + MARGIN, GRID_STEP)
    X, Y = np.meshgrid(xs, ys)
    I = np.zeros_like(X)

    for i, atoms in site_atoms.items():
        if W[i] == 0.0:
            continue
        for (xa, ya) in atoms:
            I += W[i] * np.exp(-KAPPA * np.sqrt((X - xa) ** 2 + (Y - ya) ** 2 + Z0 ** 2))

    return X, Y, I, (all_x, all_y)


# ================================================================
#  5.  Driver
# ================================================================
def main():
    here     = os.path.dirname(os.path.abspath(__file__))
    out_dir  = os.path.join(here, 'dIdV_Maps')
    os.makedirs(out_dir, exist_ok=True)           # create output folder

    n_sites = 4 * N_MON
    W, info = transition_weights(n_sites)

    print(f"Ground state : E0 = {info['E0']:.4f} meV   (S = {info['S_gs']:.2f})")
    print(f"Excited mult.: E1 = {info['E1']:.4f} meV   (S = {info['S_ex']:.2f}, "
          f"degeneracy {info['mult']})")
    print(f"Spin gap     : Delta = {info['gap']:.4f} meV   (eV_bias = Delta)")
    labels = ['P', 'Ta', 'Tb', 'P'] * N_MON
    print("\nPer-site spin-flip transition weights W_i:")
    for i in range(n_sites):
        print(f"   site {i:2d} ({labels[i]:>2s}) :  W = {W[i]:.4f}")

    site_atoms = load_xyz_to_sites()
    X, Y, I, (ax_, ay_) = simulate_didv(site_atoms, W)
    I /= I.max()                                   # normalize for display

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(14, 4.2))
    pcm = ax.pcolormesh(X, Y, I, cmap='magma', shading='auto')
    ax.scatter(ax_, ay_, s=6, c='cyan', alpha=0.35,
               edgecolors='none', label='C atoms')

    # per-site index annotations at site centroids
    for i, atoms in site_atoms.items():
        cx, cy = atoms[:, 0].mean(), atoms[:, 1].mean()
        ax.text(cx, cy, str(i), color='white', fontsize=7,
                ha='center', va='center', fontweight='bold')

    ax.set_aspect('equal')
    ax.set_xlabel(r'$x\ (\mathrm{\AA})$')
    ax.set_ylabel(r'$y\ (\mathrm{\AA})$')
    ax.set_title(r'Simulated IETS $dI/dV$ map of the first magnetic excitation '
                 rf'($\Delta = {info["gap"]:.2f}$ meV) — P-T-P polymer, $N=3$',
                 fontsize=11)
    cb = fig.colorbar(pcm, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label(r'$dI/dV$  (norm.)')
    ax.legend(loc='upper right', fontsize=8, framealpha=0.3)

    fig.tight_layout()
    out_png = os.path.join(out_dir, 'iets_didv_map_N3.png')
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    print(f"\nSaved dI/dV map to {out_png}")


if __name__ == '__main__':
    main()
