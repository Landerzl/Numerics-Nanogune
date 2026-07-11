"""
================================================================================
 Point-spectroscopy IETS  dI/dV(V)  waterfall for the Triangulene-Phenalenyl
 polymer  P-T-P-[P-T-P]_{N-1}  at  N = 3  (12 spin-1/2 sites).
--------------------------------------------------------------------------------
 For every one of the 12 spin sites we "park" the STM tip at the geometric
 centroid of the carbon atoms that build that site and simulate the local
 inelastic  dI/dV  spectrum.  The inelastic step amplitude at each parked
 position is fixed *physically* by the Tersoff-Hamann spatial sum evaluated at
 that centroid (no abstract per-site weight is used for the height directly:
 the height is the full molecule-wide decay sum seen from that tip position).

 PHYSICS OUTLINE
   1. Exact diagonalization of the 12-site Heisenberg chain
        J_PP = 38,  J_TP = 40,  J_T = -500  meV
      -> singlet ground state |0>  and 3-fold triplet {|1_m>}, gap Delta.
   2. Per-site spin-flip transition weights (basis-independent multiplet sum)
        W_j = sum_{m} sum_{mu=x,y,z} |<1_m|S_j^mu|0>|^2 .
   3. Tip parked at centroid (x_i, y_i) of site i.  The local inelastic
      amplitude collects the exponential-decay contribution of every carbon
      atom of every site, weighted by that site's W_j:
        A_i = sum_{j=0}^{11} W_j sum_{alpha in site j}
                 exp(-kappa sqrt((x_i-x_alpha)^2 + (y_i-y_alpha)^2 + z0^2)) .
   4. Spectrum at site i = flat elastic baseline + two symmetric inelastic
      steps at V = +/- Delta, each scaled by A_i and broadened by an error
      function of width Gamma (experimental/thermal resolution):
        dI/dV_i(V) = 1 + A_i/2 [ erf((V-Delta)/(sqrt(2) Gamma))
                                - erf((V+Delta)/(sqrt(2) Gamma)) ] .
      (V < -Delta and V > +Delta lie in the inelastic-open plateaus; the
       |V| < Delta window is the elastic gap.)
   5. All 12 spectra stacked with a constant vertical offset -> waterfall.

 OUTPUT
   IETS_Simulation/Waterfall_Plots/iets_waterfall_N3.png
================================================================================
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.special import erf
import qutip as qt

# ================================================================
#  0.  Model / imaging / spectroscopy parameters
# ================================================================
N_MON = 3                     # number of monomers -> 4*N_MON = 12 spin sites
J_PP  =  38.0                 # phenalenyl - phenalenyl    (AFM, meV)
J_TP  =  40.0                 # triangulene - phenalenyl   (AFM, meV)
J_T   = -500.0                # intra-triangulene Ta - Tb  (FM,  meV)

KAPPA = 1.0                   # inverse vacuum decay length      (1/Angstrom)
Z0    = 4.5                   # constant tip height              (Angstrom)

V_MIN, V_MAX = -30.0, 30.0    # bias sweep range                 (mV == meV)
N_V   = 1201                  # number of bias points
GAMMA = 1.0                   # inelastic-step broadening        (meV)

OFFSET = 1.3                  # vertical offset between stacked curves (a.u.)

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
#  2.  Transition weights  W_j  and spin gap Delta
# ================================================================
def transition_weights(n_sites):
    """Return (W, Delta, info) for the singlet -> triplet spin-flip."""
    Sx, Sy, Sz = build_spin_operators(n_sites)
    H  = build_hamiltonian(build_bonds(N_MON), Sx, Sy, Sz)
    S2 = sum(Sx) ** 2 + sum(Sy) ** 2 + sum(Sz) ** 2

    evals, evecs = H.eigenstates()          # dense ED (2**12 = 4096)

    E0, gs = evals[0], evecs[0]
    S_gs   = total_spin(gs, S2)

    E1  = evals[1]
    tol = 1e-6 * max(1.0, abs(E1))
    multiplet = [evecs[k] for k in range(len(evals)) if abs(evals[k] - E1) < tol]
    S_ex = total_spin(multiplet[0], S2)

    W = np.zeros(n_sites)
    for i in range(n_sites):
        for f in multiplet:
            for Sop in (Sx[i], Sy[i], Sz[i]):
                W[i] += abs(Sop.matrix_element(f, gs)) ** 2

    Delta = float(np.real(E1 - E0))
    info  = dict(E0=float(np.real(E0)), E1=float(np.real(E1)), Delta=Delta,
                 S_gs=S_gs, S_ex=S_ex, mult=len(multiplet))
    return W, Delta, info


# ================================================================
#  3.  Geometry:  monomer3.xyz  ->  12 spin sites
# ================================================================
#  Fragment sequence along +x for N=3 :  P - T - P - P - T - P - P - T - P
#  (6 phenalenyls of 13 C, 3 triangulenes of 22 C; 6*13 + 3*22 = 144 C).
#  The atoms in monomer3.xyz are listed fragment-by-fragment along the chain,
#  so a cumulative-count partition reproduces the physical fragments.  Each
#  phenalenyl -> one spin site; each triangulene -> two spin sites
#  (Ta = lower-x half, Tb = higher-x half), matching the [P, Ta, Tb, P] order.
#
#      frag :  P    T      P    P    T      P    P    T      P
#      site :  0   1,2     3    4   5,6     7    8   9,10   11
# ----------------------------------------------------------------
FRAG_KIND = ['P', 'T', 'P', 'P', 'T', 'P', 'P', 'T', 'P']
FRAG_SIZE = [13,  22,  13,  13,  22,  13,  13,  22,  13]     # carbons per fragment


def load_xyz_to_sites(filename=None):
    """Parse the Carbon atoms of monomer3.xyz and group them into the 12 spin
    sites 0..11 (left to right along x).  Returns
    site_atoms : dict{site_index -> (M,2) array of (x,y)}."""
    if filename is None:
        here = os.path.dirname(os.path.abspath(__file__))
        filename = os.path.join(here, 'monomer3.xyz')

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


def site_centroids(site_atoms):
    """Geometric centroid (mean x, mean y) of the carbons of each site."""
    n = len(site_atoms)
    C = np.zeros((n, 2))
    for i in range(n):
        C[i] = site_atoms[i].mean(axis=0)
    return C


# ================================================================
#  4.  Local Tersoff-Hamann amplitude at each parked centroid
# ================================================================
def local_amplitudes(site_atoms, centroids, W):
    """A_i = sum_j W_j sum_{alpha in site j} exp(-kappa*sqrt(dx^2+dy^2+z0^2)),
    evaluated with the tip parked at centroid (x_i, y_i) of site i."""
    n = len(site_atoms)
    # stack all carbons and a matching per-atom weight (its site's W_j)
    all_atoms = np.concatenate([site_atoms[j] for j in range(n)], axis=0)
    all_wj    = np.concatenate([np.full(len(site_atoms[j]), W[j])
                                for j in range(n)])
    A = np.zeros(n)
    for i in range(n):
        xi, yi = centroids[i]
        r = np.sqrt((xi - all_atoms[:, 0])**2 +
                    (yi - all_atoms[:, 1])**2 + Z0**2)
        A[i] = np.sum(all_wj * np.exp(-KAPPA * r))
    return A


# ================================================================
#  5.  Broadened inelastic dI/dV spectrum (two symmetric steps)
# ================================================================
def didv_spectrum(V, A, Delta, Gamma):
    """Flat elastic baseline (=1) + symmetric inelastic steps at V = +/-Delta,
    each of height A, broadened by an error function of width Gamma.

    Physics: the spin-flip inelastic channel is CLOSED inside the gap
    (|V| < Delta, elastic baseline) and OPENS for |V| > Delta, raising the
    conductance.  Hence conductance steps UP at V = +Delta (for V > +Delta)
    and at V = -Delta (for V < -Delta)."""
    step_pos = 0.5 * (1.0 + erf((V - Delta) / (np.sqrt(2.0) * Gamma)))   # opens for V > +Delta
    step_neg = 0.5 * (1.0 - erf((V + Delta) / (np.sqrt(2.0) * Gamma)))   # opens for V < -Delta
    # both terms are 0 inside |V| < Delta and rise to +1 on their open side
    return 1.0 + A * (step_pos + step_neg)


# ================================================================
#  6.  Driver
# ================================================================
def main():
    here    = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(here, 'Waterfall_Plots')
    os.makedirs(out_dir, exist_ok=True)            # create output folder

    n_sites = 4 * N_MON

    # --- spin physics ---
    W, Delta, info = transition_weights(n_sites)
    print(f"Ground state : E0 = {info['E0']:.4f} meV  (S = {info['S_gs']:.2f})")
    print(f"Excited mult.: E1 = {info['E1']:.4f} meV  (S = {info['S_ex']:.2f}, "
          f"degeneracy {info['mult']})")
    print(f"Spin gap     : Delta = {Delta:.4f} meV\n")

    # --- geometry & parked amplitudes ---
    site_atoms = load_xyz_to_sites()
    centroids  = site_centroids(site_atoms)
    A          = local_amplitudes(site_atoms, centroids, W)
    A_norm     = A / A.max()                        # normalize step heights

    labels = ['P', 'Ta', 'Tb', 'P'] * N_MON
    print("site  type   centroid (x, y)  [A]     A_norm")
    for i in range(n_sites):
        print(f"  {i:2d}  ({labels[i]:>2s})  "
              f"({centroids[i,0]:6.2f}, {centroids[i,1]:5.2f})  "
              f"{A[i]:.4f}  {A_norm[i]:.3f}")

    # --- spectra ---
    V = np.linspace(V_MIN, V_MAX, N_V)
    spectra = np.array([didv_spectrum(V, A_norm[i], Delta, GAMMA)
                        for i in range(n_sites)])

    # ================================================================
    #  Plot : stacked waterfall
    # ================================================================
    fig, ax = plt.subplots(figsize=(14.0, 8.0))   # landscape aspect, for wide-slide display
    colors = cm.viridis(np.linspace(0.0, 0.92, n_sites))

    for i in range(n_sites):
        y = spectra[i] + i * OFFSET
        ax.plot(V, y, color=colors[i], lw=1.8)
        ax.fill_between(V, i * OFFSET + 1.0, y, color=colors[i], alpha=0.12)
        # right-hand label with site index and building-block type
        ax.text(V_MAX + 2.0, i * OFFSET + 1.0,
                f"site {i} ({labels[i]})",
                va='center', ha='left', fontsize=11, color=colors[i],
                fontweight='bold', clip_on=False)

    # inelastic-threshold guides at +/- Delta
    for s in (-1, 1):
        ax.axvline(s * Delta, color='0.5', ls='--', lw=0.9, alpha=0.7)
    ax.text(Delta, -0.6, r'$+\Delta$', ha='center', va='top',
            fontsize=12, color='0.35')
    ax.text(-Delta, -0.6, r'$-\Delta$', ha='center', va='top',
            fontsize=12, color='0.35')

    ax.set_xlabel(r'Bias voltage  $V$  (mV)', fontsize=13)
    ax.set_ylabel(r'$dI/dV$  (offset, a.u.)', fontsize=13)
    ax.tick_params(axis='x', labelsize=11)
    ax.set_title(r'Point-spectroscopy IETS waterfall — P-T-P polymer, $N=3$'
                 + f'\n(inelastic threshold  $|eV| = \\Delta = {Delta:.2f}$ meV,'
                 + f'  $\\Gamma = {GAMMA:.1f}$ meV)', fontsize=14)
    ax.set_xlim(V_MIN, V_MAX)
    ax.set_ylim(-1.2, n_sites * OFFSET + 1.6)
    ax.set_yticks([])
    ax.grid(axis='x', alpha=0.25)

    fig.tight_layout(rect=(0.0, 0.0, 0.80, 1.0))   # reserve right margin for labels
    out_png = os.path.join(out_dir, 'iets_waterfall_N3.png')
    fig.savefig(out_png, dpi=220, bbox_inches='tight')
    print(f"\nSaved waterfall plot to {out_png}")


if __name__ == '__main__':
    main()
