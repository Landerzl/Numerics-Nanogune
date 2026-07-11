"""
Oligomero (3 monomeros) [fenalenilo-trianguleno-fenalenilo] enlazados por
enlaces simples C-C (tipo bifenilo). TODOS los bordes pasivados con H.
Molecula cerrada, sin radicales.
Requiere: numpy, scipy, matplotlib, sisl
"""
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import sisl
from build_v2 import build

acc, ach = 1.42, 1.09
N_MONOMERS = 1

units, _ = build(N_MONOMERS)

# merge units
Ps=[]; offs=[]; off=0
for u in units:
    Ps.append(u.P); offs.append(off); off+=len(u.P)
carbon = np.vstack(Ps); nc = len(carbon)

# global rotation -> horizontal
xc = carbon - carbon.mean(axis=0)
u_,s_,vt_ = np.linalg.svd(xc, full_matrices=False)
axis = vt_[0]
theta = -np.arctan2(axis[1], axis[0])
Rg = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
carbon = (Rg@xc.T).T
carbon3 = np.hstack([carbon, np.zeros((nc,1))])

# connectivity
tree = cKDTree(carbon)
adj = [set() for _ in range(nc)]
for i,j in tree.query_pairs(r=1.6):
    adj[i].add(j); adj[j].add(i)

# NO radicals: passivate ALL edge carbons
radical_idxs = set()   # empty!

h_coords=[]
for i in range(nc):
    n_cc = len(adj[i])
    if n_cc >= 3: continue
    bv = [(carbon3[j]-carbon3[i])/np.linalg.norm(carbon3[j]-carbon3[i]) for j in adj[i]]
    if n_cc == 2:
        avg = bv[0]+bv[1]; nrm=np.linalg.norm(avg)
        if nrm < 1e-8:
            perp = np.array([-bv[0][1],bv[0][0],0.0])
            h_coords.append(carbon3[i]+ach*perp)
        else:
            h_coords.append(carbon3[i]-ach*avg/nrm)
    elif n_cc == 1:
        away=-bv[0]; perp=np.array([-away[1],away[0],0.0])
        for sign in (1,-1):
            d = away*np.cos(np.radians(30)) + sign*perp*np.sin(np.radians(30))
            h_coords.append(carbon3[i]+ach*d/np.linalg.norm(d))

h_coords = np.array(h_coords) if h_coords else np.empty((0,3))
nh = len(h_coords)

all_xyz = np.vstack([carbon3, h_coords])
atom_labels = ['C']*nc + ['H']*nh
pad=15.0
mins=all_xyz.min(axis=0)-pad; maxs=all_xyz.max(axis=0)+pad
box=maxs-mins; box[2]=max(box[2],20.0)
all_xyz = all_xyz - mins

lattice = sisl.Lattice(box)
geom = sisl.Geometry(all_xyz, atoms=atom_labels, lattice=lattice)
geom.write('oligomer_closed.xyz')
print(f"C: {nc}, H: {nh}, Total: {nc+nh}")
print("Wrote oligomer_closed.xyz")

# plot
fig, ax = plt.subplots(figsize=(26,7))
bond_lines=[[all_xyz[i,:2],all_xyz[j,:2]] for i in range(nc) for j in adj[i] if j>i]
for i in range(nc):
    if len(adj[i])<3:
        for hi in range(nc, nc+nh):
            if np.linalg.norm(all_xyz[i]-all_xyz[hi]) < ach*1.15:
                bond_lines.append([all_xyz[i,:2], all_xyz[hi,:2]])
ax.add_collection(LineCollection(bond_lines, colors='#444444', linewidths=1.0, zorder=1))
ax.scatter(all_xyz[:nc,0], all_xyz[:nc,1], s=45, c='#222222',
           edgecolors='k', linewidths=0.3, zorder=3, label='C')
if nh>0:
    ax.scatter(all_xyz[nc:,0], all_xyz[nc:,1], s=16, c='#dddddd',
               edgecolors='#999999', linewidths=0.3, zorder=2, label='H')
ax.set_aspect('equal'); ax.legend(fontsize=11)
ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
ax.set_title('Oligomer (3 monomers) - fully passivated, no radicals')
ax.grid(True, alpha=0.15)
plt.tight_layout(); plt.savefig('oligomer_closed.png', dpi=190, bbox_inches='tight')
print("Plot saved.")