import numpy as np
from scipy.spatial import cKDTree

acc = 1.42
a1 = np.array([acc*np.sqrt(3.0), 0.0])
a2 = np.array([acc*np.sqrt(3.0)/2.0, acc*1.5])
delta_B = np.array([0.0, acc])

def site_pos(n1,n2,sub):
    p = n1*a1+n2*a2
    if sub=='B': p=p+delta_B
    return np.array([p[0],p[1]])
def ring_sites(n1,n2):
    return [(n1,n2,'A'),(n1,n2,'B'),(n1,n2+1,'A'),(n1+1,n2,'B'),(n1+1,n2,'A'),(n1+1,n2-1,'B')]

def make_triangle(order):
    order = {'P':2,'T':3}.get(order, order)
    """order=3 -> [3]triangulene (22 C); order=2 -> phenalenyl (13 C). Apex up, centered at centroid."""
    if order==3:
        rings=[(0,0),(1,0),(2,0),(0,1),(1,1),(0,2)]
    elif order==2:
        rings=[(0,0),(1,0),(0,1)]
    else:
        raise ValueError
    sites = sorted(set(s for r in rings for s in ring_sites(*r)))
    P = np.array([site_pos(*s) for s in sites])
    P = P - P.mean(axis=0)          # center at centroid
    tree=cKDTree(P); adj=[set() for _ in range(len(P))]
    for i,j in tree.query_pairs(r=acc*1.05): adj[i].add(j); adj[j].add(i)
    # 3 corner atoms = the 3 furthest from centroid, and each is a local maximum of radius
    r = np.linalg.norm(P, axis=1)
    order_idx = np.argsort(-r)
    corners=[]
    for idx in order_idx:
        if all(np.linalg.norm(P[idx]-P[c])>acc*2 for c in corners):
            corners.append(idx)
        if len(corners)==3: break
    # outward direction at each corner
    outdir = {c: P[c]/np.linalg.norm(P[c]) for c in corners}
    return P, adj, corners, outdir

for order in (2,3):
    P,adj,corners,outdir = make_triangle(order)
    print(f'order={order}: {len(P)} C, corners at angles',
          [round(np.degrees(np.arctan2(*outdir[c][::-1])),1) for c in corners])