import numpy as np
from scipy.spatial import cKDTree
from units import make_triangle, acc

CC_LINK = 1.48

def rot(deg):
    r=np.radians(deg); c,s=np.cos(r),np.sin(r)
    return np.array([[c,-s],[s,c]])
def ang(v): return np.degrees(np.arctan2(v[1],v[0]))
def angdiff(a,b): return abs((a-b+180)%360-180)

class Unit:
    def __init__(self, kind, apex_up=True):
        P,adj,corners,outdir = make_triangle(kind)
        if not apex_up:
            P=P*np.array([1.0,-1.0]); outdir={c:outdir[c]*np.array([1.,-1.]) for c in corners}
        self.kind=kind; self.P=P; self.adj=adj
        self.corners=list(corners); self.outdir=dict(outdir)
    def apply(self,R,t):
        self.P=(R@self.P.T).T+t
        self.outdir={c:R@self.outdir[c] for c in self.corners}

def attach(prev_pos, prev_dir, unit, in_c):
    phi = ang(-prev_dir) - ang(unit.outdir[in_c])
    unit.apply(rot(phi), np.zeros(2))
    unit.P = unit.P + (prev_pos + CC_LINK*prev_dir - unit.P[in_c])

def pick_in(unit, prev_dir):
    # corner whose outdir best points back against prev_dir
    tgt = ang(-prev_dir)
    return min(unit.corners, key=lambda c: angdiff(ang(unit.outdir[c]), tgt))
def pick_out(unit, exclude, prefer_up):
    # among remaining corners choose the one most forward (+x); tie-break by y sign
    cs=[c for c in unit.corners if c!=exclude]
    def score(c):
        d=unit.outdir[c]
        return (-d[0], -(d[1] if prefer_up else -d[1]))
    return min(cs, key=score)

def build(n_monomers=3):
    """sequence of units: P T P | P T P | ... with apex alternating EVERY unit -> horizontal zigzag."""
    kinds=[]
    for m in range(n_monomers):
        kinds += ['P','T','P']
    units=[]
    # first unit apex up, pointing right
    apex=True
    u0=Unit(kinds[0], apex_up=apex); units.append(u0)
    rad_c0 = min(u0.corners, key=lambda c: u0.outdir[c][0])  # leftmost corner = radical
    out_c  = max(u0.corners, key=lambda c: u0.outdir[c][0])  # rightmost = forward
    cur_pos=u0.P[out_c]; cur_dir=u0.outdir[out_c]
    radicals=[(0,rad_c0)]
    prefer_up=True
    for k in range(1,len(kinds)):
        apex=not apex
        u=Unit(kinds[k], apex_up=apex)
        in_c=pick_in(u,cur_dir)
        attach(cur_pos,cur_dir,u,in_c)
        units.append(u)
        if k==len(kinds)-1:
            rad = max(u.corners, key=lambda c: u.outdir[c][0])  # rightmost = radical
            radicals.append((k,rad))
        else:
            out_c=pick_out(u,in_c,prefer_up); prefer_up=not prefer_up
            cur_pos=u.P[out_c]; cur_dir=u.outdir[out_c]
    return units, radicals

if __name__=='__main__':
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    units,rad=build(3)
    Ps=[];offs=[];off=0
    for u in units: Ps.append(u.P);offs.append(off);off+=len(u.P)
    P=np.vstack(Ps);nc=len(P)
    tree=cKDTree(P);adj=[set() for _ in range(nc)]
    for i,j in tree.query_pairs(r=1.6): adj[i].add(j);adj[j].add(i)
    fig,ax=plt.subplots(figsize=(22,7))
    ax.add_collection(LineCollection([[P[i],P[j]] for i in range(nc) for j in adj[i] if j>i],colors='#444',linewidths=1))
    ax.scatter(P[:,0],P[:,1],s=12,c='#222',zorder=3)
    ridx=[offs[k]+c for k,c in rad]
    ax.scatter(P[ridx,0],P[ridx,1],s=90,c='red',zorder=5)
    ax.set_aspect('equal');ax.set_title('v2 alternating-apex biphenyl chain')
    plt.tight_layout();plt.savefig('v2_test.png',dpi=130,bbox_inches='tight')
    print('nc',nc,'yspan',round(P[:,1].max()-P[:,1].min(),1),'xspan',round(P[:,0].max()-P[:,0].min(),1))
    clј=[(i,j) for i,j in tree.query_pairs(r=1.2)]
    print('very close pairs (<1.2):',len(clј))