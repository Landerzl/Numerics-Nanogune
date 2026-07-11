import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

N_monomers = 3
J_PP, J_TP, J_T = 38.0, 40.0, -500.0
n = 4*N_monomers

I2=sp.identity(2,format='csr',dtype=complex)
sx=0.5*sp.csr_matrix([[0,1],[1,0]],dtype=complex)
sy=0.5*sp.csr_matrix([[0,-1j],[1j,0]],dtype=complex)
sz=0.5*sp.csr_matrix([[1,0],[0,-1]],dtype=complex)

def embed(op,i,n):
    m=[I2]*n; m[i]=op
    out=m[0]
    for x in m[1:]: out=sp.kron(out,x,format='csr')
    return out

Sx=[embed(sx,i,n) for i in range(n)]
Sy=[embed(sy,i,n) for i in range(n)]
Sz=[embed(sz,i,n) for i in range(n)]

def bonds(N):
    B=[]
    for m in range(N):
        b=4*m
        B+=[(b,b+1,J_TP),(b+1,b+2,J_T),(b+2,b+3,J_TP)]
        if m<N-1: B.append((b+3,b+4,J_PP))
    return B

H=sp.csr_matrix((2**n,2**n),dtype=complex)
for i,j,J in bonds(N_monomers):
    H=H+J*(Sx[i]@Sx[j]+Sy[i]@Sy[j]+Sz[i]@Sz[j])

Sztot=sum(Sz)
S2=sum(Sx)@sum(Sx)+sum(Sy)@sum(Sy)+sum(Sz)@sum(Sz)

k=12
ev,evec=eigsh(H,k=k,which='SA')
idx=np.argsort(ev); ev=ev[idx]; evec=evec[:,idx]

def texp(A,v): return np.real(np.vdot(v,A@v))
def totS(v): return 0.5*(np.sqrt(1+4*max(texp(S2,v),0))-1)

E0=ev[0]; gs=evec[:,0]
print(f"N={N_monomers}, sites={n}, dim={2**n}")
print(f"GS:  E0={E0:.4f} meV  Sz={texp(Sztot,gs):+.3f}  S={totS(gs):.3f}")
for a in range(1,6):
    print(f"  state {a}: dE={ev[a]-E0:8.4f} meV  Sz={texp(Sztot,evec[:,a]):+.3f}  S={totS(evec[:,a]):.3f}")

# spin gap = to first excited
gap=ev[1]-E0
print(f"Spin gap Delta = {gap:.4f} meV")

# ground-state local magnetization (expect ~0)
mz_gs=np.array([texp(Sz[i],gs) for i in range(n)])
print("max|<Sz_i>| in GS =",np.max(np.abs(mz_gs)))

# lowest triplet -> Sz=+1 member (reveals texture)
S_all=np.array([totS(evec[:,a]) for a in range(k)])
trip=[a for a in range(k) if S_all[a]>0.9]
lt=trip[0]
Elt=ev[lt]
mult=[a for a in range(k) if abs(ev[a]-Elt)<1e-3]
sub=[evec[:,a] for a in mult]
M=np.array([[np.vdot(sub[a],Sztot@sub[b]) for b in range(len(sub))] for a in range(len(sub))])
w,Uu=np.linalg.eigh(M)
pol=sub[0]*0
c=Uu[:,np.argmin(np.abs(w-1.0))]
for a in range(len(sub)): pol=pol+c[a]*sub[a]
pol=pol/np.linalg.norm(pol)
mz_pol=np.array([texp(Sz[i],pol) for i in range(n)])
print("lowest triplet Sz=+1: sum<Sz>=",round(mz_pol.sum(),3))
np.savez('ed_result.npz',mz_gs=mz_gs,mz_pol=mz_pol,E0=E0,gap=gap,n=n,N=N_monomers)
