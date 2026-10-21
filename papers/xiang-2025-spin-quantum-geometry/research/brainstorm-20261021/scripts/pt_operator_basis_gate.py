"""Exhaustive PT-even/odd Hermitian 4x4 basis test, independent research only.

Basis: one orbital per layer, physical spin, PT = tau_x (i sigma_y) K
at fixed k. This representation is NOT a universal theorem for all bilayers.
"""
import numpy as np

I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1,-1]).astype(complex)
P=[I,X,Y,Z]
labels='0xyz'
U=np.kron(X,1j*Y)

def main():
    assert np.linalg.norm(U@U.conj()+np.eye(4))<1e-13
    allowed=[];forbidden=[]
    for a in range(4):
      for b in range(4):
        A=np.kron(P[a],P[b])
        AA=U@A.conj()@U.conj().T
        coeff=(np.trace(A@AA)/4).real
        parity=int(round(coeff))
        assert parity in (-1,1)
        assert np.linalg.norm(AA-parity*A)<1e-13
        (allowed if parity==1 else forbidden).append('tau_'+labels[a]+' sigma_'+labels[b])
    assert len(allowed)==6 and len(forbidden)==10
    print('PT-even basis ('+str(len(allowed))+'):',', '.join(allowed))
    print('PT-odd basis ('+str(len(forbidden))+'):',', '.join(forbidden))
    assert 'tau_x sigma_x' in forbidden and 'tau_y sigma_z' in forbidden
    print('PASS: no Hermitian PT-even spin-dependent tau_x/tau_y terms in this 4-band representation')

if __name__=='__main__':main()
