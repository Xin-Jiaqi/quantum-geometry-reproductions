"""Independent contraction of the printed SM(9) diagonal term versus SM(11).

No production kernel imports. This checks a local algebraic identity only;
it does not select a correction to the paper's full DC response.
"""
import json
import numpy as np


def audit():
    rng = np.random.default_rng(1101)
    rows = []
    for bands in (2, 3, 5):
        def hermitian():
            a = rng.normal(size=(bands,bands)) + 1j*rng.normal(size=(bands,bands))
            return (a+a.conj().T)/2
        r, spin, obs = hermitian(), hermitian(), hermitian()
        e = np.arange(bands)*.7
        f = np.linspace(.9,.1,bands)
        Q = np.real(r*spin.T)
        diagonal_rho9 = np.zeros(bands)
        coefficient = 0.
        for n in range(bands):
            for m in range(bands):
                if n == m:
                    continue
                diagonal_rho9[n] -= 2*(f[n]-f[m])*Q[n,m]/(e[n]-e[m])**2
                coefficient += 2*f[n]*(obs[n,n].real-obs[m,m].real)*Q[n,m]/(e[n]-e[m])**2
        # Strip positive g*mu_B^2; Bbar/B=-g*mu_B supplies the minus.
        contracted = -np.dot(obs.diagonal().real, diagonal_rho9)
        printed11 = -coefficient
        assert abs(contracted-coefficient) < 1e-12
        assert abs(contracted-printed11) > 1e-5
        rows.append(dict(bands=bands,eq9_contracted=float(contracted),
                         plus_reindexed=float(coefficient),eq11_printed_diagonal=float(printed11)))
    return dict(source='Jia PRL136226302 SM page3 Eqs9,11; Bbar=-g mu_B B on page2',
                sm_sha256='13bada2e109d9f1faea6be9370c50bec04968b06057c262cd7ef69e29a28673a',
                conclusion='INTERNAL_DIAGONAL_SIGN_INCONSISTENCY', checks=rows,
                scope='Only the displayed diagonal term; no claim that either full equation is corrected')


if __name__ == '__main__':
    print(json.dumps(audit(),indent=2))
