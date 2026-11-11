# Stage 7c: a physical mechanism for the sliding-induced static NMEE reversal

**Scope:** the established synthetic, 8-band Kane–Mele homobilayer, gapped T=0 equilibrium/static finite-Zeeman response. This is NOT a correction of driven-DC NMEE, not a demonstrated ferroelectric, and not a DFT material claim.

## Microscopic test of the two orbital hopping channels

In a declared A/B orbital basis, $H=H_0+sV_{\rm same}+cV_{\rm cross}$. The layer-sign unitary gives $\Gamma(s,c)=\Gamma(-s,-c)$, and disconnected mirror-related layers give $\Gamma(0,0)=0$. Under analyticity,
\[
\Gamma_{xz;x}(u;s,c)=A(u)s^2+B(u)c^2+C(u)sc+O(t^4).
\]
$A,B,C$ have no universal signs. The finite-$t$ odd cross-channel coefficient $C_{\rm odd}=[\Gamma(t,t)-\Gamma(t,-t)]/(2t^2)$ tends to $C$ at weak coupling; $E_{\rm even}=[\Gamma(t,t)+\Gamma(t,-t)]/(2t^2)$ tends to $A+B$. At finite $t$ they include $t^4$ corrections. Orbital-channel decomposition is basis dependent.

For $t=5$ meV, radius 0.75 nm$^{-1}$ *K/K'* disks, 26 radial Gauss ×72 angular points, ±0.25 T, original program units divided by $t^2$:

| u | $C_{\rm odd}$ | $E_{\rm even}$ | Full | low→low | involving remote |
|--:|--:|--:|--:|--:|--:|
| .10 | -.031093 | +.003542 | -.027551 | -.023096 | -.004455 |
| .20 | -.058611 | +.026669 | -.031942 | -.102296 | +.070354 |
| .30 | -.079632 | +.083471 | +.003839 | -.274240 | +.278079 |
| .40 | -.092349 | +.179022 | +.086673 | -.545952 | +.632626 |

**Physical reading:** the AA/BB–AB/BA interference is negative, but AB/BA on its own produces a positive and more rapidly rising response. The AA/BB-only response is nearly negligible ($A\approx8.5\times10^{-6}$ versus $B\approx0.00314$ at $u=.1$ and $t=2.5$ meV). By $u=.4$, the exact low-four-only response remains negative, while virtual transitions involving remote occupied/empty bands are positive and larger. A positive $B$ includes important low-occupied→remote-empty processes: at $u=.4$, $t=2.5$ meV, cross-only LL/RL/LR/RR coefficients are approximately $(-.4169,+.1999,+.6054,-.2208)$.

## Symmetry versus an accidental coefficient hierarchy

The *tensor trace* is constrained by $C_3+m_z$ to start no lower than cubic in sliding displacement, while the nematic tensor can be linear. This **does not** force the individual $\Gamma_{xz;x}$, its mixed-channel coefficient, or the even sector to be cubic.

Weak $t=2.5$ meV disk values for $u\leq.1$ fit the odd-polynomial hierarchy
\[
C_{\rm odd}(u)\simeq -.30810u+.55467u^3-.42012u^5,\qquad
E_{\rm even}(u)\simeq .0007642u+3.20566u^3-3.66006u^5.
\]
The EVEN coefficient has a **nonzero, symmetry-allowed linear term**. Its linear term happens to be about 400 times smaller than the mixed negative linear term in this model; the apparent fast-cubic positive growth is a numerical microscopic coefficient hierarchy, not a universal order-selection rule.

The crossing is therefore physically understandable as $C_{\rm odd}(u)+E_{\rm even}(u)=0$ at $s=c=t$, with a leading-order estimate $u_c\sim\sqrt{.308/3.2}\sim.31$ (not an exact root).

## Integration and checks

The two-valley disk zero is $u_c\simeq.29332$ at $t=5$ meV, radius .75 nm$^{-1}$. The roots change to .29366 at radius .5; to .30030 at $t=2.5$ meV; and .26798 at $t=10$ meV. Disk quadrature 28×72 versus 40×96 changes the $t=5$ meV, radius .75 root by <$10^{-7}$; varying Zeeman step ±.125/.25/.5 T moves it by several $10^{-6}$. The **full BZ** is checked separately with two accurate disks plus an exterior reciprocal-cell mesh. At $u=.28$ the resulting total stays negative (~$-1.73\times10^{-7}$), and at $u=.31$ positive (~$+2.52\times10^{-7}$). Exterior terms near the root are $\sim2\times10^{-9}$, so the disk-only root is not exactly the complete BZ root. No continuum error bound is claimed. A naive 91² *uniform* BZ grid misses narrow valley hotspots and yields grossly inaccurate conclusions.

For seven tested parameter settings (baseline, Rashba ±20%, intrinsic SOC ±20%, hopping-decay length ±12.5%), at $u=.1,.4$ the total response changes from negative to positive; $C_{\rm odd}<0$, cross-only response positive, chemical-potential gap sampled positive. These restricted cases are not a universal sign theorem.

## What an analytic reduction must preserve

At each k, the static coefficient uses *all* occupied→empty virtual transitions,
\[
\chi_{x;x}=-2\mathrm{Im}\sum_{v,c}
\frac{(S_x)_{vc}(v_x)_{cv}}{(\varepsilon_v-\varepsilon_c)^2},\quad
\Gamma_{xz;x}=\int\frac{d^2k}{(2\pi)^2}\,\partial_{B_z}\chi_{x;x}.
\]
Exact spectral decomposition low→low, remote occupied→low empty, low occupied→remote empty, and remote→remote reconstructs the full response. The remote-to-remote background must not be dropped. Occupied projectors and contour resolvents give a formulation stable under internal degeneracy. If one integrates out the remote subspace, both the source-differentiated self energy and the remote occupied-band background/contact pieces must remain; a 4-band Hamiltonian or even a covariant 4-band velocity alone cannot reproduce the numerical findings.

**Next physics question:** use the specific valley symmetry-adapted orbital matrix elements to establish a testable model-specific sufficient condition for $C<0$ and $B>0$, with an error-controlled virtual-excitation amplitude. Keep eight bands until an effective response theory is validated. No DFT.

See the independent Stage 7c research archive for full numeric CSV, small-u and hopping-strength scans, adaptive BZ tests, parameter sensitivity, and figures. The earlier Stage 7 and Stage 7b reproductions remain unchanged.
