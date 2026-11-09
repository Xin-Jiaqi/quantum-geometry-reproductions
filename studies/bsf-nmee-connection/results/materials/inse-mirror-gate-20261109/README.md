# InSe mirror partner: electrical tuning and structural switching

Independent calculation is in progress; see `jobs.json`. Pre-result prediction
and numerical scope are sealed in [PLAN.md](PLAN.md). The primary observable
is the existing fixed-ion insulating spin coefficient C_l in s_x,l=C_l Ex Bz.
The perpendicular gate F is distinct from the electric and magnetic probes.

## Question and prediction

In the original A domain, a small positive F weakens both opposite layer-spin
responses. In B, the same positive F is predicted to strengthen both. Neither
claim means that a layer's spin changes its own sign. The mirror law relates
C_lower(B,F)=-C_upper(A,-F); the finite-field B values must be calculated
independently. Two meshes test tuning directions, without pursuing identical
enhancement factors. Five new SCFs are planned; no ionic switching is simulated.

## Can the gate also move the structure?

Yes in principle, but a change of electronic response does not itself establish
sliding, and a fixed-ion calculation cannot establish that a real device stays
in its original domain. There are two different derivatives: the change of C
with F at fixed atomic coordinates, and the additional change when coordinates
or domain fractions move with F. This study computes the first.

The two requested papers sharpen the distinction:

- [Ke, Liu and Liu, PRL135,046201 (2025)](https://doi.org/10.1103/jhlq-2dd7):
  in ideal C3-symmetric h-BN AB/BA domain interiors, a perpendicular field
  cannot produce a nonzero net in-plane layer force. Their simulations show
  switching through symmetry-breaking, wide domain walls, where off-diagonal
  Born-charge components provide lateral driving forces. These h-BN wall
  speeds and fields are not InSe switching thresholds.
- [Wang et al., PRL137,126801 (2026)](https://doi.org/10.1103/jq6y-y9ft):
  uniaxial strain breaks C3 in monodomain h-BN and allows directed lateral
  forces. Compressive/tensile strain changes the force direction and the
  accessible barrier; merely allowing a force does not guarantee switching.
  Their uniaxial strain is different from our previous C3-preserving biaxial
  strain. No new strain calculation is started here.

The force definition we use is unambiguous: the relevant derivative is
partial F_y / partial E_z = partial p_z / partial u_y (with a consistent
charge/unit convention). Published BEC index order is not silently copied.
Zero net lateral force at a high-symmetry point is not proof of positive
curvature, absence of thermal nucleation, or a finite experimental coercive
field. A field may change stiffness or barriers even when the initial force
vanishes. Interfaces, defects, edges and domain walls can remove the ideal
symmetry restriction.

Practical implication: distinguish a small, reversible readout within a
stationary domain from a pulse that changes the domain. Track an independent
A/B structural/domain marker during a gate scan; net spin signal or hysteresis
alone does not uniquely identify electronic tuning versus domain motion.
If C_eff=f C_A+(1-f) C_B, its gate derivative contains the additional term
(C_A-C_B) df/dF. At zero gate in ideal mirror partners, G_A=G_B, so a fixed
mixture has the same net linear electronic slope; moving domain fractions
can still add a separate signal. Ionic relaxation contributes further terms.

The current finite-field SCFs retain force tables for a limited consistency
check. They do not measure a full Born-charge tensor, a sliding barrier, or
domain-wall mobility. A future minimal mechanical check would examine
near-registry restoring forces/stiffness under the proposed readout field;
domain-wall switching would require a separate treatment.

## FQFE connection and scope

Internal sublattice switching within an InSe monolayer is a different
coordinate from relative sliding of two intact monolayers. A fractional
polarization branch does not exempt ordinary response tensors from point-group
selection rules. In an isolated D3h layer at zero gate, the presently probed
Gamma_xxz (s_x from Ex Bz) remains forbidden. Some all-in-plane EB components
are allowed; a space-inverted pair reverses these allowed coefficients.
Two axial factors (spin and magnetic field) cancel under spatial operations,
so this NMEE tensor has the same spatial transformation character as a polar
rank-three tensor, without implying identical microscopic physics to BPVE.

`symmetry_checks.json` verifies the allowed indices in the actual monolayer
Cartesian frame. No FQFE switching or additional magnetic-field-direction
calculation is launched. The user's manuscript was read locally; its PDF,
full text, numerical results and figures are not redistributed in this archive.

## Reading scope and provenance

The two published main texts were read from Zotero. Relevant complete pages,
including the domain-wall paper's End Matter, were visually inspected.
The directional paper's SM computational-method/strain-control sections were
accessible through the publisher's indexed PDF; later retrieval was blocked.
The domain-wall SM was not obtained in this pass. No unseen supplemental
table is used to claim an InSe switching threshold.

Scientific computation is restricted to PBS compute nodes. No other shared
account job is altered. Licensed PAW data and wavefunction caches remain on
compute storage; only the selected inputs, code, small results and our own
analysis are archived. Local literature extractions are not committed.
