# SCF diagnostics — idea_screen_1/soc_static_1 (vasp_ncl, single point)

Run: ~/agf2-b1-soc1/B1_g3.0/ (private WSL), start 2026-10-22 13:41:50 UTC.

## Selected log markers (vasp.stdout)

```
 running    8 mpi-ranks, on    1 nodes
 distrk:  each k-point on    8 cores,    1 groups
 distr:  one band on    4 cores,    2 groups
 WAVECAR not read
 entering main loop
```
64 electronic iterations (DAV 1..64). Convergence tail (DAV 59..64):
```
DAV:  59    -0.341952857712E+02   -0.13350E-05   -0.61287E-07 ...
DAV:  60    -0.341952863207E+02   -0.54947E-06   -0.55978E-07 ...
...
DAV:  64    -0.341952864792E+02   -0.96868E-06   -0.45746E-07 ...
   1 F= -.35550016E+02 E0= -.35550016E+02  d E =-.235045E-11  mag=     0.0000    -0.0010    -0.0378
```
OUTCAR markers: `aborting loop because EDIFF` (1), `General timing and accounting informations` (present, clean finish), TOTAL-FORCE block (1, final forces), elapsed 16769.979 s.

## Diagnostics summary

| item | value |
| --- | --- |
| electronic iterations to EDIFF=1e-6 | 64 |
| final dE (extrapolated) | -9.69e-7 eV |
| final d eps | -4.57e-8 |
| EDIFF reached in OUTCAR | yes (1 marker) |
| OSZICAR F= line | `1 F= -35.550016  E0= -35.550016  d E=-2.35e-11  mag= (0.0000, -0.0010, -0.0378)` |
| final TOTEN (vasprun e_fr_energy) | -35.55001581 eV |
| STOPCAR | never created |
| NELM reached | no (160 cap not hit) |
| convergence path | monotone energy descent to -35.195 eV (per-scstep sum), slow SCF tail (~4.2 min/iteration average; last ~40 min spent from dE~5e-6 to 1e-6) |
| max memory (rank0) | 928,708 kB |
| total CPU / user / system (sec) | 16353.2 / 15972.7 / 380.6 |

## Notes

- OSZICAR of this VASP 6.6.0 build carries the DAV electronic-iteration stream rather than the classic per-ionic-step records; the authoritative electronic record used here is OUTCAR + vasprun.xml.
- Convergence tail crossing the 4 h wall-clock cap is documented in RUN_RECEIPT.json (budget_overshoot_notice).
