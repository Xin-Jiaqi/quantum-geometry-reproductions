# Source audit and optional author question (not sent)

Jia et al., PRL **136**, 226302 (2026),
[DOI](https://doi.org/10.1103/r3x2-kb93).
The local published SM was visually checked at pp. 2–3; SHA256:
`13bada2e109d9f1faea6be9370c50bec04968b06057c262cd7ef69e29a28673a`.

## Confirmed local inconsistency (2026-11-01)

Page 2 defines Bbar = -g mu_B B and page 3 defines
s_alpha = mu_B Tr(sigma_alpha rho). The displayed diagonal term of Eq. (9) is

```
rho_nn / (E Bbar) = -2 sum_l (f_n-f_l) Q_nl / epsilon_nl^2.
```

For Hermitian r and spin, Q_nl=Q_ln. Contracting this term using the stated
magnetic-field convention gives

```
Gamma_diagonal = +2 g mu_B^2 sum_nl f_n (sigma_nn-alpha - sigma_ll-alpha)
                                            Q_nl / epsilon_nl^2.
```

Eq. (11) prints **minus** in front of the same reindexed sum. Hence these two
displayed diagonal terms cannot both follow from the definitions. This
comparison requires no frequency limit, Fourier normalization, fitting or
choice of eigenvector phases. It does not determine which full equation
should be corrected or rule out a redistribution among other terms.

`validation/check_diagonal_source.py` independently contracts the literal
Eq. (9) term and compares the Eq. (11) term for 2/3/5-band matrices. Residuals
against the plus-sign reindexing are below 1e-12; the printed term is opposite.
`results/diagonal_source_audit.json` records the numbers.

An additional gapped two-band test in `tests/test_dynamic.py` has a diagonal
observable at k=0, so the printed Eq. (13)/(14) kernels vanish at that point.
The full finite-frequency orbital response tends to a nonzero diagonal
contribution as the two switching rates decrease together. This is a local
model check, not a proof of the full integrated metallic DC response.

## Still open

The previous question about the relation between Eq. (9), Eq. (13)/(14),
and the nonequilibrium DC preparation remains open, particularly for exact
PT doublets at finite T. The new dynamic interface evaluates SM (6)/(7)
before that split; it does not silently choose a correction to Fig. 1(c).

If author discussion is later desired, the concrete question is:

> Using Bbar=-g mu_B B, the displayed diagonal term in Eq. (9) appears to give
> a plus-sign contribution after contraction, whereas Eq. (11) prints a minus.
> Could you clarify this sign and how that term enters Eqs. (13)/(14), including
> the field-limit and degenerate-subspace prescription used for Fig. 1(c)?

No correspondence has been sent.
