# Maple worksheets

## RotatingWormhole

* `MyExploration.mw`: derivation of the compactified two-dimensional operator: the scalar
  wave equation on the global metric, the asymptotic factorisation (axis regularity and the
  interpolating outgoing phase), the compactification x = (2/π) arctan(η), and the coefficient
  functions L200R … L000I of the quadratic pencil (Sec. II).
* `MatrixGenerator.mw`: Chebyshev–Gauss collocation of that operator with N_x = 90, 95, 100
  radial and N_χ = 3 angular nodes, in 200-digit arithmetic, for the spin tags k = 0..10
  (j* = k/10, j* = 0.95 for k = 10, j = j*/(1 − j*)).  Set `m` (1, 2 or 3) in the second
  execution group; the matrices are written to `data/` as
  `M2-<Nx>-3-<k>.mat`, `M1R-…`, `M1I-…`, `M0R-…`, `M0I-…` (MATLAB ASCII), i.e. the pencil
  M0 + w M1 + w² M2 with M1 = M1R + i M1I, M0 = M0R + i M0I.  Copy or link `data/` to
  `Julia/RotatingWormhole/data-m<m>` for `QNMThreads.jl`.

## Kerr

* `MyExploration.mw`: derivation of the coefficient functions of the two-dimensional Kerr
  operator (s = 0, units 2M = 1) after the factorisation of Eq. (Kerr_full_factorisation),
  including the angular prefactor e^{aωχ}.
* `MatrixGenerator_Nu.mpl`: plain-text Maple script: collocation on Chebyshev–Gauss grids at
  (N_x, N_χ) = (60, 8), (70, 10), (80, 12) for a/M = 0, 0.1, …, 0.9, 0.99, 40 digits; writes
  `M2-<Nx>-<Nu>-<k>.mat`, … to `DATADIR` (default `data`).  Run with
  `restart; read "MatrixGenerator_Nu.mpl";` (about 1.5 min and 350 MB per spin).
