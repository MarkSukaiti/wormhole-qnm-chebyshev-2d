# Kerr benchmark: Julia code

Run every script from this folder with `julia --project=.. <script>`.  `<output>` below is the
folder with `J0` … `J10` written by `QNMThreads_Nu.jl`.

| Step | Script | What it does | Paper |
|---|---|---|---|
| 1 | `QNMThreads_Nu.jl` | Solves the quadratic pencil for the matrices of `Maple/Kerr/MatrixGenerator_Nu.mpl` (in `data/`) at (N_r, N_χ) = (60, 8), (70, 10), (80, 12) and matches the roots across the three resolutions (tolerance 1e-6). `ARITHMETIC` selects Float64, Double64 (used in the paper) or BigFloat. `julia -t 12 QNMThreads_Nu.jl`. | Sec. V |
| 2 | `compare_with_qnm.jl <output>` | Identifies every converged root with a qnm mode (l ≤ 9, m = ±1), and prints Tables IV and V, the deviations of the unconverged modes at each resolution, and the comparison of δ with d_max. Output: `results/kerr_converged_identified.csv`, `results/kerr_tables.tex`. | Sec. V A, Tables IV, V |
| 3 | `angular_convergence.jl` | Double-precision rebuild of the same 2D operator; deviation from qnm against N_χ at a/M = 0.9, N_r = 60. | Sec. V A, Table VI |
| 4 | `literature_comparison.jl <output>` | Comparison with every published value in `reference/published_values.csv`: deviations, digits reproduced, rounding analysis, LaTeX rows. Output: `results/kerr_published_comparison.csv`, `results/literature_tables.tex`. | Sec. V B, Tables VII, XI, XII |

`KerrData.jl` holds the spin tags, file names and readers shared by these scripts.

## Reference data

* `reference/kerr_qnm_reference.csv` — continued-fraction values M·ω from the
  [qnm](https://github.com/duetosymmetry/qnm) package (Stein 2019) for s = 0, l = 1..9,
  m = ±1, n = 0..11, at the 11 benchmark spins; regenerate with `reference/qnm_reference.py`
  (Python, needs `pip install qnm`).
* `reference/published_values.csv` — the published values compared in Sec. V B, with the
  digits as quoted and the units of each source (Lin et al. 2017 tabulate in units 2M = 1).
