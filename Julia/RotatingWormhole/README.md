# Rotating wormhole: Julia code

Run every script from this folder with `julia --project=.. <script>`.  `<root>` below is the
folder that contains `output-m1`, `output-m2`, `output-m3` (written by `QNMThreads.jl`), and
`results` is where the analysis scripts read and write their CSV files.

| Step | Script | What it does | Paper |
|---|---|---|---|
| 1 | `QNMThreads.jl` | Solves the quadratic pencil (M0 + w M1 + w² M2) b = 0 for the Maple matrices in `data-m<m>` at N_x = 90, 95, 100 (N_χ = 3), 220 digits, companion linearisation and generalised Schur decomposition; matches roots across the three resolutions and writes the raw spectra and the three-resolution diagnostics to `output-m<m>/J<k>`. Run once per sector: `QNM_M=1 julia -t 12 QNMThreads.jl`. | Sec. III B |
| 2 | `radial_pencils_arb.jl` | Solves the 99 decoupled radial pencils (Proposition III.1) at N_x = 90 in 256-bit ball arithmetic (Arb). Use `-t 4`. Output: `radial_spectra/`. | Secs. III A, IV A |
| 3 | `classify_modes.jl <root>` | Checks that the radial spectra reproduce the raw 2D spectra, labels roots by L, applies the refined criterion (Eq. refined_criterion), classifies the rejected roots, assigns overtone numbers. Output: `results/converged_qnms.csv`, `rejected_roots.csv`, `retained_roots.csv`; prints the counts of Sec. IV A and Table II. | Sec. IV A, Table II |
| 4 | `audit_symmetries.jl <root>` | Mirror symmetry of every raw spectrum and the shared families of runs with equal μ = mj, in the solver precision. | Sec. IV A |
| 5 | `validate_direct.jl` | Recomputes every converged mode with n ≤ 1 by direct integration and compares. Use `-t 4`. | Sec. IV A |
| 6 | `continue_mu.jl` | Continues the n = 0, 1 branches of L = 1..5 in μ by direct integration up to μ = 100. Use `-t 4`. Output: `results/continuation.csv`. | Figs. 1, 3, 4; Tables I, III |
| 7 | `limit_problem.jl` | Eigenvalues w_{L,n} of the large-μ limit problem and the light-ring estimate. Output: `results/limit_eigenvalues.csv`. | Sec. IV E, Table I |
| 8 | `static_parity.jl` | Parity of the static modes under η → −η. | Sec. IV B |
| 9 | `make_tables.jl` | LaTeX tables of the converged modes. Output: `results/appendix_tables.tex`. | Tables VIII–X |
| 10 | `make_figures.jl <root>` | Figures (PDF and PNG) in `plots/`. | Figs. 1–4 |

Supporting modules: `WormholeData.jl` (spin tags, file names, readers) and
`DirectIntegration.jl` (direct integration of the radial equation along complex rays with an
adaptive eighth-order Runge–Kutta scheme, and the large-μ limit problem).

Steps 5–8 need only `results/retained_roots.csv`; step 10 also reads the raw spectra.

The folder `results/` contains the output of steps 3, 5, 6 and 7 for the computations reported in
the paper.
