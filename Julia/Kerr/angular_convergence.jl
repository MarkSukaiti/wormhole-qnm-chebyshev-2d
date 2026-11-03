# =====================================================================================
#  angular_convergence.jl
#
#  Convergence of the two-dimensional Kerr spectrum with the number N_chi of angular
#  collocation nodes (Sec. V, Table VI), in double precision.
#
#  The operator is the one of MatrixGenerator_Nu.mpl (coefficients L200R ... L000I, units
#  2M = 1, s = 0), collocated on the same Chebyshev-Gauss grids in x and chi = cos(theta),
#  with the same ordering of unknowns; the quadratic pencil (M0 + w M1 + w^2 M2) b = 0 is
#  linearised as  [0 I; -M0 -M1] v = w [I 0; 0 M2] v.  For a/M = 0.9 and N_x = 60 the
#  relative deviation of the modes (l, n) = (1,0), (1,1), (2,0), (2,1) from the qnm values
#  is listed for N_chi = 2, 4, ..., 12.
#
#  Usage:  julia --project=.. angular_convergence.jl  [a/M] [N_x] [N_chi ...]
# =====================================================================================

using LinearAlgebra
using Printf

include("KerrData.jl")
using .KerrData

"Operator coefficients of MatrixGenerator_Nu.mpl, for spin a in units 2M = 1."
function coefficients(a, m = 1)
    k = sqrt(1 - 4a^2)
    K = (2a - 1)^3 * (2a + 1)^3
    return (
        L200R = (x, t) -> -256 * (a - 0.5)^2 * ((-2a^2 - x / 4 + 0.75) * k + a^2 * x - 3a^2 - x / 4 + 0.75) * (a + 0.5)^2,
        L111I = (x, t) -> -2 * (3k * x^2 + 16a^2 * x - 2x * k + 16a^2 + x^2 - k - 6x - 3) * K / k,
        L121R = (x, t) -> -8a * (t^2 - 1) * K,
        L100I = (x, t) -> -2 * (3x * k + 8a^2 - k + x - 3) * K / k,
        L100R = (x, t) -> 512 * (a - 0.5)^3 * (((-t - 1) * m - t) * k + m * (x - 1)) * (a + 0.5)^3 * a / k,
        L012R = (x, t) -> 2K * (1 + x) * (x - 1)^2,
        L011R = (x, t) -> 2K * (x - 1) * (3x + 1),
        L011I = (x, t) -> 4K * (x - 1)^2 * m * a / k,
        L022R = (x, t) -> -4K * (t^2 - 1),
        L021R = (x, t) -> -8K * (m + 1) * t,
        L000R = (x, t) -> -2K * (2m^2 + 2m - x + 1),
        L000I = (x, t) -> 4K * m * (x - 1) * a / k,
    )
end

"Values of T_p, T_p', T_p'' (p = 0..N-1) at the points, by recurrence."
function chebyshev_tables(pts, N)
    n = length(pts)
    T = zeros(n, N); T1 = zeros(n, N); T2 = zeros(n, N)
    for (i, x) in enumerate(pts)
        T[i, 1] = 1.0
        if N > 1
            T[i, 2] = x; T1[i, 2] = 1.0
        end
        for p in 2:N-1
            T[i, p+1] = 2x * T[i, p] - T[i, p-1]
            T1[i, p+1] = 2T[i, p] + 2x * T1[i, p] - T1[i, p-1]
            T2[i, p+1] = 4T1[i, p] + 2x * T2[i, p] - T2[i, p-1]
        end
    end
    return T, T1, T2
end

chebroots(N) = [cos((2i - 1) * pi / (2N)) for i in 1:N]

"""
    spectrum(a, Nx, Nu; m = 1) -> Vector{ComplexF64}

Finite eigenvalues 2M*omega of the collocated pencil (a in units 2M = 1, a = (a/M)/2).
Rows (i, j) -> (i-1) Nu + j, columns (p, q) -> p Nu + q + 1, as in the Maple generator.
"""
function spectrum(a, Nx, Nu; m = 1)
    L = coefficients(a, m)
    xs, ts = chebroots(Nx), chebroots(Nu)
    Tx, Tx1, Tx2 = chebyshev_tables(xs, Nx)
    Tt, Tt1, Tt2 = chebyshev_tables(ts, Nu)
    F, Fx, Fxx = kron(Tx, Tt), kron(Tx1, Tt), kron(Tx2, Tt)
    Ft, Ftt = kron(Tx, Tt1), kron(Tx, Tt2)
    X = repeat(xs; inner = Nu)
    TH = repeat(ts; outer = Nx)
    g(f) = f.(X, TH)
    M2 = g(L.L200R) .* F
    M1 = (g(L.L100R) .+ im .* g(L.L100I)) .* F .+ g(L.L121R) .* Ft .+ im .* g(L.L111I) .* Fx
    M0 = (g(L.L000R) .+ im .* g(L.L000I)) .* F .+ g(L.L021R) .* Ft .+ g(L.L022R) .* Ftt .+
         (g(L.L011R) .+ im .* g(L.L011I)) .* Fx .+ g(L.L012R) .* Fxx
    n = size(M0, 1)
    Z, I_ = zeros(n, n), Matrix{Float64}(I, n, n)
    ev = eigvals([Z I_; -M0 -M1], [I_ Z; Z M2])
    return filter(isfinite, ev)
end

function main(aM = 0.9, Nx = 60, Nus = [2, 4, 6, 8, 10, 12])
    ref = read_reference(joinpath(@__DIR__, "reference", "kerr_qnm_reference.csv"))
    modes = [(1, 0), (1, 1), (2, 0), (2, 1)]
    println("a/M = $aM, N_x = $Nx: relative deviation from qnm for (l,n) = $modes")
    for Nu in Nus
        ev = spectrum(aM / 2, Nx, Nu) ./ 2                         # M*omega
        errs = map(modes) do (l, n)
            r = ref[(l, 1, n, aM)]
            minimum(abs.(ev .- r)) / abs(r)
        end
        @printf("  N_chi = %2d: %s\n", Nu, join([@sprintf("%.1e", e) for e in errs], "  "))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        main()
    else
        main(parse(Float64, ARGS[1]), parse(Int, ARGS[2]), parse.(Int, ARGS[3:end]))
    end
end
