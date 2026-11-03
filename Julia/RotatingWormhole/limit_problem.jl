# =====================================================================================
#  limit_problem.jl
#
#  Eigenvalues w_{L,n} of the large-mu limit problem (Sec. IV E, Table I),
#
#      u_yy + [(w + 2/|y|^3)^2 - L(L+1)/y^2] u = 0,      y < 0,
#
#  with transmission into the throat region as y -> 0^- and outgoing behaviour as
#  y -> -inf, solved by complex-contour integration (DirectIntegration.jl).  Each root is
#  recomputed from a perturbed start to check its independence of the iteration.  The
#  light-ring estimate of Eq. (limit_eikonal),
#
#      w_{L,n} ~ (2/(3 sqrt 6)) l_L^(3/2) - i (n + 1/2) (l_L/3)^(1/2),   l_L = sqrt(L(L+1)),
#
#  and, if results/continuation.csv exists, 10 w_- at mu = 100 are printed alongside.
#
#  Usage:  julia --project=.. limit_problem.jl  [results]
#  Output: results/limit_eigenvalues.csv
# =====================================================================================

using Printf

include("DirectIntegration.jl")
using .DirectIntegration

# starting guesses (n = 0 and n = 1; the n = 1 branch of L = 1 reaches the branch cut)
const GUESSES = Dict((1, 0) => 0.45 - 0.37im, (2, 0) => 1.04 - 0.47im, (3, 0) => 1.75 - 0.55im,
                     (4, 0) => 2.57 - 0.62im, (5, 0) => 3.49 - 0.68im,
                     (2, 1) => 0.45 - 1.38im, (3, 1) => 1.28 - 1.62im, (4, 1) => 2.16 - 1.84im,
                     (5, 1) => 3.12 - 2.03im)

eikonal(L, n) = (l = sqrt(L * (L + 1)); 2 / (3sqrt(6)) * l^1.5 - 1im * (n + 0.5) * sqrt(l / 3))

function mu100(resdir)
    fn = joinpath(resdir, "continuation.csv")
    d = Dict{Tuple{Int,Int},ComplexF64}()
    isfile(fn) || return d
    for line in readlines(fn)[2:end]
        c = split(line, ',')
        parse(Float64, c[3]) == 100 || continue
        d[(parse(Int, c[1]), parse(Int, c[2]))] = complex(parse(Float64, c[4]), parse(Float64, c[5]))
    end
    return d
end

function main(resdir = "results")
    mkpath(resdir)
    C = mu100(resdir)
    open(joinpath(resdir, "limit_eigenvalues.csv"), "w") do io
        println(io, "L,n,Re_w,Im_w,restart_difference,Re_eikonal,Im_eikonal")
        for key in sort(collect(keys(GUESSES)); by = t -> (t[2], t[1]))
            L, n = key
            w, _ = qnm_limit(L, GUESSES[key])
            w2, _ = qnm_limit(L, w * 1.02)
            e = eikonal(L, n)
            @printf(io, "%d,%d,%.9f,%.9f,%.1e,%.6f,%.6f\n", L, n, real(w), imag(w), abs(w2 - w), real(e), imag(e))
            c = get(C, key, NaN + NaN * im)
            @printf("L=%d n=%d: w_{L,n} = %.6f%+.6fi  (restart %.0e)   10 w_-(mu=100) = %.6f%+.6fi   eikonal %.3f%+.3fi\n",
                    L, n, real(w), imag(w), abs(w2 - w), real(10c), imag(10c), real(e), imag(e))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? "results" : ARGS[1])
end
