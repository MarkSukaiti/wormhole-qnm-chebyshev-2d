# =====================================================================================
#  continue_mu.jl
#
#  Continuation in mu = m j of the n = 0 and n = 1 branches of the multipole families
#  L = 1..5 by direct integration (DirectIntegration.jl), anchored on the converged
#  spectral values and extended beyond the resolved range (Secs. IV D-F, Figs. 3 and 4b,
#  Tables I and III), up to mu = 100.  Each anchored point is recomputed from the spectral value; beyond
#  the last anchor the starting guess extrapolates  mu^(1/2) w_- = a + b/mu.
#
#  Usage:  julia --project=.. -t 4 continue_mu.jl  [results]
#  Input:  results/retained_roots.csv
#  Output: results/continuation.csv  (L, n, mu, w_- from direct integration, and the
#          spectral anchor where there is one)
# =====================================================================================

using Printf
using Base.Threads

include("WormholeData.jl")
include("DirectIntegration.jl")
using .WormholeData
using .DirectIntegration

# beyond mu ~ 150 the start-up values at the first end underflow in double precision
const EXTENSION = [19, 27, 38, 57, 100]

function anchors(R, L, n)
    A = Dict{Float64,ComplexF64}()
    for o in R
        if o.converged && o.L == L && o.n == n && o.q < 1e-6
            A[round(o.mu; digits = 12)] = o.wminus
        end
    end
    return A
end

function run_branch(R, L, n)
    A = anchors(R, L, n)
    isempty(A) && return NamedTuple[]
    ladder = sort(unique(vcat(collect(keys(A)), [mu for mu in EXTENSION if mu > maximum(keys(A))])))
    out = NamedTuple[]
    for mu in ladder
        g = if haskey(A, mu)
            A[mu]
        else
            m1, m2 = out[end-1].mu, out[end].mu
            z1, z2 = sqrt(m1) * out[end-1].wminus, sqrt(m2) * out[end].wminus
            b = (z2 - z1) / (1 / m2 - 1 / m1)
            a = z2 - b / m2
            (a + b / mu) / sqrt(mu)
        end
        w, step = try
            qnm(mu, L, g + 3pi * mu)
        catch err
            @printf("L=%d n=%d mu=%g: no root (%s); branch stopped\n", L, n, mu, typeof(err))
            break
        end
        wm = w - 3pi * mu
        if !isfinite(wm)
            @printf("L=%d n=%d mu=%g: integration failed; branch stopped\n", L, n, mu)
            break
        end
        if step > 1e-6 * max(1.0, abs(w))       # secant did not converge: the branch is lost
            @printf("L=%d n=%d mu=%g: no convergence (last step %.1e); branch stopped\n", L, n, mu, step)
            break
        end
        anchor = get(A, mu, NaN + NaN * im)
        push!(out, (L = L, n = n, mu = mu, wminus = wm, step = step, anchor = anchor))
        @printf("L=%d n=%d mu=%9.4f  w_- = %.9f%+.9fi  mu^1/2 w_- = %.6f%+.6fi  |direct - spectral| = %.1e\n",
                L, n, mu, real(wm), imag(wm), real(sqrt(mu) * wm), imag(sqrt(mu) * wm), abs(wm - anchor))
    end
    return out
end

function main(resdir = "results")
    R = read_retained(joinpath(resdir, "retained_roots.csv"))
    jobs = [(L, n) for n in (0, 1) for L in 1:5]
    res = Vector{Vector{NamedTuple}}(undef, length(jobs))
    @threads for i in eachindex(jobs)
        res[i] = run_branch(R, jobs[i]...)
    end
    open(joinpath(resdir, "continuation.csv"), "w") do io
        println(io, "L,n,mu,Re_omega_minus,Im_omega_minus,Re_anchor_spectral,Im_anchor_spectral,last_secant_step")
        for b in res, r in b
            @printf(io, "%d,%d,%.15g,%.12f,%.12f,%.12f,%.12f,%.1e\n", r.L, r.n, r.mu, real(r.wminus),
                    imag(r.wminus), real(r.anchor), imag(r.anchor), r.step)
        end
    end
    println("written ", joinpath(resdir, "continuation.csv"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? "results" : ARGS[1])
end
