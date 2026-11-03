# =====================================================================================
#  validate_direct.jl
#
#  Compares every converged spectral quasinormal frequency with Re varpi >= 0 and overtone
#  number n <= 1 with the value obtained by direct integration of the radial equation
#  (DirectIntegration.jl), started from the spectral value (Sec. IV A: 166 modes, median
#  deviation ~1e-12, all within the convergence diagnostic d_max).
#
#  Usage:  julia --project=.. -t 4 validate_direct.jl  [results]
#  Input:  results/retained_roots.csv (classify_modes.jl)
#  Output: results/direct_validation.csv and a summary on stdout.
# =====================================================================================

using Printf
using Statistics
using Base.Threads

include("WormholeData.jl")
include("DirectIntegration.jl")
using .WormholeData
using .DirectIntegration

function main(resdir = "results")
    R = read_retained(joinpath(resdir, "retained_roots.csv"))
    jobs = [o for o in R if o.converged && real(o.varpi) >= -1e-9 && 0 <= o.n <= 1]
    out = Vector{Any}(undef, length(jobs))
    lk = ReentrantLock()
    @threads for i in eachindex(jobs)
        o = jobs[i]
        t0 = time()
        wd, step = try
            qnm(o.mu, o.L, o.w)
        catch err
            (NaN + NaN * im, NaN)
        end
        out[i] = (o, wd, step)
        lock(lk) do
            @printf("m=%d k=%2d L=%d n=%d  spectral %s  direct %s  |diff| %.1e  d_max %.1e  (%.1f s)\n",
                    o.m, o.k, o.L, o.n, string(o.w), string(wd), abs(wd - o.w), o.dmax, time() - t0)
        end
    end
    open(joinpath(resdir, "direct_validation.csv"), "w") do io
        println(io, "m,k,mu,L,n,Re_omega_spectral,Im_omega_spectral,Re_omega_direct,Im_omega_direct,abs_diff,d_max")
        for (o, wd, _) in out
            @printf(io, "%d,%d,%.15g,%d,%d,%.15f,%.15f,%.15f,%.15f,%.3e,%.3e\n", o.m, o.k, o.mu,
                    o.L, o.n, real(o.w), imag(o.w), real(wd), imag(wd), abs(wd - o.w), o.dmax)
        end
    end
    ok = [(abs(wd - o.w), o.dmax) for (o, wd, _) in out if isfinite(abs(wd))]
    d = first.(ok)
    @printf("modes compared: %d of %d\n", length(ok), length(jobs))
    @printf("median |w_direct - w_spectral| = %.1e, max = %.1e\n", median(d), maximum(d))
    @printf("within d_max (or the double-precision floor 1e-9): %d of %d\n",
            count(t -> t[1] <= max(t[2], 1e-9), ok), length(ok))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? "results" : ARGS[1])
end
