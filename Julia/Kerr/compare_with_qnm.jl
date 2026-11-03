# =====================================================================================
#  compare_with_qnm.jl
#
#  Comparison of the two-dimensional spectral Kerr frequencies (QNMThreads_Nu.jl) with the
#  continued-fraction values of the qnm package (Sec. V, Tables IV and V):
#
#   * every converged root is identified with the nearest reference mode (l <= 9, m = +-1,
#     the m = -1 modes as mirror images -conj(omega)); roots with no reference within a
#     relative distance 1e-5 are reported separately (near-extremal branch-cut roots);
#   * Table IV: fundamental modes l = 1, 2 at all spins, with delta = |w_SM - w_qnm|/|w_qnm|;
#   * Table V: log10(delta) for l = 1, 2 and n <= 3, '---' where no converged root matches;
#     for those modes the deviation of the nearest raw root at each resolution is printed;
#   * delta is compared with the convergence diagnostic d_max/|omega|.
#
#  Usage:  julia --project=.. compare_with_qnm.jl  <Kerr-output-folder>  [results]
#  where <Kerr-output-folder> contains J0 ... J10 written by QNMThreads_Nu.jl.
#  Output: results/kerr_converged_identified.csv, results/kerr_tables.tex
# =====================================================================================

using Printf
using Statistics

include("KerrData.jl")
using .KerrData

const MATCH_TOL = 1e-5

function identify(w, a, ref)
    best = (Inf, (0, 0, 0))
    for (l, m, n, aa) in keys(ref)
        aa == a || continue
        r = reference_value(ref, l, m, n, a)
        d = abs(w - r) / abs(r)
        d < best[1] && (best = (d, (l, m, n)))
    end
    return best
end

function main(root, resdir = "results")
    mkpath(resdir)
    ref = read_reference(joinpath(@__DIR__, "reference", "kerr_qnm_reference.csv"))
    allrows = []
    for k in TAGS
        a = a_of_tag(k)
        for r in read_diagnostics(diagnostics_file(root, k))
            d, (l, m, n) = identify(r.av, a, ref)
            push!(allrows, (k = k, a = a, w = r.av, dmax = r.dmax, delta = d, l = l, m = m, n = n,
                            matched = d < MATCH_TOL))
        end
    end
    open(joinpath(resdir, "kerr_converged_identified.csv"), "w") do io
        println(io, "a_over_M,Re_Momega,Im_Momega,dmax_over_absomega,l,m,n,delta,matched")
        for r in allrows
            @printf(io, "%.2f,%.15f,%.15f,%.2e,%d,%d,%d,%.2e,%s\n", r.a, real(r.w), imag(r.w),
                    r.dmax / abs(r.w), r.l, r.m, r.n, r.delta, r.matched)
        end
    end
    matched = [r for r in allrows if r.matched]
    @printf("converged roots: %d; identified with qnm modes: %d (max delta %.1e, l <= %d, n <= %d)\n",
            length(allrows), length(matched), maximum(getfield.(matched, :delta)),
            maximum(getfield.(matched, :l)), maximum(getfield.(matched, :n)))
    for r in allrows
        r.matched || @printf("  unidentified: a/M = %.2f  M*omega = %s\n", r.a, string(r.w))
    end

    # Tables IV and V: l = 1, 2, m = 1
    best = Dict{Tuple{Int,Int,Float64},Any}()
    for r in matched
        r.m == 1 && r.l <= 2 && r.n <= 3 || continue
        key = (r.l, r.n, r.a)
        (!haskey(best, key) || r.delta < best[key].delta) && (best[key] = r)
    end
    tex = IOBuffer()
    println("=== Table IV: fundamental modes (M*omega) and delta")
    for k in TAGS
        a = a_of_tag(k)
        cells = String[]
        for l in (1, 2)
            r = best[(l, 0, a)]
            push!(cells, @sprintf("\$%.10f%+.10f\\,\\mathrm{i}\$ & \$%s\$", real(r.w), imag(r.w), sci(r.delta)))
        end
        line = @sprintf("\$%g\$ & %s \\\\", a, join(cells, " & "))
        println(tex, line); println("  ", line)
    end
    println("=== Table V: log10(delta), l = 1, 2, n = 0..3")
    for k in TAGS
        a = a_of_tag(k)
        cells = [haskey(best, (l, n, a)) ? @sprintf("\$%.1f\$", log10(best[(l, n, a)].delta)) : "---"
                 for l in (1, 2) for n in 0:3]
        line = @sprintf("\$%g\$ & %s \\\\", a, join(cells, " & "))
        println(tex, line); println("  ", line)
    end
    write(joinpath(resdir, "kerr_tables.tex"), String(take!(tex)))

    # modes that do not pass the three-resolution test: nearest raw root at each resolution
    println("=== l = 1, 2, n <= 3 modes without a converged root: deviation of the nearest raw root")
    for k in TAGS, l in (1, 2), n in 0:3
        a = a_of_tag(k)
        haskey(best, (l, n, a)) && continue
        w0 = reference_value(ref, l, 1, n, a)
        devs = [minimum(abs.(read_raw(raw_file(root, k, res)) .- w0)) / abs(w0) for res in RESOLUTIONS]
        @printf("  a/M = %.2f (l,n) = (%d,%d): %s\n", a, l, n, join([@sprintf("%.1e", d) for d in devs], " -> "))
    end
    sel = [best[key] for key in keys(best)]
    @printf("l = 1, 2, n <= 3: %d of 88 modes converged; max delta %.1e; median delta n=0: %.0e, n=1: %.0e\n",
            length(sel), maximum(getfield.(sel, :delta)),
            median([r.delta for r in sel if r.n == 0]), median([r.delta for r in sel if r.n == 1]))
    viol = [r for r in sel if r.delta > r.dmax / abs(r.w)]
    println("modes with delta > d_max/|omega|: ", length(viol),
            join([@sprintf("  (%d,%d,a=%.2f: delta %.1e, d_max/|w| %.1e)", r.l, r.n, r.a, r.delta, r.dmax / abs(r.w)) for r in viol]))
end

sci(x) = (e = floor(Int, log10(x)); @sprintf("%.1f\\times10^{%d}", x / 10.0^e, e))

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia compare_with_qnm.jl <Kerr-output-folder> [results]")
    main(ARGS[1], length(ARGS) >= 2 ? ARGS[2] : "results")
end
