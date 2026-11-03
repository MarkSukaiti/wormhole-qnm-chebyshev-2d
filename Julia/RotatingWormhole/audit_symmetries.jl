# =====================================================================================
#  audit_symmetries.jl
#
#  Tests of the exact properties of the discrete pencil on the RAW spectra written by
#  QNMThreads.jl, in the working precision of the eigenvalue solver (Sec. IV A):
#
#   * mirror symmetry (Corollary II.3):  w  ->  3 pi m j - conj(w)  maps every raw spectrum
#     onto itself, at every resolution;
#   * dependence on (m, j) only through mu = m j (Proposition II.4): runs with equal mu share
#     their common multipole families, exactly 2 N_x eigenvalues per family, while runs with
#     different mu share none;
#   * the static runs share their common families in the same way.
#
#  These agreements test exact identities of the discrete pencil; they say nothing about the
#  accuracy of the eigenvalues as approximations of the continuum frequencies, which is
#  estimated by d_max (classify_modes.jl).
#
#  Usage:  julia --project=.. audit_symmetries.jl  <QNM-output-root>
# =====================================================================================

using Printf

include("WormholeData.jl")
using .WormholeData

set_working_precision()

function closest(B::Vector{Complex{BigFloat}}, Bf::Vector{ComplexF64}, z)
    idx = partialsortperm(abs.(Bf .- ComplexF64(z)), 1:min(3, length(Bf)))
    return minimum(abs(B[i] - z) for i in idx)
end

function mirror_test(raw, m, k)
    e2 = 3 * BigFloat(π) * m * BigFloat(j_of_tag(k))
    f = ComplexF64.(raw)
    worst = zero(BigFloat)
    worst_rel = zero(BigFloat)
    for w in raw
        d = closest(raw, f, e2 - conj(w))
        worst = max(worst, d)
        worst_rel = max(worst_rel, d / max(1, abs(w)))
    end
    return worst, worst_rel
end

function shared_count(A, B, thr)
    Bf = ComplexF64.(B)
    n = 0
    dmax = zero(BigFloat)
    for w in A
        d = closest(B, Bf, w)
        if d < thr * max(1, abs(w))
            n += 1
            dmax = max(dmax, d)
        end
    end
    return n, dmax
end

function main(root)
    raw(m, k, N) = read_raw(raw_file(root, m, k, N))
    println("=== mirror symmetry of the raw spectra (all eigenvalues, all runs)")
    overall = zero(BigFloat)
    for m in SECTORS, k in TAGS
        line = String[]
        for N in RESOLUTIONS
            w, wr = mirror_test(raw(m, k, N), m, k)
            overall = max(overall, w)
            push!(line, @sprintf("N=%d: %.1e", N, Float64(w)))
        end
        @printf("  m=%d k=%2d  %s\n", m, k, join(line, "   "))
    end
    @printf("  largest distance over all %d raw spectra: %.2e\n",
            length(SECTORS) * length(TAGS) * length(RESOLUTIONS), Float64(overall))

    println("=== shared multipole families between independent runs (threshold 1e-120 relative)")
    thr = BigFloat("1e-120")
    pairs = [((1, 0), (2, 0), 2, "j=0: L=2,3 shared"),
             ((2, 0), (3, 0), 2, "j=0: L=3,4 shared"),
             ((1, 0), (3, 0), 1, "j=0: L=3 shared"),
             ((2, 5), (3, 4), 2, "mu=2: (m=2,j=1) vs (m=3,j=2/3): L=3,4 shared"),
             ((2, 6), (3, 5), 2, "mu=3: (m=2,j=3/2) vs (m=3,j=1): L=3,4 shared"),
             ((1, 5), (2, 5), 0, "control: (m=1,j=1) vs (m=2,j=1): different mu")]
    for (a, b, nfam, label) in pairs, N in (90, 100)
        n, dm = shared_count(raw(a..., N), raw(b..., N), thr)
        @printf("  %-50s N=%d: shared %4d (expected %4d), max difference %.1e\n",
                label, N, n, 2N * nfam, Float64(dm))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia audit_symmetries.jl <QNM-output-root>")
    main(ARGS[1])
end
