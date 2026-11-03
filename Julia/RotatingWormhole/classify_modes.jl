# =====================================================================================
#  classify_modes.jl
#
#  Post-processing of the three-resolution output of QNMThreads.jl (Secs. III B, IV A):
#
#   1. checks that the union of the Arb radial spectra (radial_pencils_arb.jl) reproduces
#      every eigenvalue of the raw N_x = 90 two-dimensional spectrum (Proposition III.1);
#   2. labels each retained root by its multipole L (nearest radial eigenvalue);
#   3. applies the refined, mirror-symmetric convergence criterion, Eq. (refined_criterion):
#          d_max <= eps_q * min(|w|, |w - 3 pi m j|),   eps_q = 1e-3,
#      and discards roots on the lines Re w = 0 and Re w_- = 0;
#   4. classifies the rejected roots (zero-frequency artefacts: distance to the nearest
#      zero-frequency point shrinks like N_x^-2, i.e. by a factor 1.10-1.40 from 90 to 100);
#   5. assigns overtone numbers n within each (L, m) family, ordered by damping, stopping
#      at a gap that signals unresolved intermediate overtones;
#   6. writes CSV files and prints the counts quoted in Sec. IV A.
#
#  Usage:
#      julia --project=.. classify_modes.jl  <QNM-output-root>  [radial_spectra]  [results]
#  where <QNM-output-root> contains output-m1, output-m2, output-m3 written by QNMThreads.jl.
#  Output: results/converged_qnms.csv, results/rejected_roots.csv, results/retained_roots.csv
# =====================================================================================

using Printf
using Statistics

include("WormholeData.jl")
using .WormholeData

set_working_precision()

const EPS_Q = 1e-3          # refined criterion
const LINE_TOL = 1e-6       # roots on Re w = 0 or Re w_- = 0 are discarded
const LABEL_TOL = 1e-20     # relative distance to a radial eigenvalue for a multipole label

function read_radial(dir, m, k, L)
    fn = joinpath(dir, "radial-m$m-J$k-L$L-N90.txt")
    vals = Complex{BigFloat}[]
    for line in eachline(fn)
        startswith(line, "#") && continue
        a, b = split(line)
        push!(vals, Complex(parse(BigFloat, a), parse(BigFloat, b)))
    end
    return vals
end

"Index and BigFloat distance of the element of `pool` nearest to w (Float64 prefilter)."
function nearest(pool::Vector{Complex{BigFloat}}, poolf::Vector{ComplexF64}, w)
    wf = ComplexF64(w)
    idx = partialsortperm(abs.(poolf .- wf), 1:min(4, length(poolf)))
    best = argmin(i -> abs(pool[i] - w), idx)
    return best, abs(pool[best] - w)
end

function process(root, raddir)
    modes = NamedTuple[]
    repro = []
    for m in SECTORS, k in TAGS
        jr = j_of_tag(k)
        j = BigFloat(jr)
        e2 = 3 * BigFloat(π) * m * j                       # m * Omega_minus (hat)
        # union of the three radial spectra, with labels
        pool = Complex{BigFloat}[]
        labs = Int[]
        for L in (m, m + 1, m + 2)
            r = read_radial(raddir, m, k, L)
            append!(pool, r)
            append!(labs, fill(L, length(r)))
        end
        poolf = ComplexF64.(pool)
        # (1) reproduction of the full raw N = 90 spectrum
        raw = read_raw(raw_file(root, m, k, 90))
        d = [Float64(nearest(pool, poolf, w)[2] / max(1, abs(w))) for w in raw]
        push!(repro, (m, k, median(d), maximum(d), count(<(1e-25), d), length(d)))
        # (2)-(4) labels and classification of the retained roots
        fam = Dict{Int,Vector{Int}}()
        for r in read_diagnostics(diagnostics_file(root, m, k))
            i, dl = nearest(pool, poolf, r.w1)
            L = dl < LABEL_TOL * max(1, abs(r.w1)) ? labs[i] : 0
            w = r.av
            dedge = min(abs(w), abs(w - e2))
            q = r.dmax / max(dedge, BigFloat("1e-300"))
            online = abs(real(w)) < LINE_TOL || abs(real(w) - e2) < LINE_TOL
            conv = q <= EPS_Q && !online
            e = abs(w) <= abs(w - e2) ? zero(e2) : e2
            ratio = abs(r.w1 - e) / max(abs(r.w3 - e), BigFloat("1e-300"))
            cls = conv ? "converged" :
                  (1.10 < ratio < 1.40 && q > 0.05) ? "zero-frequency artefact" : "not converged"
            push!(modes, (m = m, k = k, j = j, mu = m * j, w = w, w90 = r.w1, w95 = r.w2, w100 = r.w3,
                          wminus = w - e2, varpi = w - e2 / 2, dmax = r.dmax, dedge = dedge,
                          q = q, L = L, class = cls, n = Ref{Int}(-1)))
        end
    end
    # (5) overtone numbers among converged members with Re varpi >= 0, per (m, k, L)
    for m in SECTORS, k in TAGS, L in (m, m + 1, m + 2)
        fam = [o for o in modes if o.m == m && o.k == k && o.L == L &&
               o.class == "converged" && real(o.varpi) >= -1e-9]
        sort!(fam; by = o -> -imag(o.w))
        prev_gap = 0.0
        for (n, o) in enumerate(fam)
            if n > 1
                gap = Float64(imag(fam[n-1].w) - imag(o.w))
                gap > max(2.5, 1.6 * prev_gap) && break        # unresolved overtones: n undefined
                prev_gap = gap
            end
            o.n[] = n - 1
        end
    end
    return modes, repro
end

nstr(x, d = 25) = @sprintf("%.*e", d - 1, x)

function write_csvs(modes, outdir)
    mkpath(outdir)
    open(joinpath(outdir, "converged_qnms.csv"), "w") do io
        println(io, "m,k,jstar,j,mu=m*j,L,n,Re_omega,Im_omega,Re_omega_minus,Im_omega_minus,Re_varpi,d_max,q=d_max/dist_to_zero_freq")
        for o in sort([o for o in modes if o.class == "converged"]; by = o -> (o.m, o.k, o.L, -imag(o.w)))
            n = o.n[] < 0 ? "" : string(o.n[])
            @printf(io, "%d,%d,%s,%.15g,%.15g,%d,%s,%s,%s,%.12f,%.12f,%.12f,%.3e,%.3e\n",
                    o.m, o.k, string(Float64(jstar_of_tag(o.k))), o.j, o.mu, o.L, n,
                    nstr(real(o.w)), nstr(imag(o.w)), real(o.wminus), imag(o.wminus),
                    real(o.varpi), o.dmax, o.q)
        end
    end
    open(joinpath(outdir, "rejected_roots.csv"), "w") do io
        println(io, "m,k,j,Re_omega,Im_omega,Re_omega_N90,Im_omega_N90,Re_omega_N100,Im_omega_N100,d_max,dist_to_zero_freq,q,class")
        for o in modes
            o.class == "converged" && continue
            @printf(io, "%d,%d,%.15g,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.3e,%.3e,%.3e,%s\n",
                    o.m, o.k, o.j, real(o.w), imag(o.w), real(o.w90), imag(o.w90),
                    real(o.w100), imag(o.w100), o.dmax, o.dedge, o.q, o.class)
        end
    end
    # everything, for the downstream scripts (tables, figures, direct validation)
    open(joinpath(outdir, "retained_roots.csv"), "w") do io
        println(io, "m,k,j,mu,L,n,class,Re_omega,Im_omega,Re_omega_N90,Im_omega_N90,Re_omega_N95,Im_omega_N95,Re_omega_N100,Im_omega_N100,d_max,q")
        for o in modes
            n = o.n[] < 0 ? "" : string(o.n[])
            @printf(io, "%d,%d,%.17g,%.17g,%d,%s,%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.6e,%.6e\n",
                    o.m, o.k, o.j, o.mu, o.L, n, o.class, nstr(real(o.w)), nstr(imag(o.w)),
                    real(o.w90), imag(o.w90), real(o.w95), imag(o.w95), real(o.w100), imag(o.w100),
                    o.dmax, o.q)
        end
    end
end

function summary(modes, repro)
    println("=== reproduction of the raw N_x = 90 spectra by the Arb radial pencils (relative distance)")
    for (m, k, med, mx, ngood, ntot) in repro
        @printf("  m=%d k=%2d: median %.1e  max %.1e  (< 1e-25: %d/%d)\n", m, k, med, mx, ngood, ntot)
    end
    tot = length(modes)
    conv = count(o -> o.class == "converged", modes)
    art = count(o -> o.class == "zero-frequency artefact", modes)
    @printf("retained by the first-stage criterion: %d; converged: %d (%d pairs); zero-frequency artefacts: %d; not converged: %d\n",
            tot, conv, conv ÷ 2, art, tot - conv - art)
    println("converged with Im w >= 0: ", count(o -> o.class == "converged" && imag(o.w) >= 0, modes))
    for o in modes
        imag(o.w) >= 0 && @printf("  retained root in the upper half-plane: m=%d k=%d w=%s q=%.2e (%s)\n",
                                   o.m, o.k, string(ComplexF64(o.w)), o.q, o.class)
    end
    cv = [o for o in modes if o.class == "converged"]
    best = argmax(o -> imag(o.w), [o for o in cv if real(o.varpi) >= -1e-9])
    @printf("least damped converged root (Re varpi >= 0 member): m=%d k=%d L=%d n=%d  w = %s\n", best.m, best.k, best.L,
            best.n[], string(ComplexF64(best.w)))
    println("unlabelled converged roots: ", count(o -> o.L == 0, cv))
    println("converged roots with Re varpi >= 0 left without overtone number (gaps): ",
            count(o -> real(o.varpi) >= -1e-9 && o.n[] < 0, cv))
    @printf("converged q: median %.1e\n", median(Float64.(getfield.(cv, :q))))
    println("=== converged mirror pairs per (m, j)   [Table II]")
    for k in TAGS
        @printf("  j=%-4s", j_label(k))
        for m in SECTORS
            c = [o for o in cv if o.m == m && o.k == k]
            mx = isempty(c) ? "   ---  " : @sprintf("%8.4f", maximum(o -> Float64(imag(o.w)), c))
            @printf("   m=%d: N_c=%3d  max Im=%s", m, length(c) ÷ 2, mx)
        end
        println()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia classify_modes.jl <QNM-output-root> [radial_spectra] [results]")
    root = ARGS[1]
    raddir = length(ARGS) >= 2 ? ARGS[2] : "radial_spectra"
    outdir = length(ARGS) >= 3 ? ARGS[3] : "results"
    modes, repro = process(root, raddir)
    write_csvs(modes, outdir)
    summary(modes, repro)
end
