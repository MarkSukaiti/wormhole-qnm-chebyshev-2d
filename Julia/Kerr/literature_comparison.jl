# =====================================================================================
#  literature_comparison.jl
#
#  Comparison of the converged two-dimensional spectral Kerr frequencies with published
#  tabulations (Sec. V B, Table VII; Appendix B, Tables XI and XII).
#
#  reference/published_values.csv lists every published massless-scalar value used, with
#  the digits as quoted: Konoplya & Zhidenko (2006), Berti & Cardoso (2006), Lin et al.
#  (2017; continued fraction and matrix method, units 2M = 1), Leaver (1985) and the third-
#  and sixth-order WKB values as tabulated by Konoplya (2003), Mamani et al. (2022) and
#  Cho et al. (2012).  For each value the nearest converged spectral root is taken (for
#  m = -1 the mirror image -conj(omega) in the computed m = 1 spectrum), and
#
#   * delta = |w_SM - w_ref| / |w_SM|;
#   * whether the spectral value, rounded in the units and to the decimals of the source,
#     reproduces every quoted digit, and otherwise by how many units of the last digit it
#     differs and how close the accurate (qnm) value lies to a rounding boundary.
#
#  Usage:  julia --project=.. literature_comparison.jl  <Kerr-output-folder>  [results]
#  Output: results/kerr_published_comparison.csv, results/literature_tables.tex
# =====================================================================================

using Printf

include("KerrData.jl")
using .KerrData

struct Entry
    source::String
    method::String
    l::Int
    m::Int
    n::Int
    a::Float64
    re::String
    im::String
    units::String
    dre::Int
    dim::Int
end

function read_published(fn)
    E = Entry[]
    for line in readlines(fn)[2:end]
        c = split(line, ',')
        push!(E, Entry(c[1], c[2], parse(Int, c[3]), parse(Int, c[4]), parse(Int, c[5]),
                       parse(Float64, c[6]), c[7], c[8], c[9], parse(Int, c[10]), parse(Int, c[11])))
    end
    return E
end

scale(e::Entry) = e.units == "2M" ? 0.5 : 1.0
published(e::Entry) = complex(parse(Float64, e.re), -parse(Float64, e.im)) * scale(e)   # M*omega

"Deviation of x from the quoted value in units of the last quoted digit."
lastdigit(x, quoted::AbstractString, d) = abs(x - parse(Float64, quoted)) * 10.0^d

"Distance of x from the nearest rounding boundary, in units of the d-th decimal."
boundary_distance(x, d) = (u = x * 10.0^d; abs(u - floor(u) - 0.5))

function main(root, resdir = "results")
    mkpath(resdir)
    ref = read_reference(joinpath(@__DIR__, "reference", "kerr_qnm_reference.csv"))
    E = read_published(joinpath(@__DIR__, "reference", "published_values.csv"))
    conv = Dict(a_of_tag(k) => [r.av for r in read_diagnostics(diagnostics_file(root, k))] for k in TAGS)
    rows = []
    for e in E
        p = published(e)
        target = e.m == 1 ? p : -conj(p)
        roots = conv[e.a]
        w = roots[argmin(abs.(roots .- target))]
        sm = e.m == 1 ? w : -conj(w)                  # spectral value for the tabulated mode
        q = ref[(e.l, e.m, e.n, e.a)]
        s = scale(e)
        ure = lastdigit(real(sm) / s, e.re, e.dre)
        uim = lastdigit(-imag(sm) / s, e.im, e.dim)
        all_digits = max(ure, uim) < 0.5
        bd = min(ure >= 0.5 ? boundary_distance(real(q) / s, e.dre) : Inf,
                 uim >= 0.5 ? boundary_distance(-imag(q) / s, e.dim) : Inf)
        push!(rows, (e = e, p = p, sm = sm, q = q, delta = abs(sm - p) / abs(sm),
                     delta_qnm = abs(sm - q) / abs(q), units = max(ure, uim), all = all_digits, bd = bd))
    end
    open(joinpath(resdir, "kerr_published_comparison.csv"), "w") do io
        println(io, "source,method,l,m,n,a_over_M,Re_Momega_published,minusIm_Momega_published,Re_Momega_SM,minusIm_Momega_SM,Re_Momega_qnm,minusIm_Momega_qnm,rel_dev_SM_vs_published,rel_dev_SM_vs_qnm,all_digits")
        for r in rows
            e = r.e
            @printf(io, "%s,%s,%d,%d,%d,%g,%.7f,%.7f,%.12f,%.12f,%.12f,%.12f,%.3e,%.3e,%s\n",
                    e.source, e.method, e.l, e.m, e.n, e.a, real(r.p), -imag(r.p), real(r.sm),
                    -imag(r.sm), real(r.q), -imag(r.q), r.delta, r.delta_qnm, r.all)
        end
    end
    println("=== summary per source (Table VII)")
    groups = unique([(r.e.source, r.e.method) for r in rows])
    for g in groups
        G = [r for r in rows if (r.e.source, r.e.method) == g]
        @printf("  %-52s %-15s N=%2d  all digits=%2d  max dev %5.2f last-digit units  max delta %.1e\n",
                g[1], g[2], length(G), count(r -> r.all, G), maximum(r -> r.units, G),
                maximum(r -> r.delta, G))
    end
    exact = [r for r in rows if !occursin("WKB", r.e.method) && r.e.method != "matrix method"]
    miss = [r for r in exact if !r.all]
    @printf("continued-fraction, pseudospectral and AIM values: %d, all digits reproduced: %d; the %d others differ by one unit in the last digit, the accurate value lying within %.3f units of a rounding boundary\n",
            length(exact), count(r -> r.all, exact), length(miss), maximum(r -> r.bd, miss; init = 0.0))
    @printf("largest deviation of a listed spectral value from qnm: %.1e\n", maximum(r -> r.delta_qnm, rows))
    write_tables(rows, joinpath(resdir, "literature_tables.tex"))
end

cpx(w, d) = @sprintf("\$%.*f-%.*f\\,\\mathrm{i}\$", d, real(w), d, -imag(w))
"Half of a quoted decimal number, exactly (one more decimal, trailing zeros removed)."
halve(s::AbstractString, d) = (t = @sprintf("%.*f", d + 1, parse(Float64, s) / 2); occursin('.', t) ? rstrip(rstrip(t, '0'), '.') : t)
pubstr(e::Entry) = e.units == "2M" ?
    "\$$(halve(e.re, e.dre))-$(halve(e.im, e.dim))\\,\\mathrm{i}\$" :
    "\$$(e.re)-$(e.im)\\,\\mathrm{i}\$"

"LaTeX rows of Tables XI (spinning) and XII (Schwarzschild limit)."
function write_tables(rows, fn)
    io = IOBuffer()
    println(io, "% Table XI rows: a/M & SM & CF (KZ or BC) & CF (Lin) & matrix method (Lin)")
    blocks = [((1, 1, 0), "KonoplyaZhidenko2006"), ((2, 1, 0), "BertiCardoso2006"),
              ((2, 1, 1), "BertiCardoso2006"), ((2, -1, 0), "BertiCardoso2006"),
              ((2, -1, 1), "BertiCardoso2006")]
    for ((l, m, n), src) in blocks
        println(io, "\\multicolumn{5}{l}{\$(\\ell,m,n)=($l,$m,$n)\$; CF values of Ref.~\\cite{$src}} \\\\")
        B = sort([r for r in rows if r.e.source == src && (r.e.l, r.e.m, r.e.n) == (l, m, n)]; by = r -> r.e.a)
        for r in B
            lin = [x for x in rows if startswith(x.e.source, "LinEtAl") && (x.e.l, x.e.m, x.e.n, x.e.a) == (l, m, n, r.e.a)]
            cf = [pubstr(x.e) for x in lin if x.e.method == "CF"]
            mm = [pubstr(x.e) for x in lin if x.e.method == "matrix method"]
            @printf(io, "\$%g\$ & %s & %s & %s & %s \\\\\n", r.e.a, cpx(r.sm, 8), pubstr(r.e),
                    isempty(cf) ? "" : cf[1], isempty(mm) ? "" : mm[1])
        end
    end
    println(io, "% Table XII rows: (l,n) & SM & CF & PS/AIM & AIM & WKB6 & WKB3")
    for (l, n) in sort(unique([(r.e.l, r.e.n) for r in rows if r.e.a == 0 && r.e.source == "Mamani2022EPJC"]))
        cell(src) = (c = [r for r in rows if startswith(r.e.source, src) && (r.e.l, r.e.n, r.e.a) == (l, n, 0.0)];
                     isempty(c) ? "---" : pubstr(c[1].e))
        sm = first(r.sm for r in rows if r.e.source == "Mamani2022EPJC" && (r.e.l, r.e.n) == (l, n))
        @printf(io, "\$(%d,%d)\$ & %s & %s & %s & %s & %s & %s \\\\\n", l, n, cpx(sm, 8),
                cell("refLeaver"), cell("Mamani"), cell("ChoEtAl"), cell("Konoplya2003PRD"), cell("IyerWill"))
    end
    write(fn, String(take!(io)))
    println("written ", fn)
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia literature_comparison.jl <Kerr-output-folder> [results]")
    main(ARGS[1], length(ARGS) >= 2 ? ARGS[2] : "results")
end
