# =====================================================================================
#  make_figures.jl
#
#  Figures of Sec. IV (PDF + PNG):
#    spectra_overview   Fig. 1  converged spectra of m = 1, 2, 3 in the second-exterior frame
#    spectrum_anatomy   Fig. 2  raw spectrum, retained roots and zero-frequency artefacts
#    mj_scaling         Fig. 3  fundamental modes against mu = m j, direct-integration curves
#    lightring_limit    Fig. 4  co-rotating light ring and the approach to the large-mu limit
#
#  Usage:  julia --project=.. make_figures.jl  <QNM-output-root>  [results]  [plots]
#  Input:  raw spectra in <QNM-output-root>, results/retained_roots.csv,
#          results/continuation.csv (continue_mu.jl), results/limit_eigenvalues.csv
#          (limit_problem.jl)
# =====================================================================================

using Plots
using LaTeXStrings
using Printf

include("WormholeData.jl")
using .WormholeData

gr()
default(fontfamily = "Computer Modern", framestyle = :box, grid = false, tickdirection = :in,
        legendfontsize = 6, guidefontsize = 9, tickfontsize = 7, titlefontsize = 8,
        titlelocation = :left, markerstrokewidth = 0.4, linewidth = 0.8)

const INK2 = colorant"#52514e"
const CAT = [colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a", colorant"#eda100", colorant"#e87ba4"]
const RAMP = cgrad([colorant"#86b6ef", colorant"#2a78d6", colorant"#0d366b"])
const MARK = Dict(1 => :circle, 2 => :rect, 3 => :utriangle)
const UNRESOLVED = Dict((1, 10) => [1], (2, 10) => [2, 3, 4], (3, 9) => [3, 4, 5], (3, 10) => [3, 4, 5])

jcolor(k) = get(RAMP, k / 10)
ell(L) = sqrt(L * (L + 1))
frame_dragging(e) = pi / 2 - atan(e) - e / (1 + e^2)

"Monotone piecewise-cubic (Fritsch-Carlson) interpolation of y(x) at xs."
function pchip(x, y, xs)
    n = length(x)
    h = diff(x)
    d = diff(y) ./ h
    s = similar(y)
    s[1], s[n] = d[1], d[n-1]
    for i in 2:n-1
        s[i] = d[i-1] * d[i] <= 0 ? 0.0 :
               (w1 = 2h[i] + h[i-1]; w2 = h[i] + 2h[i-1]; (w1 + w2) / (w1 / d[i-1] + w2 / d[i]))
    end
    map(xs) do t
        i = clamp(searchsortedlast(x, t), 1, n - 1)
        u = (t - x[i]) / h[i]
        h00, h10, h01, h11 = 2u^3 - 3u^2 + 1, u^3 - 2u^2 + u, -2u^3 + 3u^2, u^3 - u^2
        h00 * y[i] + h10 * h[i] * s[i] + h01 * y[i+1] + h11 * h[i] * s[i+1]
    end
end

function read_continuation(fn)
    C = Dict{Tuple{Int,Int},Dict{Float64,ComplexF64}}()
    for line in readlines(fn)[2:end]
        c = split(line, ',')
        key = (parse(Int, c[1]), parse(Int, c[2]))
        get!(C, key, Dict{Float64,ComplexF64}())[parse(Float64, c[3])] =
            complex(parse(Float64, c[4]), parse(Float64, c[5]))
    end
    return C
end

function read_limit(fn)
    D = Dict{Tuple{Int,Int},ComplexF64}()
    for line in readlines(fn)[2:end]
        c = split(line, ',')
        D[(parse(Int, c[1]), parse(Int, c[2]))] = complex(parse(Float64, c[3]), parse(Float64, c[4]))
    end
    return D
end

nearest_mu(C, mu) = C[argmin(x -> abs(x - mu), collect(keys(C)))]

# ------------------------------------------------------------------ Fig. 1
function fig_overview(R, C, outdir)
    panels = []
    for m in SECTORS
        p = plot(xlims = (-0.6, 6.2), ylims = (-5.2, 0.15), xlabel = L"\mathrm{Re}\,\widehat{\omega}_-",
                 ylabel = m == 1 ? L"\mathrm{Im}\,\widehat{\omega}_-" : "", title = latexstring("m = $m"),
                 legend = false)
        vline!(p, [0.0]; color = INK2, ls = :dash, lw = 0.6)
        for L in (m, m + 1, m + 2)
            pts = ComplexF64[]
            for k in TAGS
                c = [o for o in R if o.m == m && o.k == k && o.converged && o.L == L && o.n == 0]
                isempty(c) || push!(pts, c[1].wminus)
            end
            plot!(p, real.(pts), imag.(pts); color = INK2, lw = 0.6)
            isempty(pts) || annotate!(p, real(pts[1]), imag(pts[1]) + 0.18, text(latexstring("$L"), 7))
        end
        for k in TAGS
            z = [o.wminus for o in R if o.m == m && o.k == k && o.converged && real(o.varpi) >= -1e-9]
            scatter!(p, real.(z), imag.(z); ms = 2.6, mc = jcolor(k), msc = :white)
        end
        for ((mm, kk), Ls) in UNRESOLVED
            mm == m || continue
            mu = mm * Float64(j_of_tag(kk))
            for L in Ls
                wm = nearest_mu(C[(L, 0)], mu)
                scatter!(p, [real(wm)], [imag(wm)]; ms = 2.8, mc = :white, msc = jcolor(kk), msw = 0.9)
            end
        end
        push!(panels, p)
    end
    # colour bar: one cell per sampled spin, coloured as the markers
    cb = heatmap([0.0], collect(TAGS), reshape(collect(TAGS) ./ 10, :, 1); color = RAMP, clims = (0, 1),
                 colorbar = false, xticks = false, yticks = (collect(TAGS), [latexstring(j_label(k)) for k in TAGS]),
                 ymirror = true, ylabel = L"j", framestyle = :box,
                 ylims = (-0.5, 10.5), xlims = (-0.5, 0.5), title = "")
    f = plot(panels..., cb; layout = grid(1, 4; widths = [0.315, 0.315, 0.315, 0.055]), size = (760, 290),
             left_margin = 3Plots.mm, bottom_margin = 4Plots.mm)
    savefig(f, joinpath(outdir, "spectra_overview.pdf")); savefig(f, joinpath(outdir, "spectra_overview.png"))
end

# ------------------------------------------------------------------ Fig. 2
function fig_anatomy(root, R, outdir)
    m, k = 1, 5
    j = Float64(j_of_tag(k))
    Om0 = 1.5pi * j
    pa = plot(xlims = (-9, 9), ylims = (-4.2, 0.45), xlabel = L"\mathrm{Re}\,\varpi",
              ylabel = L"\mathrm{Im}\,\varpi", title = L"(a) $m = 1$, $j = 1$: spectrum in the throat frame",
              legend = :bottom)
    for N in RESOLUTIONS
        raw = ComplexF64.(read_raw(raw_file(root, m, k, N))) .- m * Om0
        scatter!(pa, real.(raw), imag.(raw); ms = 1.0, mc = colorant"#c3c2b7", msw = 0,
                 label = N == 100 ? L"raw, $N_x = 90, 95, 100$" : "")
    end
    modes = [o for o in R if o.m == m && o.k == k]
    art = [o for o in modes if o.class == "zero-frequency artefact"]
    oth = [o for o in modes if o.class == "not converged"]
    cv = [o for o in modes if o.converged]
    scatter!(pa, real.(getfield.(oth, :varpi)), imag.(getfield.(oth, :varpi)); m = :circle, ms = 3,
             mc = :white, msc = INK2, label = "retained, not converged")
    scatter!(pa, real.(getfield.(art, :varpi)), imag.(getfield.(art, :varpi)); m = :xcross, ms = 3.5,
             mc = CAT[2], label = "zero-frequency artefacts")
    scatter!(pa, real.(getfield.(cv, :varpi)), imag.(getfield.(cv, :varpi)); m = :circle, ms = 3.6,
             mc = CAT[1], msc = :white, label = "converged QNMs")
    vline!(pa, [-m * Om0, m * Om0]; color = INK2, ls = :dash, lw = 0.7, label = "")
    vline!(pa, [0.0]; color = INK2, ls = :dot, lw = 0.6, label = "")
    # (b) artefacts at omega = 0: Im omega_N against N^-2
    pb = plot(xlims = (0, 1.32), ylims = (-5.5, 0.0), xlabel = L"10^4/N_x^2",
              ylabel = L"10^3\,\mathrm{Im}\,\widehat{\omega}_{N_x}", legend = false,
              title = L"(b) artefacts at $\widehat{\omega} = 0$: $m = 1$, $j = 1$")
    xs = [90, 95, 100] .^ -2.0 .* 1e4
    seen = Float64[]
    for o in modes
        (o.converged || abs(o.w) > 0.02) && continue
        ys = imag.([o.w90, o.w95, o.w100]) .* 1e3
        any(abs(ys[1] - y0) < 1e-3 for y0 in seen) && continue
        push!(seen, ys[1])
        c = sum(xs .* ys) / sum(xs .^ 2)
        plot!(pb, [0, 1.32], [0, 1.32c]; color = CAT[2], ls = :dash, lw = 0.6)
        scatter!(pb, xs, ys; ms = 3, mc = CAT[2], msc = :white)
    end
    f = plot(pa, pb; layout = grid(1, 2; widths = [0.6, 0.4]), size = (720, 300),
             left_margin = 3Plots.mm, bottom_margin = 4Plots.mm)
    savefig(f, joinpath(outdir, "spectrum_anatomy.pdf")); savefig(f, joinpath(outdir, "spectrum_anatomy.png"))
end

# ------------------------------------------------------------------ Fig. 3
function fig_scaling(R, C, outdir)
    panels = []
    for part in (:re, :im)
        p = plot(xscale = :log10, yscale = :log10, xlims = part == :re ? (0.09, 330) : (0.09, 160),
                 ylims = part == :re ? (0.04, 6.5) : (0.035, 0.62), xlabel = L"\mu = mj",
                 ylabel = part == :re ? L"\mathrm{Re}\,\widehat{\omega}_-" : L"-\mathrm{Im}\,\widehat{\omega}_-",
                 title = part == :re ? L"(a) oscillation frequency, $n = 0$" : L"(b) damping rate, $n = 0$",
                 grid = true, gridalpha = 0.25, legend = part == :im ? :bottomleft : false,
                 xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"]),
                 yticks = part == :re ? ([0.05, 0.1, 0.2, 0.5, 1, 2, 5], ["0.05", "0.1", "0.2", "0.5", "1", "2", "5"]) :
                          ([0.04, 0.06, 0.1, 0.15, 0.2, 0.3, 0.5], ["0.04", "0.06", "0.1", "0.15", "0.2", "0.3", "0.5"]))
        val(z) = part == :re ? real(z) : -imag(z)
        for L in 1:5
            Cl = C[(L, 0)]
            cm = sort([x for x in keys(Cl) if x > 0])
            cv = [val(Cl[x]) for x in cm]
            ts = range(log(cm[1]), log(cm[end]); length = 300)
            plot!(p, exp.(ts), exp.(pchip(log.(cm), log.(cv), ts)); color = CAT[L], lw = 0.9, label = "")
            for m in SECTORS
                pts = [o for o in R if o.converged && o.n == 0 && o.L == L && o.m == m && o.mu > 0]
                scatter!(p, getfield.(pts, :mu), [val(o.wminus) for o in pts]; m = MARK[m], ms = 3,
                         mc = CAT[L], msc = :white, label = "")
            end
            part == :re && annotate!(p, cm[end] * 1.15, cv[end], text(latexstring("L = $L"), 6, :left))
        end
        amp = part == :re ? 4.6 : 0.62
        plot!(p, [8.0, 100.0], amp .* ([8.0, 100.0] ./ 8) .^ -0.5; color = INK2, ls = :dash, lw = 0.8, label = "")
        if part == :im
            for m in SECTORS
                scatter!(p, [NaN], [NaN]; m = MARK[m], mc = :white, msc = :black, label = latexstring("spectral, \$m = $m\$"))
            end
            plot!(p, [NaN], [NaN]; color = INK2, label = "direct integration")
        end
        push!(panels, p)
    end
    f = plot(panels...; layout = (1, 2), size = (720, 290), left_margin = 3Plots.mm, bottom_margin = 4Plots.mm)
    savefig(f, joinpath(outdir, "mj_scaling.pdf")); savefig(f, joinpath(outdir, "mj_scaling.png"))
end

# ------------------------------------------------------------------ Fig. 4
function fig_lightring(R, C, LIM, outdir)
    L = 2
    l = ell(L)
    pa = plot(xlims = (-9.5, 3.5), ylims = (-1.3, 2.9), xlabel = L"\eta",
              ylabel = L"V_+(\eta) - m\widehat{\Omega}_-", title = L"(a) co-rotating light ring, $L = 2$",
              legend = :topleft)
    eta = range(-9.5, 3.5; length = 1200)
    cols = [get(RAMP, t) for t in (0.1, 0.4, 0.75, 1.0)]
    for (mu, col) in zip((0.0, 0.5, 2.0, 8.0), cols)
        V(e) = 3mu * frame_dragging(e) + l / sqrt(1 + e^2) - 3pi * mu
        plot!(pa, eta, V.(eta); color = col, lw = 1.1, label = latexstring("\\mu = $(mu == round(mu) ? Int(mu) : mu)"))
        cc = 6mu / l                                   # light ring: eta sqrt(1+eta^2) = -6 mu / l_L
        es = -sqrt((-1 + sqrt(1 + 4cc^2)) / 2)
        scatter!(pa, [es], [V(es)]; ms = 3.2, mc = col, msc = :white, label = "")
        wm = nearest_mu(C[(L, 0)], mu)
        scatter!(pa, [es], [real(wm)]; m = :hline, ms = 7, mc = col, msc = col, msw = 1.2, label = "")
    end
    hline!(pa, [0.0]; color = INK2, ls = :dash, lw = 0.6, label = "")
    annotate!(pa, 3.2, -0.22, text(L"window edge $\widehat{\omega}_- = 0$", 6, :right, INK2, "Computer Modern"))
    pb = plot(xlims = (-0.05, 3.85), ylims = (-2.25, 0.08), xlabel = L"\mu^{1/2}\,\mathrm{Re}\,\widehat{\omega}_-",
              ylabel = L"\mu^{1/2}\,\mathrm{Im}\,\widehat{\omega}_-", title = L"(b) approach to the limit $\mu \to \infty$",
              legend = false)
    for L in 1:5, (n, ls) in ((0, :solid), (1, :dash))
        haskey(C, (L, n)) || continue
        Cl = C[(L, n)]
        cm = sort([x for x in keys(Cl) if x >= 0.1])
        length(cm) < 4 && continue
        z = [sqrt(x) * Cl[x] for x in cm]
        ts = range(log(cm[1]), log(cm[end]); length = 400)
        plot!(pb, pchip(log.(cm), real.(z), ts), pchip(log.(cm), imag.(z), ts); color = CAT[L], ls = ls,
              lw = 0.9, label = "")
        for o in R
            if o.converged && o.L == L && o.n == n && o.mu >= 0.1
                zz = sqrt(o.mu) * o.wminus
                scatter!(pb, [real(zz)], [imag(zz)]; m = MARK[o.m], ms = 2.6, mc = CAT[L], msc = :white, label = "")
            end
        end
        if haskey(LIM, (L, n))
            scatter!(pb, [real(LIM[(L, n)])], [imag(LIM[(L, n)])]; m = :star5, ms = 5.5, mc = :black,
                     msc = :white, label = "")
        end
        if n == 0
            eik = 2 / (3sqrt(6)) * ell(L)^1.5 - 0.5im * sqrt(ell(L) / 3)     # Eq. (limit_eikonal)
            scatter!(pb, [real(eik)], [imag(eik)]; m = :diamond, ms = 3.2, mc = :white, msc = INK2, label = "")
            annotate!(pb, real(z[1]), imag(z[1]) + 0.08, text(latexstring("$L"), 7))
        end
    end
    # legend, drawn by hand in the empty lower-left corner
    x0, x1, xt, ys = 0.08, 0.30, 0.36, (-1.66, -1.80, -1.94, -2.08)
    plot!(pb, Shape([0.02, 1.18, 1.18, 0.02], [-1.57, -1.57, -2.17, -2.17]); fillcolor = :white,
          linecolor = :black, lw = 0.6)
    plot!(pb, [x0, x1], [ys[1], ys[1]]; color = INK2, lw = 0.9)
    plot!(pb, [x0, x1], [ys[2], ys[2]]; color = INK2, ls = :dash, lw = 0.9)
    scatter!(pb, [(x0 + x1) / 2], [ys[3]]; m = :star5, ms = 5.5, mc = :black, msc = :white)
    scatter!(pb, [(x0 + x1) / 2], [ys[4]]; m = :diamond, ms = 3.2, mc = :white, msc = INK2)
    for (y, s) in zip(ys, (L"$n = 0$", L"$n = 1$", L"limit $w_{L,n}$", L"light ring, $n = 0$"))
        annotate!(pb, xt, y, text(s, 6, :left, "Computer Modern"))
    end
    f = plot(pa, pb; layout = (1, 2), size = (720, 300), left_margin = 3Plots.mm, bottom_margin = 4Plots.mm)
    savefig(f, joinpath(outdir, "lightring_limit.pdf")); savefig(f, joinpath(outdir, "lightring_limit.png"))
end

function main(root, resdir = "results", outdir = "plots")
    setprecision(BigFloat, 760)
    mkpath(outdir)
    R = read_retained(joinpath(resdir, "retained_roots.csv"))
    C = read_continuation(joinpath(resdir, "continuation.csv"))
    LIM = read_limit(joinpath(resdir, "limit_eigenvalues.csv"))
    fig_overview(R, C, outdir)
    fig_anatomy(root, R, outdir)
    fig_scaling(R, C, outdir)
    fig_lightring(R, C, LIM, outdir)
    println("figures written to ", outdir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("usage: julia make_figures.jl <QNM-output-root> [results] [plots]")
    main(ARGS[1], length(ARGS) >= 2 ? ARGS[2] : "results", length(ARGS) >= 3 ? ARGS[3] : "plots")
end
