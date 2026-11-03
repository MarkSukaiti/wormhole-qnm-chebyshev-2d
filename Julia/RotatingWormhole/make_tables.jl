# =====================================================================================
#  make_tables.jl
#
#  LaTeX tables of the converged quasinormal frequencies of the sectors m = 1, 2, 3
#  (Appendix A, Tables VIII-X): one row per converged root with Re varpi >= 0 and a defined
#  overtone number, labelled (L, n), quoted to the number of decimals (3-6) supported by
#  d_max, in three columns per table.
#
#  Usage:  julia --project=.. make_tables.jl  [results]
#  Input:  results/retained_roots.csv
#  Output: results/appendix_tables.tex
# =====================================================================================

using Printf

include("WormholeData.jl")
using .WormholeData

function fmt(w, dmax)
    dec = Int(min(6, max(3, floor(-log10(max(dmax, 1e-300))))))
    s = imag(w) < 0 ? "-" : "+"
    return @sprintf("\$%.*f%s%.*f\\,\\mathrm{i}\$", dec, real(w), s, dec, abs(imag(w)))
end

function table_for_m(R, m)
    rows = Tuple{Int,Bool,Any}[]
    for k in TAGS
        modes = [o for o in R if o.m == m && o.k == k && o.converged &&
                 real(o.varpi) >= -1e-9 && o.n >= 0]
        sort!(modes; by = o -> (o.L, o.n))
        for (i, o) in enumerate(modes)
            push!(rows, (k, i == 1, o))
        end
    end
    # three columns; a split point is moved to the start of the next j-block when that block
    # is short, otherwise the block continues in the next column with an empty j cell
    ncol = 3
    target = cld(length(rows), ncol)
    cols = Vector{Vector{Tuple{Int,Bool,Any}}}()
    start = 1
    for c in 1:ncol
        if c == ncol
            push!(cols, rows[start:end])
            break
        end
        cut = min(start + target, length(rows) + 1)
        nxt = cut
        while nxt <= length(rows) && !rows[nxt][2] && nxt - cut < 4
            nxt += 1
        end
        if nxt <= length(rows) && rows[nxt][2]
            cut = nxt
        end
        push!(cols, rows[start:cut-1])
        start = cut
    end
    io = IOBuffer()
    println(io, "\\begin{table}[H]")
    println(io, "\\caption{\\textit{Converged scalar quasinormal frequencies for \$m=$m\$, labelled by the multipole \$L\$ and the overtone number \$n\$ within each \$(L,m)\$ family. Only the member of each mirror pair with \$\\operatorname{Re}\\varpi\\geq0\$ is shown; its partner is \$3\\pi jm-\\widehat{\\omega}^{*}\$ (Corollary~\\ref{cor:mirror}). Entries are the three-resolution averages of Eq.~\\eqref{frequency_average_wormhole}, quoted to the number of decimals (\$\\leq6\$) supported by \$d_{\\max}\$. Roots that do not satisfy the criterion \\eqref{refined_criterion} are omitted, as are the few converged roots with \$\\operatorname{Im}\\widehat{\\omega}<-7.5\$ that are separated from their overtone sequence by unresolved members (Sec.~\\ref{subsec:checks}).}}")
    println(io, "\\label{tab:qnm-m$m}")
    println(io, "\\centering\n\\begingroup\n\\scriptsize")
    println(io, "\\setlength{\\tabcolsep}{3.0pt}\n\\setlength{\\arrayrulewidth}{0.4pt}\n\\renewcommand{\\arraystretch}{1.08}")
    for (ci, col) in enumerate(cols)
        println(io, "\\begin{minipage}[t]{0.325\\textwidth}\n\\centering")
        println(io, "\\begin{tabular}[t]{|c|c|c|}\n\\hline")
        println(io, "\$j\$ & \$(L,n)\$ & \$\\widehat{\\omega}\$ \\\\\n\\hline")
        prev = -1
        for (k, first, o) in col
            prev >= 0 && k != prev && println(io, "\\hline")
            jcell = k != prev ? "\$$(j_tex(k))\$" : ""
            println(io, "$jcell & \$($(o.L),$(o.n))\$ & $(fmt(o.w, o.dmax)) \\\\")
            prev = k
        end
        println(io, "\\hline\n\\end{tabular}")
        println(io, "\\end{minipage}", ci < length(cols) ? "\\hfill" : "")
    end
    println(io, "\\endgroup\n\\end{table}")
    return String(take!(io)), length(rows)
end

function main(resdir = "results")
    R = read_retained(joinpath(resdir, "retained_roots.csv"))
    parts = String[]
    for m in SECTORS
        t, n = table_for_m(R, m)
        push!(parts, t)
        println("m=$m: $n table rows")
    end
    write(joinpath(resdir, "appendix_tables.tex"), join(parts, "\n"))
    println("written ", joinpath(resdir, "appendix_tables.tex"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? "results" : ARGS[1])
end
