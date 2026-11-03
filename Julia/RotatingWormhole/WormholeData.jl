"""
    WormholeData

Shared definitions for the post-processing scripts of the rotating-wormhole spectra:
spin tags, file locations and readers for the output of `QNMThreads.jl`.

Conventions (as in the Maple matrix generator and the paper):

* azimuthal sectors m = 1, 2, 3, run separately (folders output-m1, output-m2, output-m3);
* spin tag k = 0..10 encodes j* = k/10 (k = 10: j* = 0.95) and the physical spin
  j = J/r0^2 = j*/(1 - j*), i.e. j in {0, 1/9, 1/4, 3/7, 2/3, 1, 3/2, 7/3, 4, 9, 19};
* frequencies are w = r0*omega (first-exterior frame); w_- = w - 3 pi m j is the frequency
  in the frame of the second exterior and varpi = w - (3/2) pi m j the throat-frame value;
* three radial resolutions N_x = 90, 95, 100 and N_chi = 3.
"""
module WormholeData

using Printf

export SECTORS, TAGS, RESOLUTIONS, jstar_of_tag, j_of_tag, j_label, j_tex,
       read_diagnostics, read_raw, output_dir, diagnostics_file, positive_file, raw_file,
       set_working_precision

const SECTORS = (1, 2, 3)
const TAGS = 0:10
const RESOLUTIONS = (90, 95, 100)
const N_CHI = 3

"Working precision used for reading the 220-digit output files."
set_working_precision() = setprecision(BigFloat, 760)

jstar_of_tag(k::Integer) = k == 10 ? 95 // 100 : k // 10
j_of_tag(k::Integer) = (s = jstar_of_tag(k); s / (1 - s))           # exact Rational

const J_LABEL = Dict(0 => "0", 1 => "1/9", 2 => "1/4", 3 => "3/7", 4 => "2/3", 5 => "1",
                     6 => "3/2", 7 => "7/3", 8 => "4", 9 => "9", 10 => "19")
const J_TEX = Dict(0 => "0", 1 => "\\tfrac19", 2 => "\\tfrac14", 3 => "\\tfrac37",
                   4 => "\\tfrac23", 5 => "1", 6 => "\\tfrac32", 7 => "\\tfrac73", 8 => "4",
                   9 => "9", 10 => "19")
j_label(k) = J_LABEL[k]
j_tex(k) = J_TEX[k]

output_dir(root, m, k) = joinpath(root, "output-m$m", "J$k")
diagnostics_file(root, m, k) = joinpath(output_dir(root, m, k), "stable-modes-diagnostics-J$k.txt")
positive_file(root, m, k) = joinpath(output_dir(root, m, k), "stable-positive-imaginary-diagnostics-J$k.txt")
raw_file(root, m, k, N) = joinpath(output_dir(root, m, k), "raw-modes-N$N-M$N_CHI-J$k.txt")

"""
    read_diagnostics(path) -> Vector{NamedTuple}

Rows of a three-resolution diagnostics file: the average `av`, the three matched roots
`w1, w2, w3` (N_x = 90, 95, 100), the pairwise distances, `dmax` and the tolerance.
"""
function read_diagnostics(path::AbstractString)
    rows = NamedTuple[]
    isfile(path) || return rows
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        c = split(s)
        length(c) < 13 && continue
        v = parse.(BigFloat, c[1:13])
        push!(rows, (av = Complex(v[1], v[2]), w1 = Complex(v[3], v[4]),
                     w2 = Complex(v[5], v[6]), w3 = Complex(v[7], v[8]),
                     d12 = v[9], d23 = v[10], d13 = v[11], dmax = v[12], tol = v[13]))
    end
    return rows
end

"""
    read_raw(path) -> Vector{Complex{BigFloat}}

All eigenvalues of one raw spectrum.
"""
function read_raw(path::AbstractString)
    vals = Complex{BigFloat}[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        a, b = split(s)
        push!(vals, Complex(parse(BigFloat, a), parse(BigFloat, b)))
    end
    return vals
end

"""
    read_retained(path) -> Vector{NamedTuple}

Rows of results/retained_roots.csv written by classify_modes.jl, in double precision:
m, k, j, mu, L (0 if unlabelled), n (-1 if undefined), class, w, w90, w95, w100, dmax, q,
and
the derived wminus = w - 3 pi mu and varpi = w - 3 pi mu / 2.
"""
function read_retained(path::AbstractString)
    rows = NamedTuple[]
    lines = readlines(path)
    for line in lines[2:end]
        c = split(line, ',')
        mu = parse(Float64, c[4])
        w = complex(parse(Float64, c[8]), parse(Float64, c[9]))
        push!(rows, (m = parse(Int, c[1]), k = parse(Int, c[2]), j = parse(Float64, c[3]),
                     mu = mu, L = parse(Int, c[5]), n = isempty(c[6]) ? -1 : parse(Int, c[6]),
                     class = c[7], w = w,
                     w90 = complex(parse(Float64, c[10]), parse(Float64, c[11])),
                     w95 = complex(parse(Float64, c[12]), parse(Float64, c[13])),
                     w100 = complex(parse(Float64, c[14]), parse(Float64, c[15])),
                     dmax = parse(Float64, c[16]), q = parse(Float64, c[17]),
                     wminus = w - 3pi * mu, varpi = w - 1.5pi * mu,
                     converged = c[7] == "converged"))
    end
    return rows
end

export read_retained

end # module
