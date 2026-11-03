"""
    KerrData

Shared definitions for the Kerr benchmark (Sec. V): spin tags, readers for the output of
`QNMThreads_Nu.jl`, and the continued-fraction reference values of the qnm package.

Conventions: the solver works in units 2M = 1 and writes 2M*omega; everything here is
converted to M*omega.  Spin tag k = 0..9 is a/M = k/10 and tag 10 is a/M = 0.99.  The
computed spectrum has azimuthal number m = 1; it also contains the mirror images
-conj(omega) of the counter-rotating (m = -1) modes.
"""
module KerrData

using Printf

export TAGS, SPINS, a_of_tag, RESOLUTIONS, reslabel, read_diagnostics, read_raw,
       read_reference, diagnostics_file, raw_file, reference_value

const TAGS = 0:10
a_of_tag(k::Integer) = k == 10 ? 0.99 : k / 10
const SPINS = a_of_tag.(TAGS)
const RESOLUTIONS = ((60, 8), (70, 10), (80, 12))
reslabel(r) = "N$(r[1])-M$(r[2])"

diagnostics_file(root, k) = joinpath(root, "J$k", "stable-damped-modes-diagnostics-J$k.txt")
raw_file(root, k, r) = joinpath(root, "J$k", "raw-modes-$(reslabel(r))-J$k.txt")

"""
    read_diagnostics(path) -> Vector{NamedTuple}

Converged roots (M*omega): average `av`, the three matched roots, and `dmax` (also M*omega).
"""
function read_diagnostics(path)
    rows = NamedTuple[]
    isfile(path) || return rows
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        c = parse.(Float64, split(s)[1:13])
        push!(rows, (av = complex(c[1], c[2]) / 2, w1 = complex(c[3], c[4]) / 2,
                     w2 = complex(c[5], c[6]) / 2, w3 = complex(c[7], c[8]) / 2,
                     dmax = c[12] / 2))
    end
    return rows
end

"All eigenvalues of one raw spectrum, as M*omega."
function read_raw(path)
    vals = ComplexF64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        a, b = split(s)
        push!(vals, complex(parse(Float64, a), parse(Float64, b)) / 2)
    end
    return vals
end

"""
    read_reference(path) -> Dict{Tuple{Int,Int,Int,Float64},ComplexF64}

qnm-package values M*omega keyed by (l, m, n, a/M); see reference/qnm_reference.py.
"""
function read_reference(path)
    ref = Dict{Tuple{Int,Int,Int,Float64},ComplexF64}()
    for line in readlines(path)[2:end]
        c = split(line, ',')
        ref[(parse(Int, c[1]), parse(Int, c[2]), parse(Int, c[3]), parse(Float64, c[4]))] =
            complex(parse(Float64, c[5]), parse(Float64, c[6]))
    end
    return ref
end

"""
    reference_value(ref, l, m, n, a)

Reference frequency as it appears in the computed m = 1 spectrum: omega for m = 1 and the
mirror image -conj(omega) for m = -1.
"""
reference_value(ref, l, m, n, a) = (w = ref[(l, m, n, a)]; m == 1 ? w : -conj(w))

end # module
