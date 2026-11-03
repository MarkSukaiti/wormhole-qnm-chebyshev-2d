# =====================================================================================
#  static_parity.jl
#
#  Parity of the static (j = 0) quasinormal modes under eta -> -eta (Sec. IV B).  For
#  j = 0 the radial equation is even in eta, so each mode is either even (u'(0) = 0) or
#  odd (u(0) = 0).  The solution continued from eta = -inf is evaluated at the throat by
#  direct integration; |u(0)|/|u'(0)| is large for even and small for odd modes.  The
#  overtones of each family alternate in parity.  The test is decisive for n <= 3; for more
#  strongly damped overtones the ratio is O(1) in double precision and is not used.
#
#  Usage:  julia --project=.. static_parity.jl  [results]
# =====================================================================================

using Printf

include("WormholeData.jl")
include("DirectIntegration.jl")
using .WormholeData
using .DirectIntegration

function main(resdir = "results")
    R = read_retained(joinpath(resdir, "retained_roots.csv"))
    seen = Set{Tuple{Int,Int}}()
    for o in sort([o for o in R if o.k == 0 && o.converged && 0 <= o.n <= 3 && real(o.w) > 0];
                  by = o -> (o.L, o.n))
        (o.L, o.n) in seen && continue            # the same static mode appears for every m <= L
        push!(seen, (o.L, o.n))
        r = parity_ratio(o.w, o.L)
        @printf("L=%d n=%d  w = %.6f%+.6fi   |u(0)|/|u'(0)| = %.2e   %s\n", o.L, o.n,
                real(o.w), imag(o.w), r, r > 1 ? "even" : "odd")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? "results" : ARGS[1])
end
