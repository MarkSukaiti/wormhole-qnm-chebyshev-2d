# =====================================================================================
#  radial_pencils_arb.jl
#
#  Independent arbitrary-precision solution of the decoupled radial pencils of
#  Proposition III.1 (Eq. radial_pencils of the paper),
#
#      R0 + L(L+1) B_x + w R1 + w^2 R2,      L = |m|, |m|+1, |m|+2,
#
#  at N_x = 90, on the same Chebyshev-Gauss grid and with the same normalisation as the
#  Maple matrix generator, in midpoint-radius (ball) arithmetic with 256-bit precision
#  (Arb, via Arblib.jl).  The companion linearisation of each pencil is diagonalised with
#  acb_mat_approx_eig_qr.  By Proposition III.1 the union of the three radial spectra must
#  coincide with the raw two-dimensional spectrum computed by QNMThreads.jl; the
#  comparison and the multipole labels are made in classify_modes.jl.
#
#  Usage:
#      julia --project=.. -t 4 radial_pencils_arb.jl [output_folder]
#  Output: <output_folder>/radial-m<m>-J<k>-L<L>-N90.txt   (default folder: radial_spectra)
#  Run time: about 15 s per pencil, 99 pencils.
# =====================================================================================

using Arblib
using Printf
using Base.Threads

include("WormholeData.jl")
using .WormholeData

const PREC = 256
const NX = 90

function radial_matrices(m::Integer, k::Integer, L::Integer, N::Integer; prec = PREC)
    jr = j_of_tag(k)
    j = Arb(numerator(jr); prec) / Arb(denominator(jr); prec)
    piA = Arb(π; prec)
    A0 = AcbMatrix(N, N; prec)
    A1 = AcbMatrix(N, N; prec)
    A2 = AcbMatrix(N, N; prec)
    for i in 1:N
        th = (2i - 1) * piA / (2N)
        x = cos(th)
        st = sin(th)
        c = cos(piA * x / 2)
        s = sin(piA * x / 2)
        c2 = c * c
        Om = Arb(3; prec) / 2 * j * (piA * (1 - x) - sin(piA * x))   # m-independent profile
        Om0 = Arb(3; prec) / 2 * piA * j                             # throat value
        H = s * (1 + c2)
        l2 = 1 - c2 - c2 * c2
        l1 = Acb(2m / c2 * (Om + H * (1 - H) * Om0), c * (1 - 3c2); prec)
        l0 = Acb(m^2 / c2 * ((1 - H)^2 * Om0^2 - Om^2) + c2 + L * (L + 1),
                 -c * (1 - 3c2) * m * Om0; prec)
        d10 = Acb(2 / piA * sin(piA * x), -4 / piA * (1 - H) * m * Om0; prec)
        d11 = Acb(0, -4 / piA * H; prec)
        d2 = -4 / piA^2 * c2
        for p in 0:N-1
            T = cos(p * th)
            Tp = p * sin(p * th) / st
            Tpp = (x * Tp - p^2 * T) / (st * st)
            A2[i, p+1] = l2 * T
            A1[i, p+1] = l1 * T + d11 * Tp
            A0[i, p+1] = l0 * T + d10 * Tp + d2 * Tpp
        end
    end
    return A0, A1, A2
end

"""
    radial_spectrum(m, k, L; N = NX) -> AcbVector

Eigenvalues of the quadratic pencil A0 + w A1 + w^2 A2 from its companion linearisation.
"""
function radial_spectrum(m, k, L; N = NX, prec = PREC)
    A0, A1, A2 = radial_matrices(m, k, L, N; prec)
    A2i = inv(A2)
    B0 = -(A2i * A0)
    B1 = -(A2i * A1)
    C = AcbMatrix(2N, 2N; prec)
    for i in 1:N
        C[i, N+i] = Acb(1; prec)
        for q in 1:N
            C[N+i, q] = B0[i, q]
            C[N+i, N+q] = B1[i, q]
        end
    end
    ev = AcbVector(2N; prec)
    Arblib.approx_eig_qr!(ev, C; prec)
    return ev
end

function main(outdir = "radial_spectra")
    mkpath(outdir)
    jobs = [(m, k, L) for m in SECTORS for k in TAGS for L in (m, m + 1, m + 2)]
    lk = ReentrantLock()
    @threads for job in jobs
        m, k, L = job
        t0 = time()
        ev = radial_spectrum(m, k, L)
        fn = joinpath(outdir, "radial-m$m-J$k-L$L-N$NX.txt")
        open(fn, "w") do io
            println(io, "# radial pencil eigenvalues, m=$m, spin tag $k (j=$(j_label(k))), L=$L, N_x=$NX, prec=$PREC bits")
            for z in ev
                re = string(Arblib.midref(Arblib.realref(z)); digits = 40, remove_trailing_zeros = false)
                im = string(Arblib.midref(Arblib.imagref(z)); digits = 40, remove_trailing_zeros = false)
                println(io, re, " ", im)
            end
        end
        lock(lk) do
            @printf("m=%d k=%2d L=%d: %d eigenvalues (%.1f s)\n", m, k, L, length(ev), time() - t0)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? "radial_spectra" : ARGS[1])
end
