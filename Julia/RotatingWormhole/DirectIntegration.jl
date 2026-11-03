"""
    DirectIntegration

Independent computation of the scalar quasinormal modes of the spinning Morris-Thorne
wormhole by direct integration of the radial equation (Sec. IV A of the paper).  No
compactification, no asymptotic factorisation and no spectral discretisation is used.

After the substitution  Psi = (1+eta^2)^(-1/2) u(eta) P_L^m(chi),  the radial equation is

    u'' + Q(eta) u = 0,
    Q(eta) = (w - m*Omega(eta))^2 - L(L+1)/(1+eta^2) - 1/(1+eta^2)^2,
    m*Omega(eta) = 3 mu [pi/2 - arctan(eta) - eta/(1+eta^2)],      mu = m j,

with  eta = ell/r0,  w = r0*omega,  and the outgoing conditions

    u ~ exp(+i w eta)        (eta -> +inf),
    u ~ exp(-i w_- eta)      (eta -> -inf),     w_- = w - 3 pi mu.

Each outgoing solution is started from the exact Riccati-Hankel function at a complex
point on a ray along which it is subdominant (so the start-up data are exponentially
insensitive to the neglected far-field terms) and is integrated back to the real axis and
then to the matching point with an adaptive eighth-order Runge-Kutta scheme (DP8).  The
quasinormal condition is the equality of the two logarithmic derivatives at the matching
point; roots are found with a damped secant iteration.

The module also solves the large-mu limit problem of Sec. IV E,

    u_yy + [(w + 2/|y|^3)^2 - L(L+1)/y^2] u = 0,   y < 0,

with transmission into the throat region as y -> 0^- and outgoing behaviour as y -> -inf;
its eigenvalues w_{L,n} give  w_- ~ w_{L,n} mu^(-1/2).

All computations are in double precision; fundamental modes and first overtones are
accurate to about 1e-11.  For strongly damped overtones keep the rays away from the
branch points eta = +-i of arctan.
"""
module DirectIntegration

using OrdinaryDiffEqHighOrderRK: ODEProblem, solve, DP8

export qnm, qnm_limit, parity_ratio, mismatch, frame_dragging

const RTOL = 1e-12
const ATOL = 1e-30

# ---------------------------------------------------------------------------------------
# Special functions and the radial potential
# ---------------------------------------------------------------------------------------

"""
    riccati_hankel(L, z) -> (xi, dxi)

Riccati-Hankel function xi_L(z) = z h_L^(1)(z) (up to a constant factor, which cancels in
logarithmic derivatives) and its derivative, from the exact finite sum.
"""
function riccati_hankel(L::Integer, z::Complex)
    s = zero(z)
    ds = zero(z)
    for k in 0:L
        a = factorial(big(L + k)) / (big(2)^k * factorial(big(k)) * factorial(big(L - k)))
        c = (1.0im)^(k - L - 1) * Float64(a)
        s += c * z^(-k)
        ds += c * (1im * z^(-k) - k * z^(-k - 1))
    end
    e = exp(1im * z)
    return e * s, e * ds
end

"""
    frame_dragging(eta)

pi/2 - arctan(eta) - eta/(1+eta^2), continued analytically into the complex eta-plane
(branch cuts on the imaginary axis, |Im eta| >= 1).  m*Omega_hat = 3 mu * frame_dragging(eta).
"""
frame_dragging(eta) = pi / 2 - atan(eta) - eta / (1 + eta * eta)

function Q_full(eta, w, mu, L)
    k = w - 3.0 * mu * frame_dragging(eta)
    r2 = 1 + eta * eta
    return k * k - L * (L + 1) / r2 - 1.0 / (r2 * r2)
end

function Q_limit(y, w, L)
    k = w - 2.0 / y^3
    return k * k - L * (L + 1) / (y * y)
end

# ---------------------------------------------------------------------------------------
# Integration along a polygonal path in the complex plane
# ---------------------------------------------------------------------------------------

"""
    integrate_path(Q, pts, u0, du0, args) -> (u, du)

Integrate u'' = -Q(z) u along the straight segments pts[1] -> pts[2] -> ... and return
u and du/dz at the end point.  The state is renormalised on every segment.  The complex
system is integrated as a real system for (Re u, Im u, Re u', Im u'), so that the error
control acts on each real component (as in DOP853); with a complex state vector the
step-size control is not reliable for the rapidly growing solutions met at large |w|.
"""
function integrate_path(Q, pts, u0, du0, args)
    u, du = ComplexF64(u0), ComplexF64(du0)
    for s in 1:length(pts)-1
        a, b = ComplexF64(pts[s]), ComplexF64(pts[s+1])
        d = b - a
        sc = max(abs(u), abs(du), 1e-300)
        function f(Y, p, t)
            y, dy = complex(Y[1], Y[2]), complex(Y[3], Y[4])
            r1 = d * dy
            r2 = -d * Q(a + t * d, args...) * y
            return [real(r1), imag(r1), real(r2), imag(r2)]
        end
        y0 = [real(u / sc), imag(u / sc), real(du / sc), imag(du / sc)]
        sol = solve(ODEProblem(f, y0, (0.0, 1.0)), DP8(); reltol = RTOL, abstol = ATOL,
                    save_everystep = false, maxiters = 10^7)
        Y = sol.u[end]
        u, du = complex(Y[1], Y[2]), complex(Y[3], Y[4])
    end
    return u, du
end

"""
    best_ray(arg, lo, hi, sign) -> (angle, s)

Ray direction in [lo, hi] on which exp(i*sign*|w| e^{i arg} t e^{i angle}) decays, i.e.
s = sign*sin(arg + angle) > 0; among sufficiently decaying directions, the one closest
to the real axis is chosen.  Returns (nothing, smax) if there is none.
"""
function best_ray(arg, lo, hi, sign)
    angs = range(lo, hi; length = 721)
    sv = sign .* sin.(arg .+ angs)
    smax = maximum(sv)
    smax <= 0 && return nothing, smax
    target = min(0.5, 0.8 * smax)
    axis = sign > 0 ? 0.0 : pi
    ok = findall(>=(target), sv)
    i = ok[argmin(abs.(angs[ok] .- axis))]
    return angs[i], sv[i]
end

"""
    paths(w, mu; eta_m = 0.5, supp = 30.0)

Integration paths for the outgoing solutions at the two ends, frozen for a root search.
The start points are chosen so that the subdominant solution is suppressed by e^{-supp}
relative to the dominant one.
"""
function paths(w, mu; eta_m = 0.5, supp = 30.0)
    wm = w - 3pi * mu
    d = 0.005
    phi, sR = best_ray(angle(w), -pi / 2 + d, pi / 2 - d, +1)      # first end
    psi, sL = best_ray(angle(wm), pi / 2 + d, 1.5pi - d, -1)        # second end
    (phi === nothing || psi === nothing) && return nothing, nothing, -1.0, -1.0
    cR, cL = abs(w) * sR, abs(wm) * sL
    TR, TL = max(supp / cR, 1.0), max(supp / cL, 1.0)
    aR = abs(w) < 30 ? 3.0 : 1.0
    right = ComplexF64[aR + TR * cis(phi), aR, eta_m]
    left = ComplexF64[-3.0 + TL * cis(psi), -3.0, eta_m]
    return right, left, cR, cL
end

"""
    solutions(w, mu, L, P) -> (uR, duR, uL, duL)

Outgoing solutions at the first (R) and second (L) end, continued to the matching point.
"""
function solutions(w, mu, L, P)
    wm = w - 3pi * mu
    right, left = P
    o, do_ = riccati_hankel(L, w * right[1])
    uR, duR = integrate_path(Q_full, right, o, w * do_, (w, mu, L))
    o, do_ = riccati_hankel(L, -wm * left[1])
    uL, duL = integrate_path(Q_full, left, o, -wm * do_, (w, mu, L))
    return uR, duR, uL, duL
end

"""
    mismatch(w, mu, L, P)

Difference of the logarithmic derivatives of the two outgoing solutions at the matching
point; it vanishes at a quasinormal frequency.
"""
function mismatch(w, mu, L, P)
    uR, duR, uL, duL = solutions(w, mu, L, P)
    return duR / uR - duL / uL
end

"""
    secant(F, w0; w1, tol, maxit) -> (w, last_step)

Damped secant iteration for a complex root of F.
"""
function secant(F, w0; w1 = w0 * (1 + 1e-4) + 1e-4im, tol = 1e-11, maxit = 60)
    f0, f1 = F(w0), F(w1)
    for _ in 1:maxit
        f1 == f0 && break
        w2 = w1 - f1 * (w1 - w0) / (f1 - f0)
        lim = 0.3 * max(1.0, abs(w1))
        if abs(w2 - w1) > lim                 # damp very large steps
            w2 = w1 + lim * (w2 - w1) / abs(w2 - w1)
        end
        w0, f0 = w1, f1
        w1, f1 = w2, F(w2)
        abs(w1 - w0) < tol * max(1.0, abs(w1)) && return w1, abs(w1 - w0)
    end
    return w1, abs(w1 - w0)
end

"""
    qnm(mu, L, guess; eta_m = 0.5) -> (w, last_step)

Quasinormal frequency w = r0*omega of the multipole L at mu = m j, starting from `guess`.
The frequency in the frame of the second exterior is w - 3 pi mu.
"""
function qnm(mu::Real, L::Integer, guess::Number; eta_m = 0.5)
    right, left, cR, cL = paths(ComplexF64(guess), mu; eta_m = eta_m)
    (cR <= 0 || cL <= 0) && error("no subdominant direction for w = $guess")
    P = (right, left)
    return secant(w -> mismatch(w, mu, L, P), ComplexF64(guess))
end

"""
    parity_ratio(w, L)

Static case only: |u(0)|/|u'(0)| of the solution continued from eta = -inf.  A small
value means an odd mode (u(0) = 0), a large value an even mode (u'(0) = 0).
"""
function parity_ratio(w::Number, L::Integer)
    right, left, _, _ = paths(ComplexF64(w), 0.0; eta_m = 0.0)
    o, do_ = riccati_hankel(L, -w * left[1])
    u, du = integrate_path(Q_full, left, o, -w * do_, (ComplexF64(w), 0.0, L))
    return abs(u) / abs(du)
end

# ---------------------------------------------------------------------------------------
# Large-mu limit problem (Sec. IV E)
# ---------------------------------------------------------------------------------------

function mismatch_limit(w, L; rs = 0.2, rho = 1.0, supp = 30.0)
    # near end: the transmitted solution is subdominant on arg y = 3 pi / 4 as |y| -> 0
    th = 0.75pi
    ys = rs * cis(th)
    k = w - 2.0 / ys^3
    kp = 6.0 / ys^4
    u0 = k^(-0.5) * exp(1im * (w * ys + 1.0 / ys^2))
    du0 = u0 * (1im * k - 0.5 * kp / k)
    arc = [rho * cis(t) for t in range(th, pi; length = 6)]
    near = vcat(ComplexF64[ys], arc)
    uN, duN = integrate_path(Q_limit, near, u0, du0, (w, L))
    # far end: outgoing exp(-i w y), subdominant along direction psi
    psi, sL = best_ray(angle(w), pi / 2 + 0.005, 1.5pi - 0.005, -1)
    psi === nothing && error("no subdominant direction for w = $w")
    c = abs(w) * sL
    T = max(supp / c, 10.0)
    yf = -rho + T * cis(psi)
    o, do_ = riccati_hankel(L, -w * yf)
    uF, duF = integrate_path(Q_limit, ComplexF64[yf, -rho], o, -w * do_, (w, L))
    return duN / uN - duF / uF
end

"""
    qnm_limit(L, guess) -> (w_{L,n}, last_step)

Eigenvalue of the large-mu limit problem, Eqs. (limit_equation)-(limit_bc) of the paper.
"""
qnm_limit(L::Integer, guess::Number) = secant(w -> mismatch_limit(w, L), ComplexF64(guess))

end # module
