# =====================================================================================
#  MatrixGenerator_Nu.mpl  --  Kerr 2D spectral matrices (s = 0), units 2M = 1
#
#  Rewrite of MatrixGenerator.mw.  Same operator (the L-coefficients derived in
#  MyExploration.mw), same Chebyshev-root grids in x and theta (= cos of the polar
#  angle), same ordering of unknowns and collocation rows, same file names
#      data/M2-<Nx>-<Nu>-<tag>.mat, M1R-..., M1I-..., M0R-..., M0I-...
#  so the pencil is still (M0 + w*M1 + w^2*M2) b = 0 with M1 = M1R + I*M1I etc.
#
#  Changes:
#   * the angular resolution N_u is a parameter and is raised together with N_x
#     (pairs in RES), so the three-resolution test also checks angular convergence;
#     N_u = 2 is exact only at a = 0 and caused the 0.04-3.5 % errors at a > 0;
#   * Chebyshev polynomials and their derivatives are evaluated numerically by
#     recurrence at the collocation points (no symbolic expansion, no coeff calls),
#     which is much faster for large N_x*N_u;
#   * Digits is lowered: Kerr does not need 200 digits (double precision already
#     reaches ~1e-10).  Set DECIMAL_DIGITS in the Julia solver to match.
#
#  Checked against the existing files: with RES = [[90, 2]] it reproduces the entries of
#  data/M*-90-2-3.mat to all 40 digits.
#
#  Usage (Maple):   restart;  read "MatrixGenerator_Nu.mpl";
#  or paste into a worksheet.  Run time: about 1.5 min per spin for the default RES;
#  the files take about 350 MB per spin at 40 digits (about 3.8 GB for all 11 spins).
# =====================================================================================
Digits := 40:      # enough for the Float64 or Double64 solver; 20 halves the file sizes for Float64

M := 0.5:          # 2M = 1
m := 1:

# (N_x, N_u) pairs of the three-resolution test.  The a/M = 0.99 overtones need N_x of
# about 80-100; N_u = 10 already resolves the angular direction to ~1e-10.
RES := [[60, 8], [70, 10], [80, 12]]:

# spin tags: a/M = tag/10, except tag 10 -> a/M = 0.99 (as in MatrixGenerator.mw)
SPINTAGS := [seq(k, k = 0 .. 10)]:

# Output folder for the matrices (absolute path, or relative to Maple's current directory);
# point it at the data folder next to QNMThreads_Nu.jl.  The files carry N_x and N_u in
# their names: data/M2-<Nx>-<Nu>-<tag>.mat, etc.
DATADIR := "data":

# ------------------------------------------------------------------------------------
#  Operator coefficients (verbatim from MatrixGenerator.mw / MyExploration.mw).
#  They use the globals a and m;  theta is the polar variable cos(vartheta).
# ------------------------------------------------------------------------------------
L200R := (x, theta) -> -256*(a - 1/2)^2*((-2*a^2 - x/4 + 3/4)*sqrt(-4*a^2 + 1) + a^2*x - 3*a^2 - x/4 + 3/4)*(a + 1/2)^2:
L111I := (x, theta) -> -2*(3*sqrt(-4*a^2 + 1)*x^2 + 16*a^2*x - 2*x*sqrt(-4*a^2 + 1) + 16*a^2 + x^2 - sqrt(-4*a^2 + 1) - 6*x - 3)*(2*a - 1)^3*(2*a + 1)^3/sqrt(-4*a^2 + 1):
L121R := (x, theta) -> -8*a*(theta^2 - 1)*(2*a - 1)^3*(2*a + 1)^3:
L100I := (x, theta) -> -2*(3*x*sqrt(-4*a^2 + 1) + 8*a^2 - sqrt(-4*a^2 + 1) + x - 3)*(2*a - 1)^3*(2*a + 1)^3/sqrt(-4*a^2 + 1):
L100R := (x, theta) -> 512*(a - 1/2)^3*(((-theta - 1)*m - theta)*sqrt(-4*a^2 + 1) + m*(-1 + x))*(a + 1/2)^3*a/sqrt(-4*a^2 + 1):
L012R := (x, theta) -> 2*(2*a - 1)^3*(1 + x)*(-1 + x)^2*(2*a + 1)^3:
L011R := (x, theta) -> 2*(2*a - 1)^3*(2*a + 1)^3*(-1 + x)*(3*x + 1):
L011I := (x, theta) -> 4*(2*a - 1)^3*(2*a + 1)^3*(-1 + x)^2*m*a/sqrt(-4*a^2 + 1):
L022R := (x, theta) -> -4*(2*a - 1)^3*(theta^2 - 1)*(2*a + 1)^3:
L021R := (x, theta) -> -8*(2*a - 1)^3*(m + 1)*theta*(2*a + 1)^3:
L000R := (x, theta) -> -2*(2*a - 1)^3*(2*a + 1)^3*(2*m^2 + 2*m - x + 1):
L000I := (x, theta) -> 4*(2*a - 1)^3*(2*a + 1)^3*m*(-1 + x)*a/sqrt(-4*a^2 + 1):

# ------------------------------------------------------------------------------------
#  Chebyshev roots and the values of T_n, T_n', T_n'' there (column n+1 = degree n)
# ------------------------------------------------------------------------------------
ChebRoots := proc(N::posint)
    [seq(evalf(cos((2*i - 1)*Pi/(2*N))), i = 1 .. N)]
end proc:

ChebTables := proc(pts::list, N::posint)
    local T, T1, T2, i, n, x, np;
    np := nops(pts);
    T  := Matrix(np, N);  T1 := Matrix(np, N);  T2 := Matrix(np, N);
    for i to np do
        x := pts[i];
        T[i, 1] := 1.0;  T1[i, 1] := 0.0;  T2[i, 1] := 0.0;
        if N > 1 then
            T[i, 2] := x;  T1[i, 2] := 1.0;  T2[i, 2] := 0.0;
        end if;
        for n from 2 to N - 1 do
            T[i, n + 1]  := 2*x*T[i, n] - T[i, n - 1];
            T1[i, n + 1] := 2*T[i, n] + 2*x*T1[i, n] - T1[i, n - 1];
            T2[i, n + 1] := 4*T1[i, n] + 2*x*T2[i, n] - T2[i, n - 1];
        end do;
    end do;
    return T, T1, T2;
end proc:

# ------------------------------------------------------------------------------------
#  Collocation matrices.  Row (i, j): x_i, theta_j  ->  r = (i-1)*Nu + j.
#  Column (p, q): T_p(x) T_q(theta)                 ->  c = p*Nu + q + 1.
# ------------------------------------------------------------------------------------
BuildMatrices := proc(Nx::posint, Nu::posint)
    global a;
    local X, TH, Tx, Tx1, Tx2, Tt, Tt1, Tt2, Ntot, Mat2, Mat1R, Mat1I, Mat0R, Mat0I,
          i, j, p, q, r, c, xx, tt, c200R, c100R, c100I, c121R, c111I, c000R, c000I,
          c021R, c022R, c011R, c011I, c012R, tx, tx1, tx2, t0, t1, t2, F;
    X  := ChebRoots(Nx);  TH := ChebRoots(Nu);
    Tx, Tx1, Tx2 := ChebTables(X, Nx);
    Tt, Tt1, Tt2 := ChebTables(TH, Nu);
    Ntot := Nx*Nu;
    Mat2  := Matrix(Ntot, Ntot);  Mat1R := Matrix(Ntot, Ntot);  Mat1I := Matrix(Ntot, Ntot);
    Mat0R := Matrix(Ntot, Ntot);  Mat0I := Matrix(Ntot, Ntot);
    for i to Nx do
        xx := X[i];
        for j to Nu do
            tt := TH[j];
            r := (i - 1)*Nu + j;
            c200R := evalf(L200R(xx, tt));  c100R := evalf(L100R(xx, tt));  c100I := evalf(L100I(xx, tt));
            c121R := evalf(L121R(xx, tt));  c111I := evalf(L111I(xx, tt));
            c000R := evalf(L000R(xx, tt));  c000I := evalf(L000I(xx, tt));
            c021R := evalf(L021R(xx, tt));  c022R := evalf(L022R(xx, tt));
            c011R := evalf(L011R(xx, tt));  c011I := evalf(L011I(xx, tt));  c012R := evalf(L012R(xx, tt));
            for p from 0 to Nx - 1 do
                tx := Tx[i, p + 1];  tx1 := Tx1[i, p + 1];  tx2 := Tx2[i, p + 1];
                for q from 0 to Nu - 1 do
                    c := p*Nu + q + 1;
                    t0 := Tt[j, q + 1];  t1 := Tt1[j, q + 1];  t2 := Tt2[j, q + 1];
                    F := tx*t0;
                    Mat2[r, c]  := c200R*F;
                    Mat1R[r, c] := c100R*F + c121R*tx*t1;
                    Mat1I[r, c] := c100I*F + c111I*tx1*t0;
                    Mat0R[r, c] := c000R*F + c021R*tx*t1 + c022R*tx*t2 + c011R*tx1*t0 + c012R*tx2*t0;
                    Mat0I[r, c] := c000I*F + c011I*tx1*t0;
                end do;
            end do;
        end do;
    end do;
    return Mat2, Mat1R, Mat1I, Mat0R, Mat0I;
end proc:

ExportAll := proc(Nx, Nu, tag, mats)
    local names, k, fname;
    names := ["M2", "M1R", "M1I", "M0R", "M0I"];
    for k to 5 do
        fname := sprintf("%s/%s-%d-%d-%d.mat", DATADIR, names[k], Nx, Nu, tag);
        ExportMatrix(fname, mats[k], target = MATLAB, mode = ascii);
    end do;
end proc:

# ------------------------------------------------------------------------------------
#  Main loop
# ------------------------------------------------------------------------------------
if not FileTools:-Exists(DATADIR) then FileTools:-MakeDirectory(DATADIR) end if:

for tag in SPINTAGS do
    if tag = 10 then a := evalf(0.99*M) else a := evalf(tag/10*M) end if;
    for nr in RES do
        t_start := time[real]();
        mats := [BuildMatrices(nr[1], nr[2])];
        ExportAll(nr[1], nr[2], tag, mats);
        printf("a/M = %.2f  N_x = %d  N_u = %d  done in %.1f s\n", a/M, nr[1], nr[2], time[real]() - t_start);
    end do;
end do:
