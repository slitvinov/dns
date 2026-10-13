FMM derivations (our own; papers used only as hints and cross-checks)
====================================================================

Every formula below was checked numerically against an independent direct
computation. The test that checked it is named in [brackets].

1. Harmonics and the addition theorem
-------------------------------------
Operators. Let xi = (x+iy)/2 and eta = -(x-iy)/2. Then
d_xi = d_x - i d_y and d_eta = -d_x - i d_y, so d_xi d_eta = -(d_x^2 + d_y^2).
On harmonic functions this gives d_xi d_eta = d_z^2.

Define D_n^m = (-d_z)^(n-|m|) d_eta^|m| for m >= 0, and d_xi^|m| in place
of d_eta^|m| for m < 0. On harmonic functions D_n^m D_k^l = D_(n+k)^(m+l):
every mixed d_xi d_eta pair turns into d_z^2.

Irregular harmonics: Thet_n^m(r) = D_n^m (1/r).
Regular harmonics: Ups_n^m(r) = r^(2n+1) Thet_n^m(r) / ((n-m)! (n+m)!), a
homogeneous harmonic polynomial of degree n (Kelvin inversion keeps
harmonicity).

Recurrences. These follow from D_1^m applied once more, and from
r^2 = z^2 - 4 xi eta:
  Thet_0^0 = 1/r
  Thet_n^n = (2n-1) (x+iy)/r^2 Thet_(n-1)^(n-1)
  r^2 Thet_n^m = (2n-1) z Thet_(n-1)^m - ((n-1)^2 - m^2) Thet_(n-2)^m
  Ups_0^0 = 1
  Ups_n^n = (x+iy)/(2n) Ups_(n-1)^(n-1)
  (n^2 - m^2) Ups_n^m = (2n-1) z Ups_(n-1)^m - r^2 Ups_(n-2)^m
Symmetry: Thet_n^(-m) = (-1)^m conj(Thet_n^m), and the same for Ups.
On the z axis: Thet_n^m(r e_z) = delta_m0 n!/r^(n+1) and
Ups_n^m(r e_z) = delta_m0 r^n/n!. [harm3 test, exact]

Addition theorem. Taylor-expand 1/|r - x| in x. Then
(x.grad)^n/n! = sum_m Ups_n^m(-x) conj(D_n^m) on harmonic functions, so
  1/|r - x| = sum_{n,m} Ups_n^m(x) conj(Thet_n^m(r)),   |x| < |r|.
[2e-16 at 50 random pairs, n <= 40]

Derivatives:
  d_z Ups_n^m = Ups_(n-1)^m
  d_eta Ups_n^m = Ups_(n-1)^(m+1)
  d_xi Ups_n^m = Ups_(n-1)^(m-1)
  so d_x Ups = (Ups_(n-1)^(m-1) - Ups_(n-1)^(m+1))/2
  and d_y Ups = i (Ups_(n-1)^(m-1) + Ups_(n-1)^(m+1))/2.
  For Thet: d_z Thet_n^m = -Thet_(n+1)^m, d_eta Thet_n^m = Thet_(n+1)^(m+1),
  d_xi Thet_n^m = Thet_(n+1)^(m-1).
  [gradients against the direct sum and finite differences]

2. Expansions and translations
------------------------------
P2M:  M_n^m(z) = sum_a q_a Ups_n^m(x_a - z)
M2P:  phi(x) = sum M_n^m conj(Thet_n^m(x - z))
P2L:  F_n^m(s) = sum_a q_a Thet_n^m(s - x_a)
L2P:  phi(x) = sum conj(Ups_n^m(s - x)) F_n^m
These all follow from the addition theorem with the right roles of x and r.

Translations. Apply exp(y.grad) = sum Ups_k^l(-y) conj(D_k^l) together with
D D = D:
M2M:  M_n^m(z') = sum_{k<=n,l} Ups_k^l(z - z') M_(n-k)^(m-l)(z)
      (exact for a truncated M) [2.5e-16]
L2L:  F_n^m(s') = sum_{k<=p-n,l} conj(Ups_k^l(s - s')) F_(n+k)^(m+l)(s)
      (exact for a polynomial) [1e-15]
M2L:  F_n^m(s) = sum_{k<=p,l} conj(M_k^l) Thet_(n+k)^(m+l)(s - z)
      Full truncation: n, k <= p, Thet up to degree 2p.
      For a real field conj(M_k^l) = (-1)^l M_k^(-l), so
      F = T M with T[(n,m),(k,l)] = (-1)^l Thet_(n+k)^(m-l)(t), t = s - z.
      That is the dense GEMM operator.
Scaling: Thet_j(h u) = h^-(j+1) Thet_j(u). So with M~_k = M_k h^-k and
F~_n = F_n h^(n+1), the operator T~ does not depend on the level.

3. Rotations
------------
Harmonics of degree n are closed under rotation:
Ups_n(R x) = G_n Ups_n(x), with G_n of size (2n+1) x (2n+1).
- z-rotation by a: G = diag(exp(i m a)), because phi -> phi + a.
- y-rotation by b: the generator is the angular-momentum J_y in the basis
  normalised by c_m = sqrt((n-m)!(n+m)!), so G_n = c^-1 d(b)^T c with
  d(b) = exp(-i b J_y).
  The convention (transpose, sign of b) was measured by least squares on
  random points, not taken from a paper. [1e-14 up to n = 40, b in (-pi, pi]]
- Rotation taking v to |v| e_z: R = R_y(-theta) R_z(-phi),
  G = G_y(-theta) diag(e^(-i m phi)).
Fast d(b): R_y(b) = R_x(-pi/2) R_z(b) R_x(pi/2), and R_x is R_y(pi/2)
conjugated by z-rotations of +-pi/2. So
  d(b)_{m m'} = Re( i^(m'-m) sum_k Delta_{k m} Delta_{k m'} e^(i k b) ),
  Delta = d(pi/2).
Delta is computed once per n, by a scaled-and-squared Taylor series of the
real antisymmetric -i J_y. [1.5e-14 at n = 29, including b = pi]
This replaces one matrix exponential per angle: 4.2 s -> 0.73 s for the
tables at p = 29.

Rotated translation along +z by rho:
  M2M:  M'_n^m = sum_{k=|m|}^n rho^(n-k)/(n-k)! M_k^m
  L2L:  F'_n^m = sum_{k=n}^p rho^(k-n)/(k-n)! F_k^m
  M2L:  F'_n^m = (-1)^m sum_{k=|m|}^p (n+k)!/rho^(n+k+1) M_k^m
For a multipole M -> G M. For a local, the coefficients transform as
F = G^H F' (from invariance of sum conj(Ups(s - x)) F), so F' = (G^-1)^H F.
For a real field only m >= 0 needs computing: X_n^(-m) = (-1)^m conj(X_n^m).
That is the half-spectrum kernel.

4. Mirror symmetries (option, off)
----------------------------------
z -> -z: P_n^m(-c) = (-1)^(n+m) P_n^m(c), so the coefficients pick up
(-1)^(n+m).
x -> -x: phi -> pi - phi, so c -> (-1)^m conj(c).
y -> -y: c -> conj(c).
Holds for multipole and local coefficients. [exact, 0.0]

5. Normalised, scaled harmonics (no factorials at run time)
-----------------------------------------------------------
Problem (measured): Thet_n^m carries (n-m)!(n+m)! ~ (2n)!. Our earlier double
factorial was formed with an int64 np.prod and overflowed at j = 18. A 20-level
tree at p = 30 gives h^60 ~ 1e-361, which underflows: the old code returned NaN
on that test (section 8).

Normalisation. Let c_nm = sqrt((n-m)!(n+m)!) and let Y_n^m be the Schmidt
semi-normalised surface harmonic (|Y| <= 1). Then
  Ups_n^m = r^n Y / c_nm  and  Thet_n^m = c_nm Y / r^(n+1).
For a box of size h define
  R_n^m(x) = c Ups_n^m(x) / h^n  = (r/h)^n Y,
  I_n^m(x) = h^(n+1) Thet_n^m(x) / c = (h/r)^(n+1) Y.
The c factors cancel in the addition theorem:
  1/|r - x| = (1/h) sum R_n^m(x) conj(I_n^m(r)),  any h.
  [2e-16 to 3.5e-16 at h = 1, 1e-2, 1e-6; degree 60]

Recurrences follow by substituting into those of section 1 and dividing out
c_n/c_(n-1) = sqrt(n^2 - m^2):
  R_0^0 = 1
  R_m^m = (x+iy)/h sqrt((2m-1)/(2m)) R_(m-1)^(m-1)
  R_n^m = [(2n-1)(z/h) R_(n-1)^m - sqrt((n-1)^2-m^2) (r/h)^2 R_(n-2)^m]
          / sqrt(n^2-m^2)
  I_0^0 = h/r
  I_m^m = (x+iy) h/r^2 sqrt((2m-1)/(2m)) I_(m-1)^(m-1)
  I_n^m = [(2n-1)(z h/r^2) I_(n-1)^m - sqrt((n-1)^2-m^2)(h/r)^2 I_(n-2)^m]
          / sqrt(n^2-m^2)
Every factor is O(1). [R = Ups c to 1.7e-16, I = Thet/c to 1.0e-16]

Scaled coefficients: M^ = sum q R_h(x - z) (that is M c / h^n) and
F^ = F h^(n+1) / c.
  M2P  phi = (1/h) sum M^ conj(I_h(x - z))
  L2P  phi = (1/h) sum conj(R_h(s - x)) F^
Gradients come from the derivative rules of section 1 times the c ratios:
  d_z R_n^m   =  sqrt(n^2-m^2)/h R_(n-1)^m
  d_eta R_n^m =  sqrt((n-m)(n-m-1))/h R_(n-1)^(m+1)
  d_xi R_n^m  =  sqrt((n+m)(n+m-1))/h R_(n-1)^(m-1)
  d_z I_n^m   = -sqrt((n+1)^2-m^2)/h I_(n+1)^m
  d_eta I_n^m =  sqrt((n+m+1)(n+m+2))/h I_(n+1)^(m+1)
  d_xi I_n^m  =  sqrt((n-m+1)(n-m+2))/h I_(n+1)^(m-1)
[direct sums: M2P/L2P potential and gradient to round-off or truncation, at
h = 1 and h = 1e-7]

Rotations in the normalised basis: G_hat = c G c^-1 = d^T diag(e^(-i l phi)).
This is unitary with entries bounded by 1, so multipoles and locals both
rotate with D, and back with D^H.

Axial translations (rho along +z; h_c/h_p = 1/2, rho/h_p = sqrt(3)/4 between
parent and child):
  M2M  M^'_n^m = sum_{k=|m|..n} sqrt(C(n-m,n-k) C(n+m,n-k))
                 (rho/h_p)^(n-k) (h_c/h_p)^k M^_k^m
  L2L  F^'_n^m = sum_{k=n..p} sqrt(C(k-m,k-n) C(k+m,k-n))
                 (rho/h_p)^(k-n) (h_c/h_p)^(n+1) F^_k^m
  M2L  F^'_n^m = (-1)^m sum_k sqrt(C(n+k,n+m) C(n+k,n-m))
                 (h_t/rho)^(n+1) (h_s/rho)^k M^_k^m
Each sqrt-binomial is <= 2^(n+k) and is multiplied by a power <= 2^-(n+k+1),
so every term is O(1).
Dense M2L: T^[(n,m),(k,l)] = (-1)^l sqrt(C(n+k-m+l, n-m) C(n+k+m-l, n+m))
I_(n+k)^(m-l)(-o), with I at unit scale. h cancels, so the operator does not
depend on the level.
[numpy fmmn.py at h = 1 and 1e-7: M2M and L2L exact to 1e-15, M2L to
truncation; CUDA = numpy FMM to round-off in all gate cases]

6. High-precision tables (computed once, rounded once)
-------------------------------------------------------
tables.py, mpmath at 50 digits:
- Delta_n = d(pi/2) = expm((pi/2) A_n), A_n = -i J_y real antisymmetric;
- the sqrt-binomial factors, from exact Python integers;
- the dense-M2L factor table.
fmm3.cu reads these and forms d(b) = Re(i^(m'-m) sum_k Delta Delta e^(ikb))
in double. That is a sum of 2n+1 terms of size O(1).

7. Periodic box
---------------
Lattice sums by Ewald splitting and Hobson's theorem. For a harmonic
polynomial P and a radial f, P(grad) f(r) = P(x) ((1/r) d/dr)^j f.
With Ups(grad)(1/r) = (-1)^j (2j-1)!! Ups(x)/r^(2j+1), this gives
  I_j^q(y) = (-1)^j/(2j-1)!! R_j^q(grad)(1/r).
Split 1/r = erfc(ar)/r + erf(ar)/r. The erf part is smooth; its lattice sum
is a reciprocal-space sum (the k = 0 mode dropped: neutralising background,
the same convention as the spectral solver).

Far sum for a target root cube and a source root offset d, excluding the
near images (cube offsets in [-1,1]^3):
  S_far = sum_far R(y) Bt_j(|y|)                                   (erfc)
        + (-1)^j i^j/(2j-1)!! (4pi/V) sum_k e^(-k^2/4a^2)/k^2
          R(k) e^(ik.x)                                           (recip)
        - sum_near R(y) E_j(|y|^2)                                 (erf)
        - pi/(V a^2) delta_j0
  Bt_j = B_j/(2j-1)!!, with
    Bt_j = (Bt_(j-1) + g_j e^(-a^2 r^2))/r^2,
    g_j = g_(j-1) 2a^2/(2j-1),  g_0 = 1/(a sqrt(pi)),
    Bt_0 = erfc(ar)/r.
  E_j(u) = (2a/sqrt(pi)) (2a^2)^j/(2j-1)!! sum_i (-a^2 u)^i/(i!(2j+2i+1)),
    from ((1/r)d/dr) = 2 d/du with u = r^2.
Precision: the earlier route "all images minus near images" cancels about
2^j relative (at j = 50 about 15 digits): measured floor near 1e-8 at
p >= 20. The split above never forms the near 1/r terms.
[far sums: brute-force direct sums agree to 2e-16 .. 1e-14 for j = 20..58;
alpha 1.0 vs 1.3 agree to <= 8e-14 relative for all j <= 58]

Non-harmonic part. P(x) = sum_R 1/|x - R| (regularised) satisfies
lap P = -4 pi sum delta + 4 pi/V, so P = P_h + (2 pi/(3V))|x|^2 with P_h
harmonic away from the lattice points. Solid-harmonic expansions represent
only P_h. So
- the tables use P_h: subtract (2pi/3V)|x|^2 at j = 0, and add
  (4pi/3V) R_1(x) at j = 1;
- the exact quadratic part is added once:
    phi += (2pi/3V) sum_j q_j |x - y_j|^2,
    grad += (4pi/3V)(Q x - d),  d = sum q_j y_j.
For a neutral cell this is the uniform "tinfoil" field from the cell dipole,
derived here rather than assumed.
[without it: cube potential error 43%, gradient 3e-3, not decreasing with p;
with it, the far field converges geometrically (numpy)]

Non-cubic cell 6pi x 4pi x 2pi: a forest of 3x2x1 root cubes. Image shifts
are carried in every list entry (box*27 + shift code). There is one far
operator per root offset d (15 of them).

Madelung constant (NaCl), our own Ewald:
  double precision:  1.7475645946331821  (alpha = 4, 6 agree)
  mpmath 32 digits:  1.747564594633182190636212036
                     (two independent (alpha, cutoff) sets agree to 28 digits)

8. Measured before/after (hal A100, GPU 3 idle; N = 5e5, s = 180)
-----------------------------------------------------------------
Unnormalised (v3) against normalised (fmm3n), p = 29:
  tables                       0.73 s -> 0.33 s
  total, cube (dense M2L)      1.02 s -> 0.61 s
  total, blob (dense M2L)      1.38 s -> 1.02 s
  error                        identical to 7 digits (2.5087e-11 cube,
                               5.1799e-11 blob)
  M2L, rotation kernel         0.236 s -> 0.281 s (table lookups; to tune)
20-level adaptive tree, p = 30, s = 8 (half the points in a 1e-6 blob):
  old: potential relL2 NaN (overflow/underflow)
  new: relL2 2.0e-11, max 2.2e-10, gradient max 3.2e-10;
       max|M| 37, max|L| 24, 0 non-finite, 0 subnormal (debug mode -g)
Periodic gate (Ewald reference):
  before (all-minus-near tables, unnormalised): floor ~1e-8 for p >= 20
  after (split tables, normalised):
    p   Madelung(shifted)  cube relL2  3x2x1 relL2  Biot 3x2x1 max
    20  6.6e-14            1.3e-8      1.5e-8       3.6e-8
    25  2.2e-15            1.0e-9      1.8e-9       6.7e-9
    29  2.2e-15            1.9e-10     3.1e-10      1.8e-9
  The Madelung cell with ions at the corners converges like 0.76^p
  (1.3e-6 at p = 29): source radius sqrt(3)/2 against far images at
  distance >= 2. Shifting the cell by (1/4, 1/4, 1/4) gives round-off.
Compensated (Neumaier) summation in P2P and in the direct-sum reference:
  measured effect <= 1e-17 absolute on the errors (N = 5e5, 2e6;
  p = 16, 29). Not needed; kept as compile options -DKP and -DKC.
Rotation matrices: double assembly from the 50-digit Delta against pure
mpmath expm(beta A): <= 1.7e-15 for n <= 40, including beta -> pi.
Integer safety: explicit checks for N, box*27 list codes, Morton key bits
(roots << 3D), and list sizes (sum taken in 64 bit before the 32-bit
scan). Combinatorics only in Python integers or mpmath.

9. Tube field against the spectral velocity; near-region size; M2L kernel
-------------------------------------------------------------------------
Setup (specref.py):
- tubes.py vorticity on 96 x 64 x 32 nodes of 6pi x 4pi x 2pi
  (h = 2pi/32, tube core radius 0.666 = 3.4 h);
- particles at the nodes with alpha = omega h^3 (mean removed);
- Gaussian blob cores, eps = 2h; kernel K q(r/eps),
  q = erf(rho) - 2 rho e^(-rho^2)/sqrt(pi).
Reference: the exact periodic velocity of the same blob sum by FFT. At the
nodes,
  u = sum_k0 e^(i k0.x) i A(k0) x alpha^(k0)/V,
  A(k0) = sum_aliases (k0+K) e^(-|k0+K|^2 eps^2/4)/|k0+K|^2,
with 125 alias shells; the first omitted term is ~e^(-9 pi^2).
This isolates the FMM error. The blob smoothing itself (a discretisation
choice) differs from the unregularised spectral velocity by 27% relL2 at
eps = 2h.
Leaves are 16h = 8 eps wide (s = 4096), so the point-kernel far field
adds an error ~e^(-64).
  p    max rel     rel L2      (K = 1 and K = 2 identical)
  20   6.302e-8    9.664e-9
  25   2.801e-9    3.516e-10
  29   6.101e-10   3.154e-11
  32   1.279e-10   7.933e-12
Random charges, periodic, against Ewald (gateK1/K2): K = 2 equals K = 1 to
3 digits at p = 20..32.
  Biot 3x2x1 at p = 32: max 5.8e-10, L2 3.5e-10
  potential cube at p = 32: L2 4.9e-11
Conclusion. The error is set by the one-box-separated M2L inside the tree
(worst pair ratio about sqrt(3)/2 over 2 - sqrt(3)/2), not by the
periodic far field. Widening the near region only moves offset-2 cubes from
the far table to explicit level-0 M2Ls with the same truncation. The
corner-ion Madelung cell, limited by root-level offset-2 M2L, is unchanged
too.
M2L rotation kernel (idle hal GPU 1, N = 5e5, s = 180), unnormalised v3
against normalised:
  table [m][n][k] (f0) or [n][k][m] (f1):  +18% at p = 29
  ratio recurrence with sqrt and divide (f2): +15%
  ratio table 1/sqrt((k-m+1)(k+m+1)) in L1 (f3, default):
    p 20 +17%, p 25 +15%, p 29 +14%
The old kernel ran an unbounded running product (n+k)! ir^(n+k+1) and
absorbed the c factors into huge or tiny coefficients. Matching it would
reintroduce those numbers. Dense-GEMM M2L (FMM_M2L=2) is unaffected:
0.207 s at p = 29. All variants give identical results (gate3: round-off).

10. Lattice near field, order-8 core, cost per particle
--------------------------------------------------------
3D order-8 core (derived with sympy, not taken from vm2):
  zeta(r) = pi^(-3/2) e^(-r^2) (105/16 - 63x/8 + 9x^2/4 - x^3/6), x = r^2
  radial moments: r^0 = 1; r^2, r^4, r^6 = 0; r^8 = -945/16
  Biot-Savart factor:
    q(rho) = 4 pi int_0^rho zeta r^2 dr
           = erf(rho) - (2/sqrt(pi)) rho e^(-rho^2)
             (1 - 89 rho^2/24 + 5 rho^4/3 - rho^6/6)
  Fourier symbol: zeta^(k) = e^(-kap^2)(1 + kap^2 + kap^4/2 + kap^6/6),
    kap = k delta/2 (mpmath quadrature agrees to 1e-31). This is 1 - O(k^8).
Core + quadrature error (grid sum against the spectral velocity of the true
vorticity, at the nodes; exact by FFT plus alias sum), delta = 2h:
  grid          h        max rel    rel L2
  96x64x32      0.196    3.9e-2     2.8e-2
  192x128x64    0.098    4.9e-3     2.4e-3
  384x256x128   0.049    3.0e-4     1.3e-4
  768x512x256   0.0245   1.3e-5     3.9e-6
  The ratio grows 8 -> 16 -> 23 (the tube core is 3.4h .. 27h; heading for
  256). delta = 1.5h: 1.6e-2, 1.8e-3, 7.8e-5. Gaussian (order 2): 3.5e-1,
  1.3e-1.
Lattice near field. After remeshing, the particles of a leaf block are the
B^3 cell-centred nodes. The near sum from a neighbour block at block offset o
is u = sum_c G_o^c alpha, with fixed B^3 x B^3 matrices per (o, component)
for every level, in units 1/h^2:
  G_o(i, j) = K_delta((l_i - oB - l_j) h)
The kernel values on the integer lattice come from mpmath (neartab.py).
  GEMM: 27 offsets x 3 DGEMMs (B^3 x B^3 times B^3 x 3 n_leaf), plus a
        cross-product combine.
  FFT:  per block, the (3B)^3 source patch in a (4B)^3 cyclic box (no wrap
        into the centre block), D2Z, multiply by the transformed kernel,
        Z2D, extract the centre.
  Both equal direct P2P with the same core: 4.9e-15 (GEMM), 9.2e-16
  (FFT vs GEMM).
  Offsets 2..Kc: the difference kernel (q - 1) K corrects the point-kernel
  far field for pairs inside the core tail. At B = 8 (gap 9h,
  rho = 4.5) it matters only below ~5e-8.
Timing (s per call; far = up + m2l + down + l2p; GEMM near; B = 8,
leaf = 512; tree build excluded, about 1 ns/particle at 1e8):
  A100 (hal, before contention), 192x128x64, N = 1.57e6:
    p 8 84 ms (53 ns/p), 12 117 ms (74), 16 190 ms (121), 20 328 ms (209)
    previous code (P2P near): 1.12 .. 1.46 s; the near field alone is 0.99 s
    against 0.027 s (37x)
  GH200 (hg), 384x256x128, N = 1.26e7:
    p 8 61 ns/p, 12 86, 16 124, 20 187 (Kc = 2: 234)
  GH200, 768x512x256, N = 1.01e8:
    p 8 26 ns/p, 12 37, 16 55, 20 84 (Kc = 2: 104)
    near 5.8 ns/p (GEMM), 9.2 ns/p (FFT)
  FFT is faster only at small N (96 grid: 3.6 ms vs 4.5 ms). GEMM wins from
  N ~ 1e6, where the DGEMMs have many columns.
  B = 16 is worse (near 115 ns/p at 192); B = 4 needs Kc ~ 4 with small
  GEMMs.
  Dense M2L (FMM_M2L=2) is slower than rotation for p <= 20.
FMM error (against the exact regularised-particle velocity), the same on all
grids:
  p 8 ~2e-4, 12 ~1e-5, 16 ~6e-7, 20 ~5e-8 (L2 ~9e-9)
Total error = core + FMM. p is chosen last (vm/APPROX.md): at 768^3-class
grids the core error is ~1e-5, so p = 12-14 is the balanced choice and
p = 16 makes the FMM negligible.
Cost model (GH200, N ~ 1e8, one velocity evaluation):
  "good" (p 12-15)          37-50 ns/particle
  "converged" (p 20)        84 ns   (104 with Kc = 2)
  FMM alone at 1e-8 (p 20-21, Kc 2)   ~0.1 us
  FMM alone at 1e-10 (p ~30, Kc 2)    estimate ~0.25 us
  The plan video assumed 2.4 us/particle: these are 25-65x lower ("good")
  and ~29x lower ("converged").

11. Lattice operators as fixed matrices (fmm3.cu -F, -G; near field on shifted views)
-------------------------------------------------------------------------------------
P2M and L2P on the lattice (leaftab.py, mpmath 40 digits, each entry rounded
once). A leaf of size h holds B^3 nodes at x_i - c = h (2l_i + 1 - B)/(2B),
l_i the local Morton-decoded node index. The scaled expansions (section 5) are
in units of h, so the maps do not depend on the level:
  P2M  M^_n^m = sum_i R_n^m((x_i - c)/h) q_i,  m < 0 from the mirror rule.
       As a real matrix: (2 (p+1)^2) x B^3, applied as one DGEMM to all
       leaves and all 3 components (columns leaf x component). The output is
       the stored M array itself (leaves are consecutive boxes, components
       consecutive per box).
  L2P  phi = (1/h) Re sum conj(R(s - x)) F^; the gradient uses the derivative
       rules of section 5, then u = curl(phi_1, phi_2, phi_3)/(4 pi). Each
       velocity component at node i is a real linear form in (Re, Im) of
       the 3 (p+1)^2 stored local coefficients, so L2P plus the curl is one
       (3 B^3) x (6 (p+1)^2) matrix; the factor 1/h^2 is the GEMM alpha.
       It adds straight into the near-field velocity array.
M2L on uniform levels (m2ltab.py: I_j^q(-o) at the 316 offsets o of
[-3,3]^3 minus [-1,1]^3, j <= 2p, mpmath 40 digits).
  Dense operator (section 5): T[(n,m),(k,l)] = (-1)^l fe I_(n+k)^(m-l)(-o).
  Fold the l < 0 columns with M_k^(-l) = (-1)^l conj(M_k^l): with
  T1 = T(nm,kl), T2 = (-1)^l T(nm,k,-l) (l > 0), M = a + ib,
    L = (T1 + T2) a + i (T1 - T2) b,
  a real (2 nh) x (2 nh) map, nh = (p+1)(p+2)/2: 0.54 of the complex flops.
  Children of one parent are stacked (8 x 2nh rows per parent cell and
  component). A target child q of parent P meets the source child q' of
  parent P + d at offset o = 2d + q' - q; it is a V-list pair iff
  max|o| >= 2. So for each of the 26 parent offsets d != 0 there is one
  (16 nh) x (16 nh) block matrix A_d (blocks T_o or 0; 1512 of 1664 blocks
  are nonzero), and
    L_stack(P) = sum_d A_d M_stack(P + d).
  With the parent grid padded by one cell (periodic images or zeros) and
  stored lexicographically, P + d is a constant column shift: 26 GEMMs per
  level with beta = 1 accumulation. The largest axis is innermost and the
  outer axis is the GEMM batch, so the halo columns cost (D+2)/D of the
  inner size only. One pack and one unpack per level.
  Requires every level to be full (uniform tree); level 0 with K > 1 stays
  on the rotation kernel.
Near field (lattice). u = sum_o (G_o^z alpha_y - G_o^y alpha_z, ...): each
kernel component multiplies only the two other alpha components (6 of 9
products). Leaves are packed once into a padded leaf grid (halo Kc,
[cell][component][B^3]); offset o is a column shift, the component choice is
a pointer offset with ldb = 3 B^3. The 6 products accumulate over all offsets
in the GEMMs (beta = 1) and one pass forms the cross product at the end.
[gates: velocity vs the previous code <= 4.8e-15 relative at 1.26e7, <= 9.7e-15 at 1e8 (
p 12-20); uniform non-lattice trees with -G <= 1e-17; errors vs the exact
regularised velocity unchanged]
