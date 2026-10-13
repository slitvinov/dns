vm2: 2D remeshed vortex particles, periodic FMM, adaptive blocks, MRI LTS
=========================================================================

Files: vm2.cu (solver), sp2.cu (spectral reference), gen_vm2.py (op tables,
high-precision coefficient tables, tests), Makefile rules sp2 vm2 vm2_g8.bin.
Build on hal: /usr/local/cuda-12.8/bin/nvcc -O3 -arch=sm_80 -default-stream per-thread
-Xcompiler -ffp-contract=off -o vm2 vm2.cu -lcufft (the last flag keeps host and GPU
regrid bitwise equal on hosts with FMA, e.g. aarch64)
Problem: double shear layer on [0, 2 pi)^2, omega = eps cos x - u'(y), u' a
periodic sum of sech^2 layers (n = -3..3), rho = pi/15, eps = 0.05.

1. Derivations (first principles, each checked numerically)
-----------------------------------------------------------
Remeshing kernel. W(y) even, piecewise degree 7 on |y| in [j, j+1], j < 4,
conditions: W(j) = delta_j0 (interpolating), C^3 at the knots, and
sum_k (t - k)^m W(t - k) = delta_m0 for m < 7 as polynomial identities in t.
Unique solution, rational coefficients (table K8 in vm2.cu). Moment errors
<= 2e-13 in double; remesh reproduces polynomials of degree <= 6 (test1.py).
Remeshing every step loses one order: global order 6 (1D: 5.8, 2D gate 1: 5.9).
Exact conservation in floating point: the partition of unity sum_k W holds only
to ~1e-15, so each particle's weights are divided by (sum_x W)(sum_y W),
computed separably (16 evaluations).

Wavelet prediction (DD8). Fine midpoint value = Lagrange interpolant through
8 coarse points: weights (-5, 49, -245, 1225, 1225, -245, 49, -5)/2048, exact
for degree <= 7; tensor product in 2D by separable x then y passes. Weights
sum exactly to 1 in binary, so every ghost row of the composed fill matrix G
sums to 1.

Conservation by transposed tables. Lab of block b = own nodes plus ghosts at
b's level, ghost = G * (neighbour nodes) with G composed from COPY (same
level), INJ (coarse ghost from fine node), PX/PY (fine ghost from coarse
buffer). Own particles deposit masses onto the whole lab; ghost masses go back
with G^T. Then sum_n (G^T m)_n = sum_k m_k sum_n G_kn = sum_k m_k. Circulation is
conserved to round-off by construction (gate 3).

Node volumes. Integral of the block-wise DD8 interpolant: V = G^T(h^2 q_i q_j)
with q_i = int_0^B phi(x - i) dx, phi the DD8 fundamental function. q is exact
(rational) from the refinement equation phi(x) = sum_k c_k phi(2x - k):
J(a) = int_-inf^a phi satisfies J(a) = 1/2 sum_k c_k J(2a - k), J = 0 for
a <= -7, 1 for a >= 7 (gen_vm2.py -q). sum q = 16 and q integrates degree <= 7
exactly. Checks: sum V = area to 1e-12; quadrature of smooth trig products to
1e-12 on two-level meshes. Gamma = omega V, omega = Gamma / V.

Velocity (near field). Plain point vortices on the grid (punctured trapezoid of
the singular Biot-Savart integral): measured O(h^2). Regularized kernel
K_delta = K * zeta_delta with zeta = exp(-r^2/delta^2) * poly(r^2): Fourier
symbol of the error is 1 - zeta^(k delta), so the smoothing error is
O((k delta)^order) when the first order/2 - 1 radial moments of zeta vanish;
the grid sum (trapezoid of a smooth integrand) adds only
exp(-pi^2 (delta/h)^2) = 3e-18 at delta = 2h. q(rho) = 1 - exp(-x) R(x),
x = rho^2, R from the moment equations (mpmath, gen_vm2.py -l; order 8:
R = 1 - 3x + 3/2 x^2 - x^3/6). Far field: point-kernel expansions are exact up
to exp(-x) poly(x) with x >= 64 because every well-separated source box is at
least 16 particle spacings away and delta = 2h.
Measured against the FFT velocity of the same nodal vorticity (max rel error):
    n      point   gauss(2)  4        6        8
    64     1.3e-2  1.2e-1    3.4e-2   1.3e-2   5.5e-3
    128    3.4e-3  3.8e-2    4.3e-3   6.7e-4   1.3e-4
    256    8.4e-4  1.0e-2    3.5e-4   1.8e-5   1.2e-6
    512    2.1e-4  2.6e-3    2.3e-5   3.3e-7   6.4e-9
Ratios per doubling tend to 4, 4, 16, 64, 256: orders 2, 2, 4, 6, 8. Positive
(Gaussian) cores are second order. Order 8 is used by default.

Periodic FMM. Tree = block quadtree; lists with image shifts (3x3 root images
covered by U, V, W, X lists); images beyond via a root local expansion with
lattice sums S_s = sum_{far} w^-s. In the zeta-regularised sum the s = 1, 2
terms vanish; S_s = 0 unless s = 0 mod 4; G_s from G4 = Gamma(1/4)^8/(960 pi^2)
by the Weierstrass recurrence (g3 = 0), ring of the 8 near images subtracted
exactly (mpmath; S8 agrees with a direct sum to 1.6e-17). The periodic kernel
is zeta(z) - (pi/A) zbar, which gives the extra uniform term
(i/2A)(Gamma_tot zbar - conj D), D = sum Gamma z. The resulting velocity has
zero mean and matches the spectral velocity (table above). FMM vs direct sum
over the 3x3 images plus root local: 1.5e-11 (uniform), 1.8e-8 (two levels).
p = 16 is enough below 1e-8 (p = 12: 1.2e-5 floor at n = 256).

Diffusion. PSE inside the ODE (no splitting): dGamma_i/dt =
nu eps^-4 sum_j (V_i Gamma_j - V_j Gamma_i) eta(r/eps), eta = exp(-x) poly(x)
with int rho^2 eta = 4 and higher even moments zero (order 8, mpmath).
eps_ij = 2 max(h_i, h_j). Antisymmetric pairs conserve circulation exactly:
reverse X-list pairs are added for PSE only (otherwise circulation drifted
1e-10 per unit time on two-level meshes).

Time stepping. Heun RK3 (c = 0, 1/3, 2/3; b = 1/4, 0, 3/4). MRI-GARK ERK33a
(coefficients G0, G1, dc = 1/3 as in ~/dns/vm/lts2.py; scalar test: local error
O(H^4) for any number of inner steps). Fast/slow split: fast = blob + PSE pairs
between fast particles whose ORIGIN blocks are U-neighbours (fixed during the
macro step); slow = full FMM minus that. A position-based split (current
leaves) is not smooth in time, because pairs switch between the fast and the
frozen slow part, and gives an O(H) error. Rung per block from
max |grad u| > f max (central differences on the lab) plus a one-block buffer.

2. Gates (hal GPU 3, A100)
--------------------------
G0 tables: gen_vm2.py -t: 1168 random neighbour configurations per ghost width
   (levels, parities, tensor degrees <= 7), every ghost cell exact (Fraction).
G1 uniform vm2 vs sp2 (N = 1024, IFRK4 dt 1e-4), t = 1, nu = 1e-4,
   H = 0.005, M = 2, relative L2 of omega:
       n = 64 4.9e-3, 128 1.75e-4, 256 3.55e-6, 512 5.93e-8
   ratios 28, 49, 60 -> order 6 (remeshing order). Halving H: 5.932e-8 ->
   5.941e-8 (time error negligible). Circulation 5e-14.
G2 adaptive base 256 with every block refined once == uniform 512: every
   output column bitwise identical (except wall time).
G3 two-level mesh (base 128, band 0.4 refined): circulation 3e-14 -> 2e-13
   after 200 steps (round-off) with transposed tables + symmetric PSE.
G4 MRI with all blocks on one rung vs plain Heun RK3: Z and E equal to 5e-15.
G5 thin layer rho = 0.1047, n = 256, t = 0.2, M = 12, error vs same-H
   reference (Heun, M = 96):
       H      shared M=12  MRI M=12   shared M=1
       0.04   1.5e-10      5.08e-6    2.6e-7
       0.02   1.9e-11      6.32e-7    3.3e-8
       0.01   2.7e-12      7.88e-8    4.1e-9
       0.005  1.5e-12      9.82e-9    5.1e-10
   MRI order 3 (ratio 8). No speedup on this problem: the fast group is 50% of
   the particles (layers span the domain; 1-block buffer) and the slow part
   limits the error; MRI's error is 19x that of plain Heun with the same macro
   step. Cost (pair evaluations) MRI/shared(M=12) = 0.52.

3. Follow-up items 1-5
---------------------
3.1 Level-border fix (not solved)
  Smooth 8th-order restriction for coarse ghosts (default; -j = old injection):
  coarse ghost = sum (c_a/2)(c_e/2) F(2i+a, 2j+e), c = DD8 refinement mask,
  exact to degree 7 (sum c_k k^m = 2 delta_m0, m <= 7). F = fine nodes in fine
  neighbours, DD8 prediction from exact coarse-lattice samples elsewhere. All
  composed ghost rows exact for tensor degree <= 7 on real meshes (-T: 3e-12 of
  160; row sums exactly 1).
  Smooth-field velocity error (rho = 0.5, base 128 + one band level; uniform 128:
  3.38e-7):
      injection (-j)              1.45e-3     with -b  5.0e-5 (border), 4.7e-3 max
      smooth restriction          2.04e-3     with -b  1.0e-4 (border), 4.7e-3 max
  Root cause: node weights V oscillate on the node scale next to a level border
  (V/h^2 = 1.69, 1.09, 0.86, 0.98 on the fine side, 1.10, 0.96, 1.01 coarse
  side). Two sources: (a) block-wise high-order integration weights q (DD8 limit
  function, also the remeshing-kernel integrals) have Gregory-like end weights
  (1/2, 1.065, 0.986, ...) whose complement a coarser neighbour cannot supply;
  (b) any injection in the composed rows (cself, coarse samples of fine data)
  concentrates the transposed weight on single nodes. A smooth partition of
  unity (C^7 smoothsteps, half width 7 h_min) removed (a) only partly
  (V/h^2 in [0.79, 1.69], velocity 1.5e-3) and lowered quadrature accuracy
  (1e-7), so it was removed. -b violates the FMM far-field core assumption
  (4.7e-3 elsewhere). Conclusion: lab + transposed fold conserves exactly but
  the remeshed masses at a border are consistent with an effective volume
  ~ own-particle share, not with any smooth high-order V. A consistent design
  needs a blended overlap zone (coarse and fine particles both present with a
  partition of unity) or lifted (moment-preserving) transfers throughout.

3.2 Fair LTS test (co-rotating Gaussian pair, -P -r 0.15 -e 0.25, n = 256,
  t = 0.5, fast blocks |grad u| > 0.1 max + buffer: 13% of particles), error vs
  same-H Heun M = 96:
      H      Heun M=1  M=2     M=4     M=8     MRI M=3  M=6     M=12    M=24
      0.1    1.2e-2    1.6e-3  2.1e-4  2.6e-5  5.1e-4   8.1e-5  3.5e-5  3.2e-5
      0.05   1.6e-3    2.0e-4  2.5e-5  3.2e-6  6.3e-5   9.7e-6  3.9e-6  3.5e-6
      0.025  2.0e-4    2.5e-5  3.2e-6  4.0e-7  7.9e-6   1.2e-6  4.7e-7  4.1e-7
  MRI order 3 (slow floor 3.2e-5, 3.5e-6, 4.1e-7: ratios 9.2, 8.4). Equal
  error: 3.5e-5: MRI M=12 2.3 s, 6.3e9 pairs vs Heun M=8 5.2 s, 1.9e10 pairs;
  4e-6: 4.6 s / 1.2e10 vs 10.3 s / 3.8e10; 4.4e-7: 3.7 s / 2.5e10 vs 6.3 s /
  7.6e10. Speedup 2-2.3x wall, 3x pairs. Model 1/(f + (1-f)/M) = 4.9 at f = 0.13,
  M = 12; the gap: MRI reaches the slow-part error floor at M ~ 12, and the
  3 slow evaluations per macro step cost a full FMM each.

3.3 Scaled expansions: a~_k = a_k/s^k, b~_l = b_l s^l, s = box size. All
  translation tables then depend only on integer offsets (one set for all
  depths: M2M/L2L 4 each, M2L 49, far field with unit-lattice sums), computed in
  mpmath (gen_vm2.py fmm_tables) and rounded once. -v ranges: |M_k| <= 1,
  |L_k| <= 23 at every depth (before: 1e-29 .. 1e11). FMM vs direct 8e-11.

3.4 Dynamic regridding (-R k: every k macro steps, sync points only):
  detail = |omega - DD8 prediction from even nodes| (lab ghosts), max per block
  over max|omega|, times exp(K lambda tau) (lambda = principal strain rate,
  tau = regrid interval; -K 0 = c1, -K > 0 = predictive c3); refine if > eps_r,
  coarsen (all 4 siblings) if < eps_r/2^8; one buffer block; 2:1 balance;
  transfer old -> new by copy / injection / recursive DD8 prediction composed
  from the old ghost rows (any level jump); circulation restored exactly by a
  global multiple of V. Gate result: FAILS. Adaptive runs (base 64, 3 levels,
  eps_r 1e-3..1e-7) start sparse (pair: 44k nodes = 17% of uniform 512; shear
  layer 115k) but refine to the full uniform mesh within one time unit; errors
  are worse than uniform 512 at the same node count:
      pair (d = 1.0, t = 1): uniform 128/256/512: 5.6e-3, 3.5e-4, 6.5e-6;
          adaptive eps_r 1e-4..1e-7: 8.5e-4 .. 4.0e-4 at 2.4e5 nodes
      shear layer (t = 4): uniform 128/256/512: 2.0e-2, 7.4e-4, 1.4e-5;
          adaptive: 1.8e-3 .. 8.5e-5 at 2.6e5 nodes
  c1 vs c3 (K = 1, 8) give the same meshes once refinement saturates.
  Cause: the node-scale border oscillation of 3.1 appears in omega = Gamma/V
  as wavelet detail, so every border triggers refinement of its neighbours.
  Regridding works mechanically; it becomes useful only after 3.1 is solved.
  (Superseded: with 3.6 the gate passes.)

3.5 CUDA graphs (-Gr; build with -default-stream per-thread): thrust replaced
  by cub radix sort with preallocated storage and deterministic custom kernels
  (iota, leaf starts, sums for the periodic correction, max displacement); one
  macro step captured and replayed, recaptured after each regrid. Results
  bitwise identical. Speed gain ~1% (n = 128..512): the step is compute bound in
  P2P (phase timing VM2TIME=1: P2P 85-95%). P2P with 256 threads per leaf (one
  per target) instead of 64: 2.5x faster at n = 128, 1.2x at n = 512.

3.6 Level borders (~/dns/vm/BORDER.md; 1D study ~/dns/vm/border1.py). Cause of
  3.1: a core width that jumps or varies between sources near a target gives
  O(1) and O(|grad delta|) errors. Node-volume quadratures are exact only for
  polynomials over the DD stencil, not for integrands that vary on the core
  scale. Implemented:
  - velocity: partition of unity of the level populations (C^11 product ramps
    over 16 coarse cells, coarse samples = even fine nodes), one core delta(x)
    per target, extra near list C for leaves within 7 delta_max (core
    correction of the point far field);
  - transport: per-level gather remeshing from own + ghost particles (omega
    and displacement from the ghost tables as values), with a global
    circulation correction;
  - PSE: per-target delta;
  - state is omega;
  - mesh: the fine region keeps nbuf = ceil((2 ramp + 8 + band)/B) block
    widths of level >= l - 1 around level l, and regrid keeps the same buffer
    around every flagged block.
  Smooth field: 3.41e-7 (uniform 3.38e-7; was 2.04e-3). G1 and G2 reproduce.
  Two-level shear layer: 7.4e-8 with 62% of the uniform-512 nodes (uniform
  512: 5.9e-8). The regridding gate (3.4) now PASSES: the pair reaches the
  uniform-512 error 6.47e-6 with 87k nodes (33%), the shear layer 1.43e-5 with
  172k (65%); c3 helps only at eps_r >= 1e-4. Tables in BORDER.md section 5;
  Maxima/sympy cross-checks of the ramps, core and PSE moments and the buffer
  rule in section 6 (~/dns/vm/maxima/).

3.7 Regrid on the GPU (default; -H = host reference, -X = GPU checked against
  the host at every mesh, bitwise). Before, every regrid copied all values to
  the host and rebuilt there: mesh, node weights, ghost rows (CSR), FMM tree
  and lists, then uploaded them. Now the whole regrid runs on the GPU and only
  scalars come back (sizes for allocation, one flag per balance pass).
  Steps, in the order of regrid():
  a. Detail and strain per block: lab fill, detail/strain kernels, block
     maximum (one CUDA block of 256 threads per mesh block), max|omega|.
     d = detail / max|omega| * exp(K lambda tau); flags and coarsening (all
     four siblings, looked up in the fine map) as before.
  b. Targets on the fine lattice (finest block size, nf = n/B 2^L per side):
     own new level, then each flagged block writes its level over nbuf
     new-level block widths (atomicMax).
  c. Balance (buffer rule + 2:1, BORDER.md 5): fixpoint on the fine lattice.
     Per pass: max pyramid of the targets -> leaf level of every cell, one
     thread per leaf checks its window. The host rule is order dependent (leaves
     scanned in row-major order of the corner, a later write replaces an
     earlier one); the GPU reproduces it exactly with a 64-bit atomicMax on
     (2 corner + 1 + own) << 32 | level, then applies the winners. Taking the
     maximum instead (order free) gives different meshes at coarse thresholds
     (pair c1: 30.5k nodes, L2 8.2e-4 at eps_r 1e-3, instead of 28.6k, 1.05e-3),
     so the old rule was kept.
  d. Leaves in Morton order (x bit above y bit): leaf flags in Morton order of
     the fine lattice, exclusive scan (CUB) -> leaf index; fine map cell ->
     leaf (it is also the FMM leaf map). Was level-major row-major; only
     round-off changes.
  e. Neighbour cases: the ghost rows of a block depend only on the relation of
     its 8 neighbours (same, coarse, fine) and, with a coarse neighbour, its
     parity. With 2:1 balance there are exactly 528 cases (enumerated: coarse
     neighbours come in groups that share one coarse block; no fine next to a
     coarse neighbour). Rows are compiled once at startup from the gen_vm2.py
     op tables with symbolic columns (slot k*4+q, own, coarse group), 5.4 s,
     22.1M entries (265 MB) on the GPU, independent of the mesh. Per block a
     kernel computes rel, the case and a slot table (45 block ids); the lab fill
     reads the case rows and maps slots to blocks. No per-mesh CSR exists.
  f. FMM tree: box (d, X, Y) exists iff the leaf covering its corner has depth
     >= d. The host BFS numbering is (depth, Morton) order, so one exclusive
     scan over all depths gives every box index; parent and first child come
     from the same index array.
  g. Lists U, W, V, X, R, C: count pass, scan, fill pass, one thread per leaf
     or box, same traversal and order as the host build (U/W stack walk; X = X0
     of the ancestors, shallow first; R = stable radix sort of the X entries
     by source leaf; C from the fine lattice around the leaf, sorted (leaf,
     image) keys, then the host test).
  h. PoU weights and core widths: one CUDA block per mesh block, candidates
     (coarser blocks within the ramp reach) from the fine lattice, ascending.
     Shared host/device code with correctly rounded operations (__dmul_rn ...;
     host built with -ffp-contract=off) and the smoothstep powers by products
     (no pow), so GPU and host give the same bits.
  i. Transfer old -> new: per new node, the old block covering it (old fine
     map), copy/injection or recursive DD8 prediction from the old lab (own +
     ghost values, the existing ghost machinery); global circulation
     correction.
  Also on the GPU now: the initial condition, the diagnostics (circulation,
  enstrophy, energy, max|omega|: 5 doubles per output) and the pair counter
  (was a copy of the leaf starts every FMM call). -c, -T, -v and -l still use
  the host mesh (they read the host lists).
  Gates (hal A100, hg GH200):
   G-R1 -X: every mesh of the pair and shear-layer regrid runs (e.g. shear
     K 8 eps_r 1e-6, t = 4: 104 meshes; pair K 0 eps_r 1e-3: 29 meshes) (targets,
     leaves, fine map, rel, boxes, U/W/V/X/R/C, every ghost row, weights, core
     widths, transferred values) bitwise equal to the host path, on the A100
     (x86 host) and on the GH200 (aarch64 host).
   G-R2 regrid runs (base 64, -L 3, -R 4, H 0.01, M 2), GPU = host bitwise in
     every output column; vs the old code (hal; shear K 8 old from glados):
       case  K  eps_r  old L2 / nodes             new L2 / nodes
       pair  0  1e-6   6.4675371560603e-6 87040   6.4675371560494e-6 87040
       pair  8  1e-6   6.4675436328434e-6 89283   6.4675436325753e-6 89283
       pair  0  1e-3   1.0458477465866e-3 28580   1.0458477465867e-3 28580
       pair  8  1e-4   1.0239320474968e-3 56228   1.0239320474960e-3 56228
       shear 0 1e-6    1.4321080300379e-5 171628  1.4321080300488e-5 171628
       shear 8 1e-6    1.4331106147481e-5 180191  1.4331106146919e-5 180191
     All 20 (pair, shear) x (K 0, 8) x (eps_r 1e-3..1e-7) runs: the same mean
     node counts, L2 equal to 1e-12 relative (round-off).
   G-R3 regrid time per call (GH200, shear layer with rho = (pi/15) 64/n,
     -L 3, -E 1e-6, -K 0, regrid every step; mean of the in-run calls):
       n     nodes    GPU ms   host ms   ratio
       64    1.1e5    4.6      488       106
       256   5.6e5    11.9     4635      390
       512   1.25e6   16.9     9769      578
       2048  8.1e6    57       89473     1570
     The host time grows like (blocks)^2 (PoU candidates and C lists scan all
     blocks); the GPU time is ~3 ms of fixed cost (about 120 cudaMalloc/Free,
     ~10 scalar reads) plus ~7 ns per node. With -K > 0 one extra velocity
     evaluation (FMM) is part of every regrid in both paths (pair at 9e4
     nodes: 43 ms GPU, 679 ms host per call).
   G-R4 nsys (hg, short pair run, profiler range = time loop): no
     host-to-device copies, 78 device-to-host copies of 4, 8, 28 or 40 bytes
     (list sizes, balance flag, error flag, per-depth box offsets, displacement,
     pair count, diagnostics). Old code on the same run: 212 device-to-host
     copies, 33 MB (189 above 1 KB), and 241 MB host-to-device.

3.8 Mapping to 3D (vm3)
  Data structures carry over with one more coordinate: Morton codes with 3
  interleaved bits (x above y above z, matching the octree child order), the
  fine lattice of finest blocks (nf^3 cells), box (d, X, Y, Z) exists iff the
  leaf covering its corner has depth >= d, BFS index = (depth, Morton) rank
  from one scan. Lists: U = 27 neighbours, V up to 6^3 - 3^3 = 189, X/W as in
  2D, R = transposed X, C from the fine lattice; the count/scan/fill kernels
  are the same with larger per-thread buffers (or one CUDA block per leaf).
  What does not carry over:
  - Neighbour cases. 26 neighbours; with parity the coarse slots form 7
    groups of 1, 2, 2, 2, 4, 4, 4 slots that share a coarse block, 7 slots
    sit in the parent. Balanced cases: 73,084,160 (2D: 528). Precompiled rows
    are out of the question. 3D lab fill executes the gen_vm2-style section
    programs on values instead of composed rows: per block, per neighbour
    slot, the section for its relation (same: copy; fine: injection or the
    smooth restriction; coarse: fill the coarse buffer from cself/csame/cfine/
    ccoarse, then the DD8 prediction passes in x, y, z). Sections depend only
    on (slot, relation, parity): about 26 x 5 + 26 x 8 + 1 = 339 of them. The
    per-block state is 26 relations (2 bits each) + 3 parity bits + 27
    neighbour ids. Composed rows were needed only for transposed (mass)
    folds; per-level transport and the PoU use values, so nothing needs them.
  - Balance window cost. Per cell, a leaf of size s checks (s + 2 nbuf s)^3 =
    729 s^3 cells; checked per level-l block instead (the window is a union
    of level-l blocks and a block is either inside one coarser leaf or
    refined), it is 729 lookups per leaf. Same result.
  - Transfer. DD8 prediction is a tensor product of 8^3 = 512 weights per
    refined level; do it as three 1D passes per block (24 per node and level)
    instead of the 2D recursion.
  - Fine lattice size. nf^3 ints: n/B = 16, L = 4 gives 256^3 = 1.7e7 cells
    (67 MB per int array, ~0.3 GB for targets, levels, write keys, fine map).
    Beyond nf ~ 512 use the sparse form: leaves sorted by Morton key, the leaf
    containing a cell = last key <= its code (binary search); balance on the
    leaf list with the per-block window above.
  Measured 2D phases (GH200, 8.1e6 nodes, 31.7k blocks, per regrid, nsys):
  PoU weights and core widths 21.8 ms, transfer 3.9 ms, circulation sum
  (one-block sums3) 3.8 ms, max|omega| (one block) 3.4 ms, C lists 1.3 ms,
  lab fill 0.8 ms per field, detail 0.65 ms, balance 0.3 ms (5 passes over a
  1M-cell lattice), tree + U/V/W/X/R lists + leaves + cases < 0.5 ms; the rest
  of the 58 ms is allocation and the scalar reads.
  Expected 3D cost, 1e8 particles (B = 16: 24k blocks of 4096), scaling the
  2D phases by nodes and by the 3D work per node: PoU ~0.3-0.4 s (three ramp
  products per term; the smoothstep is evaluated as a Bernstein sum with
  explicit powers, ~600 multiplies; a Horner form would cut this ~10x, it
  changes only round-off), lab fill ~70 ms per field (7 ghosts per node
  instead of 3, longer rows), transfer ~50 ms (separable prediction), detail
  ~10 ms, global sums ~50 ms (make them multi-block), balance, tree and lists
  < 20 ms, allocations ~20 ms (a pool would remove them). Total ~0.5 s per
  regrid, plus one velocity evaluation if -K > 0. One velocity evaluation at
  1e8 in 3D is minutes of P2P (27 x 4096 sources per particle), so the GPU
  regrid is well below 1% of a step. The host path would be minutes per
  regrid (2D at 1.25e6 nodes already 9.8 s, growing like blocks^2) plus
  copying 1e8 values each way. Peak memory: old and new mesh coexist during
  the transfer (state + lab of both).

4. Remaining open problems
--------------------------
- Level borders: solved (3.6). The cost is the ramp: fine regions must be
  >= 2 x 56 fine cells wide, so very thin refined features pay a wide buffer.
- LTS (fastpair) uses the new per-target cores but was not re-gated with
  per-level transport.
- Pair at d = 0.5 filaments quickly: sp2 512 vs 1024 differs 6e-3 at t = 1, so
  any 512 run is under-resolved there.
- P2P throughput 3.8e10 pairs/s (FP64): next optimisation target.
