Level borders in remeshed particle methods: 1D study and 2D check
==================================================================

Code: ~/dns/vm/border1.py (1D, numpy, explicit operator matrices), results
out/border1.txt; 2D checks with ~/dns/gpu/vm2.cu (-c, -U). Setup: periodic
[0, 2 pi), fine lattice h = 2 pi/Nf, fine region [L/4, 3L/4) at h, rest at 2h,
nested node-centred lattices, m8 remeshing kernel (order 6 when remeshing every
step), DD8 interpolating prediction. "Velocity" in 1D = periodic Hilbert
transform (the 1D Biot-Savart analogue) with an order-8 core:
K_delta(r) = cot(r/2)/(2 pi) (1 - exp(-x) P(x)), x = r^2/delta^2, P degree 4,
P(0) = 1, int exp(-rho^2) P rho^(2k) drho = 0, k = 0..3 (error O(delta^8)).

1. What has to hold
-------------------
A particle carries a mass m_j (circulation). Two different things use it:
(a) the velocity u(x) = sum_j K_delta(x - x_j) m_j, a quadrature of
    int K_delta(x - y) omega(y) dy with weights m_j = omega(x_j) V_j;
(b) the field omega_j = m_j / V_j (diagnostics, wavelet details, regridding).
For (a) to be high order three conditions are needed:
 C1 sum_j V_j g(x_j) = int g + small for every g that varies on the scale of
    the core delta (not only for polynomials over the stencil width);
 C2 the core width seen by a target must be the same for all sources near
    it: K_delta - K has vanishing moments over the whole line; if delta jumps at
    a border the moments of the truncated difference do not vanish and the
    error is O(1) within a few delta of the border, for any h;
 C3 delta must not vary from source to source near a target: with a smooth
    delta(y) the error has a term ~ |d delta / dy| omega that does not decrease
    with h.

2. What the current vm2 scheme does (S0 in border1.py)
------------------------------------------------------
Lab deposition plus transposed fold: fine particles deposit on the fine
lattice, the masses on fine cells beyond the border go to coarse nodes through
the transpose of the DD8 prediction; coarse particles deposit on the coarse
lattice, masses on coarse cells inside the fine region go to fine nodes through
the transpose of injection (or the smooth restriction). This conserves sum m
exactly (rows of every ghost operator sum to 1). Node volumes are V = R^T(h)
(R = reconstruction of the fine-lattice field from the nodes).
Why V oscillates at the border (the 1.69, 1.09, 0.86, 0.98 of 2D): R contains
injections (coarse samples of fine nodes) inside the DD8 stencils of predicted
fine points near the border; R^T puts their weight on single (even) fine
nodes. With the q-weights (integrals of the DD8 limit function) the Gregory-like
end weights q = 1/2, 1.065, 0.986, ... of a block are not complemented by a
coarser neighbour. Measured V/h in [0.75, 1.5] (1D).
Consequences (1D, max error, Nf = 128 .. 1024):
                      Hilbert velocity        transport (pulse once around)
   uniform fine       1.5e-8 .. 2e-13         6.3e-4 .. 2.5e-9 (order 6)
   S0                 0.096 .. 0.101          0.094 .. 0.050  (no convergence)
   lifted N.Nt=8.2    0.098 .. 0.101          0.11  .. 0.058
   lifted 8.4         0.098 .. 0.100          0.11  .. 0.059
   MURPHY substitution 0.096 .. 0.101         0.051 .. 0.022
   blend (PoU, pair delta) 6.8e-3 .. 3.4e-3   2.1e-2 .. 5.8e-4
   per-level (Bergdorf) 0.080 (nodal volumes) 2.2e-2 .. 1.3e-7 (order ~6)
With one common core width (delta = 4h everywhere) S0 still gives 2.8e-3 (C1
fails: the DD8-based V is exact only for polynomials over ~16 fine cells, while
K_delta varies on 4 cells) and the blend 2.0e-4.
Spurious wavelet detail next to the border after transport (Nf = 512):
S0 4.4e-2, lifted 3.3e-2 / 3.0e-2, MURPHY 1.6e-2, blend 1.5e-4, per-level 2.5e-9
(= far-field level 3.5e-9). Only per-level transport is free of spurious detail;
with the others every border triggers refinement (this is what made the 2D
regridding gate refine the whole box).
Locally fine with one-sided extrapolation instead of injection: unstable
(V < 0, transport blows up); dropped.

3. What works
-------------
Velocity: partition of unity of the level populations + per-target core.
 - fine population: fine nodes, weight phi(x) h; coarse population: coarse
   nodes outside, plus coarse samples under the fine region (the even fine
   nodes), weight (1 - phi) 2h; phi = C^k smoothstep from 0 at the border to 1
   inside, over w coarse cells. Each population is a uniform-lattice trapezoid
   of a smooth integrand (C1 holds).
 - every target uses one core width delta(x) for all sources (C2, C3):
   delta = 2 h_c within the ramp (+8 cells), easing to 2 h_f deeper inside.
 - fine particles at even nodes in the ramp simply carry weight
   phi + 4 (1 - phi), odd ones phi: no extra particles are needed.
 1D, max Hilbert-velocity error (relative to O(1) field):
   C^7 ramp, w = 8: 7e-6; w = 16: 2-3e-8 (Nf 256..2048); C^11, w = 16: 5-7e-9;
   C^inf (exp) ramps are worse (steep middle); C^23 worse (same reason).
   The floor is set by w in cells (ramp spectrum at the coarse Nyquist), not by
   h; near the border it equals uniform-coarse accuracy (uniform coarse at the
   same delta: 1.5e-8 at Nf 256), which is inherent: the border sits where the
   coarse level resolves the field.
 2D (vm2 -c -U, smooth field rho = 0.5, base 128 + one band level, velocity at
 all nodes vs the spectral velocity): 2.0e-3 (current) -> 3.9e-7 (PoU, ramp
 12-16 coarse cells, C^11); uniform 128: 3.4e-7. Target met.
 Evaluation note: u = u_FMM(original) + D_chi(PoU) - D_chi(original), D_chi =
 direct sums with a smooth radial window (r < 1.0, C^7 to 0 at 1.5). Plain 3x3
 image truncation is wrong for this: it depends on second moments of the source
 set, which differ between quadrature rules through non-periodic moments at the
 domain seam (1.6e-4 / 1.2e-3 relative) and produced a spurious strain 1e-4.

Transport: per-level particles with buffer bands (candidate e).
 Each level keeps its own particles on its own lattice over its active region
 plus a buffer of (kernel half support + displacement) cells; buffer values come
 from the coarser level (DD8 prediction); coarse values under the fine region
 are replaced by injection of the fine values after each step; every level
 remeshes only same-level particles. No mixed volumes ever enter a remesh.
 1D: order ~6 (2.2e-2, 5.0e-4, 8.1e-6, 1.3e-7 for Nf 128..1024), no spurious
 detail (2.5e-9 vs far 3.5e-9 at 512). Circulation is not exactly conserved:
 relative drift 1.3e-5, 5.3e-9, 1.6e-11, 7.2e-13 (high order, shrinks with h);
 a global correction (constant times V) restores it to round-off, local first
 moment is not exactly conserved (high-order error).
 Trade-off: exact circulation conservation (S0 and transposes) costs O(1)
 errors at borders; per-level transport gives up exactness (high-order drift,
 globally corrected) and gains high-order accuracy and clean wavelet details.
 It composes with dimension-split remeshing: every level is a uniform lattice
 inside each block plus ghosts at the same level, so 1D passes across faces are
 ordinary same-level passes.

4. Recommendation for 2D/3D (vm2, vm3)
--------------------------------------
- Transport: per-level gather remeshing (each block: own nodes + lab ghost
  particles at its level, ghost values from the ghost tables as values, not
  masses); ghost particles advected with the same velocity field (extra FMM
  targets in a band of kernel+displacement cells, or velocity interpolated from
  the lab); global circulation correction.
- Velocity sources: PoU weights at level borders (C^11 ramp over 12-16 coarse
  cells inside the finer block), per-target core delta(x) in P2P; the FMM near
  field of leaves with wide-core targets extended to the 5x5 ring (offset-2
  boxes moved from V-lists to direct), so the far-field point expansion keeps
  rho >= 8.
- PSE: per-target delta and target-side quadrature (not conservative pairwise;
  conservation comes from the global correction like transport).
- Regridding criterion on node values (no division by mixed volumes).
- Not adopted: transposed folds for transport (O(1) border error, spurious
  detail), lifted/MURPHY variants (do not fix C1-C3), one-sided extrapolation
  (unstable).

5. Implementation in vm2 (time stepping)
----------------------------------------
State: node values omega (no longer Gamma). Per node a PoU weight P and a core
width delta, computed after each mesh build (on the GPU, kpou; host
reference pou() with -H; ~/dns/gpu/vm2.md 3.7):
  phi_k(x) = prod over blocks c of level < k of
             1 - (1 - S23(gx_c/wr_k)) (1 - S23(gy_c/wr_k)),
  gx_c, gy_c = periodic gaps to block c, wr_k = ramp h_(k-1), S23 = C^11
  smoothstep (degree 23). The product form is C^11 at block corners (a plain
  L-infinity distance has kinks there).
  P = sum over k <= l, node on lattice k, of h_k^2 phi_k (1 - phi_(k+1)).
  delta = cd (h_0 + sum_k (h_k - h_(k-1)) s_k), s_k = same product with a C^7
  ramp starting wr_k + 8 h_k from the coarse block, over band h_k.
Velocity: FMM sources Gamma = omega P, every target uses its own delta in P2P.
  Extra near list C (per leaf): leaves within 7 delta_max that are not
  in U or X get the correction -(exp(-x) R(x)) K_point, i.e. the core that the
  point-kernel far field (M2L, M2P) leaves out; it is empty in uniform
  interiors. PSE uses the same per-target delta:
  d omega_i/dt = nu/delta_i^4 sum_j P_j (omega_j - omega_i) eta(r/delta_i).
Transport: each block gathers from its lab: own particles plus ghost
  particles at the block's level, with omega and displacement taken from the
  ghost tables (as values: same-level copy, DD8 prediction from coarser,
  smooth restriction from finer). One kernel (remesh), normalised m8 weights.
  The global correction omega += (Gamma_0 - sum omega P)/sum P runs every
  step.
Mesh rule: the fine region must contain the ramp and the delta band. A
  level-l block keeps all blocks within nbuf = ceil((2 ramp + 8 + band)/B)
  level-l block widths at level >= l - 1 (build). Regrid puts the same buffer
  around every flagged block, including those already at the maximum level.
  Without it the vortex cores sat in the delta band (delta = 1.5 h_f) and the
  adaptive error stalled at 4.6e-5 for every eps_r (pair).
Flags: -a ramp (coarse cells, default 16), -A band (fine cells, default 16).
  -U, -b and the env vars VM2WHERE/VM2SAME/VM2RAMP/VM2DB are removed. The
  helper cuda() was renamed ck() (the CUDA 13 headers declare a namespace
  cuda).
GPU regrid (vm2.md 3.7): the same rules run on the GPU. The buffer/2:1
  fixpoint keeps the host's write order (row-major leaf scan, last write
  wins) through a 64-bit atomicMax key; an order-free maximum gives
  different meshes at coarse eps_r. Smoothstep powers are products (no pow)
  so host and GPU weights agree bitwise; this changed P and delta at the
  round-off level only.

Gates (glados A100):
 smooth field velocity (rho = 0.5, -c), FMM path:
   uniform 128 3.38e-7; base 128 + one band level 3.41e-7 (old 2.04e-3).
   FMM vs direct with per-target cores 1.9e-11.
 G1 uniform unchanged: 128 1.7513e-4, 256 3.5463e-6 (old 1.75e-4, 3.55e-6).
 G2 base 256 fully refined == uniform 512: all columns bitwise equal.
 two-level shear layer (t = 1, nu 1e-4, H 0.005, M 2; ref sp2 1024):
   uniform 256 3.55e-6, uniform 512 5.93e-8;
   base 256 + band level (164k nodes, 62%): 7.44e-8, Linf in border blocks 1.3e-9.
   base 128 + band (41k): 8.2e-6. The fine band is 64 fine cells wide,
   narrower than the ramp + delta zone (2 x 56), so the cores see delta > h_f.
   Circulation stays at 1e-13 (round-off).
 regridding (base 64, -L 3, -R 4, H 0.01, M 2; mean nodes; c1 = -K 0,
 c3 = -K 8):
   pair (-P -r 0.15 -e 0.5, t = 1): uniform 128/256/512 5.57e-3, 3.46e-4,
   6.47e-6 at 16k/66k/262k nodes
     eps_r   c1 L2     nodes   c3 L2     nodes
     1e-3    1.05e-3   28.6k   9.93e-4   34.3k
     1e-4    1.09e-3   50.2k   1.02e-3   56.2k
     1e-5    1.06e-5   65.0k   7.46e-6   68.8k
     1e-6    6.47e-6   87.0k   6.47e-6   89.3k
     1e-7    6.47e-6   96.4k   6.47e-6   103.7k
   shear layer (t = 4): uniform 128/256/512 1.97e-2, 7.38e-4, 1.43e-5
     1e-3    5.62e-4   57.0k   4.87e-4   66.1k
     1e-4    1.56e-4   98.6k   9.33e-5   122.0k
     1e-5    1.46e-5   143.3k  1.44e-5   151.1k
     1e-6    1.43e-5   171.6k  1.43e-5   180.2k
     1e-7    1.43e-5   224.9k  1.43e-5   226.4k
   The adaptive runs reach the uniform-512 error with 33% (pair) and 65%
   (shear) of its nodes. At the uniform-256 node count the error is 33x
   (pair) and 5x (shear) below uniform 256. The old scheme refined the whole
   box and stalled at 4e-4 to 2e-3. Errors are independent of the regrid
   interval (-R 1, 2, 8, 25: 4.6e-5 in the stalled case, so the transfer adds
   nothing). c3 (predictive) beats c1 only at coarse thresholds (eps_r 1e-4:
   1.7x lower error for 1.24x the nodes) and is equal below. Level-border
   blocks are no longer where the error is (pair eps_r 1e-6: 1.6e-11 there,
   1.4e-5 in level 3).

6. Cross-checks (sympy and Maxima, ~/dns/vm/maxima/)
-----------------------------------------------------
border.mac (maxima --very-quiet < border.mac) and border_sympy.py give the
same results:
 - ramps S15, S23: S(0) = 0, S(1) = 1, derivatives 1..7 / 1..11 vanish at
   both ends, the next ones do not (259459200, 647647525324800),
   S(t) + S(1 - t) = 1.
 - order-8 core q = 1 - exp(-x) R(x), x = r^2/delta^2: R = 1 - 3x + 3/2 x^2
   - x^3/6 (the gen_vm2.py linear system gives the same coefficients);
   moments int x^k dq = 1, 0, 0, 0, -24 for k = 0..4 (order 8).
 - PSE order 8: E = (40 - 40x + 10x^2 - 2/3 x^3)/pi; pi int exp(-x) E x^k
   = 4, 0, 0, 0, -480 for k = 1..5.
 - buffer rule: delta zone = 2 ramp + 8 + band fine cells (56 for 16/16),
   coarse-block reach ramp + 4 + band/2 coarse cells (the candidate cut-off in
   pou()), nbuf = 4 for (16, 16, B = 16).
