YAO40K_PLAN: Yao & Hussain 2020a (JFM 883 A51) at Re 40000 -- run plan (nothing queued until Sergey says go)
========================================================================================================
Spectral agent, 2026-10-31. Option A: spectral boxsymmp (this agent). Option B: vm3 (vm3-yao agent, section B, merged from
~/dns/vm/YAO40K_VM3.md). (m) = measured, (e) = estimate or model. Details behind every number: ~/dns/mpi/YAO40K.md.

HEADLINE
  Recommendation: A, spectral boxsymmp on a 3072^3-equivalent grid (quarter box 1537 x 3072 x 1537), 16 H200 (4 nodes).
  Cost A (to t* 3, C 2): ~21k steps; ~8 h wall on 16 H200 (10 h with +30% steps); ~125-160 H200-h (e: the 16-H200
  s/step ~1.5 is a model; one-node speeds are measured; the multi-node H200 timing is still queued). Go/no-go first:
  G1 L3 Re 20000 (~20 H200-h) and G2 a 2048/2560/3072 ladder at Re 40000 near the peak (~70 H200-h). Total with the
  checks ~220-250 H200-h.
  Cost B (vm3, from its agent): ~28k steps, 20-30 h on 16 A100, 320-470 A100-h (e; scaling past one node unmeasured).
  Earliest start: the checks can start as soon as a 1-node H200 slot is free (hours). A 4-node H200 allocation for ~12 h:
  scheduler estimates today put 4 nodes at ~21:00 for a 30-min job; a 12-24 h job realistically 2026-11-01/02.
  Why A: it is the only option whose every ingredient is validated here against the paper (Re 5000 and 9000), whose
  single-node speed is measured on H200, whose memory is fixed and known (no node-count growth), and whose resolution
  criterion comes from our own ladders. B is attractive for its memory-adaptive mesh but its multi-node scaling and the
  |w|max resolution at 2pi/4096 are open (its own estimate: |w|max 5-15% low).

A. Spectral boxsymmp
====================
A1. Code
  ~/dns-yspec, branch yao-spec (from main 06632cf); commits up to c2f6d5e (benchmarks: binary b3 = commit 3a6121e).
  mpi/boxsymmp.cu: MPI + NCCL slab pseudo-spectral, RK4 (5-stage low storage), Hou-Li filter, quarter box by the two
  mirror symmetries of the case: u(-x,y,z) = (-u,v,w) (collision plane x = 0), u(x,y,-z) = (u,v,-w) (symmetry plane z = 0);
  x = collision normal, y = propagation, z = tube axis. Stored 1537 x 3072 x 1537 for 3072^3 (= Yao's 4096 x 2049^2 layout
  at 3/4 the size). Options used: -I yao (IC on the GPU), -c C (adaptive dt), -F (filter as a rate), -S (float all-to-all),
  -k / -r k.STEP -G (spectral checkpoints, padded restart), -e (spectra), -L (lambda2 float crops). Driver scripts on FASRC
  netscratch yspec/: seg.sh (one resumable segment), prod40k.sh (phases with padding).

A2. Grid: 3072 (2560 minimum), not 4096
  Criterion (m, our ladders): keta = kc eta, kc = (n-1)/3, box-mean eta at the enstrophy peak.
    Re 5000:  384 keta 2.96 -> Omega -0.3%, |w| peak -6%;  512 keta 3.96 -> Omega 0.07%, |w| peak -1.8% (vs 768)
    Re 9000:  768 keta 4.02 -> Omega peak -0.15%, |w| peak -1.4% (vs 1024)
  keta ~ 4 gives Omega to 0.1-0.2% and the |w| peak to ~2% at both Re. At Re 40000 (m, r40k1024): keta 1.77 at 1024 with
  its own (too low) enstrophy, 1.69 with the paper's peak Omega/Omega0 11.8 -> keta 4 at n = 2420, 5 at n = 3030.
  2560: keta 4.2 -- the minimum. 3072: keta 5.1 (4.8 if the true peak enstrophy is 15% above the paper's) -- chosen, the
  margin covers the unknown: the Re 40000 peak is more intermittent than at 5000/9000, so the box-mean criterion may be
  optimistic. Yao's 4096 with 2/3 truncation is keta 6.8 in this measure; our filter keeps ~0.45 n (vs n/3), so the
  kc = n/3 measure is conservative for us.
  1024 at Re 40000 (m): Omega peak 9.74 (paper 11.81), |w| ratio 124 (paper 569): far from enough, already 7% off in |w|
  at t* 1.0.

A3. IC and its scaling
  -I yao builds yao.py's field on the GPU at any grid (nearest-point Newton, then u = i k x omega / k^2): no IC file.
  Checked against yao.py to all printed digits at 64 and 128 (E, Omega, max|w|); Gz/G0 = 1 - 1e-6 at 128.
  Validated with this IC: Re 5000 at 512 (Omega peak 1.5400 at t* 1.790; paper 1.5401 at 1.785; Omega within 0.3% to
  t* 2.2) and Re 9000 at 1024 (2.271 at 1.835; paper 2.288 at 1.829; within 1% to t* 2.2). Re 40000 changes only nu
  (Gamma0 / 40000 = 1.920096e-4); the IC is the same field.

A4. Phasing
  1536 from t = 0 to t* 1.0 (t 2.147) on 4 H200 (one node), k. dump; zero-padded to 3072 at t* 1.0 and run to t* 3
  (t 6.441) on 16 H200. Checked (m, Re 5000): 512 to t* 1.4 padded to 768 vs 768 from t = 0: Omega 1.8e-6, Gz 5e-7.
  The padding point is conservative: at Re 40000, r40k1024 is still within 0.4% of the paper's Omega at t* 1.0; 1536 has
  1.5x its keta. The ladder (G2) re-checks this.

A5. Time step
  -c 2 -F F with F = 0.00125 x 512 / n (one full filter application per initial step: 2.1e-4 at 3072, 4.2e-4 at 1536).
  (m, Re 5000) against C 0.5: C 1 and C 2 move Omega by 1.2e-5 / 4.4e-5 and Gz by 6e-6 / 2e-5; the grid max |w| in the
  reconnection window (t* 1.5-2) moves 8.5% pointwise for both (spike timing), its peak value ~2%. C 2 = Courant ~1 in
  the h-based measure (C 1 = 0.48), the same footing as vm3's -dqe 1.0.
  Steps (velocity history U(t) of r40k1024; U = max(|u|+|v|+|w|) 4.8 early, 12.5 at t* 1.74, 9-11 after):
    1536 t* 0-1.0: 2680 (C 2);  3072 t* 1.0-3.0: 17990 (C 2) (35975 at C 1). dt 3072, C 2: 3.9e-4 early, 1.4e-4 at peak.
  U at 3072 near the peak is probably higher than at 1024 (sharper sheets): +10-30% steps is possible.

A6. Speed, wall time, GPU-h
  One node, 4 H200 (m, b3, -S): 1024 0.156, 1536 0.633, 2048 1.236 s/step (old code 0.203 / 0.781 / 1.620).
  Inside a node the step is compute-bound; across nodes it is network-bound. A100 4 nodes (one 100 Gb/s port per node, m):
  2048 / 16 ranks 5.79 s/step -S vs 13.08 old; step = exchanges alone + 5%: the overlap works.
  H200 multi-node: not yet measured: yspec-q2 (2 nodes, gpu_requeue, 40 min) and yspec-b2/b4/b8 are queued with
  --exclude=holygpu8a12204; they write to netscratch yspec/bench.txt.
  Model for 16 H200 (4x NDR400 per node): step = max(compute, exchange) x 1.05-1.1, exchange = off-node float bytes at
  ~40 GB/s per GPU: 3072: compute ~1.06 s, exchange ~1.36 s -> ~1.5 s/step.
    phase              GPUs   steps    s/step          wall      GPU-h   (C 2)
    1536, t* 0-1       4      2680     0.633 (m)       0.5 h     2
    3072, t* 1-3       16     17990    ~1.5 (e)   7.5 h   120
    total                                              8 h   ~125
  At C 1: twice the steps. 2560 instead of 3072: ~15000 steps at ~0.95 s (e): ~4 h, ~65 GPU-h.

A7. Memory, checkpoints, I/O
  Per rank (m at 2048/4: 79 GB, scaled by the slab size): 3072 on 16 H200 ~67 GB, on 8 ~133 GB (too tight); 2560 on 16
  ~39 GB; 1536 on 4 ~36 GB (m). Fixed in time (no growth after the collision).
  k. checkpoint = the exact restart (bitwise at the same grid): 174 GB at 3072 (101 GB at 2560, 22 GB at 1536), written
  one pwrite per rank and component (512: 1.6 s for 2 GB; 3072 estimated 2-5 min on netscratch). Every 0.25 time units
  (~12 per run), keeping the newest two (prod40k.sh): 350 GB peak.
  Lab netscratch: 48.8 of 50 TiB used -> ~1.2 TiB free for everyone. Peak use of the run ~600 GB (two checkpoints, the
  kept fields, crops) -- half the lab's free space: clear with the lab before the start or keep one checkpoint (174 GB,
  restart from the previous segment's end if a write is interrupted).

A8. Outputs to keep
  out (every 10 steps): t, E, Omega, C, wmax (grid max |w|), keta, tails, Gz (circulation through z = 0, x < 0, Stokes,
  spectrally exact in y), dt -> Omega(t), |w|max(t), Gamma_z(t) for figs 22a/23a/23c. ~2 MB.
  e.STEP spectra (shell and 1D) every 0.05 in time: ~130 files x 0.4 MB.
  Fields near the peak: k. dumps at t* 1.80, 1.85, 1.88 (174 GB each -> keep 1.85 only if space is short), and lambda2
  float crops (-L, -H 1.2 around the pair centre, ~17 GB each) at t* 1.75-1.95 every 0.05 for rendering.

A9. Risks (A)
  - Resolution criterion extrapolated from Re 5000/9000 (keta 4); the ladder G2 tests it at Re 40000.
  - Multi-node s/step: 16 H200 modelled (8 and 16 H200 pending in the queue).
  - Steps: U near the peak at 3072 may exceed the 1024 value (+10-30% steps).
  - Storage: ~600 GB peak of a lab-shared 1.2 TB margin.
  - Allocation: 4 H200 nodes for ~12 h; gpu_requeue may preempt (k. checkpoints every 0.25 make that a restart, not a loss);
    holygpu8a12204 is broken (always --exclude).
  - |w|max is a pointwise grid max of rolling-up sheets: 5-10% pointwise sensitivity to dt and grid in t* 1.6-2 (m at
    Re 5000/9000); compare peak values and envelopes, not single samples.

B. vm3 (from the vm3-yao agent, ~/dns/vm/YAO40K_VM3.md)
======================================================
B1. Code: branch vm3-yao (~/dns-yao): vm3-mgpu (multi-GPU, -dc 2 -dcc 1, -ac fix) + vm3-symbfs (native quarter bfs) + the
  cases2 outputs. Multi-GPU "mode Q" (-sym 1 with native quarter bfs, distributed charges, PoU weights from each rank's own
  part of the mirrored mesh) is new today; gates vs single GPU: 1 rank bitwise, 2 and 4 ranks 1.5-3e-15 (m). bfs SC 8 on a
  peak-time state: 1.2e-10 vs fs (m).
B2. Grid: n0 128 + 5 levels = finest 2pi/4096 (= Yao's grid) (e). Omega converged < 0.2% at 2pi/256 (Re 3000), 2pi/512
  (Re 5000) (m). |w|max: Re 3000 converged at 2pi/512; Re 5000 not converged at 2pi/1024 (512 -> 1024 +5% / +9% at the two
  peaks) (m); Re 9000 at 2pi/1024: Omega -0.7% vs paper, |w|max +12% (m). Expected at Re 40000: |w|max 5-15% low (e).
B3. Nodes, memory, steps: Re 40000 L 5 on 4 A100 from t 0: 8e6 nodes at t* 0.10, 3.1e7 at 0.48, 6.1e7 at 0.79 (m);
  ~1e8 at t* 1-2, 1.5-3e8 at 2.2-3 (e). 3.2 kB per own node per rank (m) -> at least 16 A100 or 8 H200 by t* 3.
  dt 3.2e-4 at 2pi/4096 (stability, -dqe 1.0) (m); ~28000 steps to t* 3 (e).
B4. Cost: 4 A100, 5.8e7 nodes: 2.5 s/step (m). Model from the tube (16 A100, m): 0.84 s + 1.47e-8 s x nodes ->
  20-30 h on 16 A100, 320-470 A100-h (e); H200 ~1.5-1.8x faster per GPU (e). 32 A100 not better than 16 today (16 ranks
  were 1.2-1.35x faster than 4 on the tube; the fix, rank-aligned patches, is not in vm3 yet).
B5. Risks (B): scaling past one node unmeasured for mode Q (if 1.2-1.35x from 4 to 16: 2-3 days on 16 A100); memory
  tight on A100 at 3e8 nodes; |w|max under-resolved at 2pi/4096; relative refinement threshold (eps 1e-4 x max|w|)
  untested at Re 40000 (at Re 9000 it halved the nodes at the peak).

C. A and B on equal terms
=========================
                          A spectral 3072, 16 H200            B vm3 2pi/4096, 16 A100
  steps to t* 3           ~21k (C 2) (e, U from 1024 run)    ~28k (e)
  s/step                  ~1.5 (e; 8 H200 pending (m))   2.6-3.8 (e; 2.5 at 5.8e7 on 4 A100 (m))
  wall on 16 GPUs         ~8 h                            20-30 h (2-3 days if scaling stays 1.2-1.35x)
  GPU-h                   ~~125 H200-h                         320-470 A100-h (~200-300 H200-h equivalent at 1.5-1.8x)
  memory                  fixed ~67 GB/GPU (m-scaled)         grows to 1.5-3e8 nodes, 16 x 60 GB at 3e8 (e)
  Omega at the peak       ~0.1-0.2% (keta 5 criterion)        < 0.2% (ladder)
  |w|max peak             ~2% (keta 4-5 criterion, ladders)   5-15% low (e)
  validated vs paper      Re 5000, 9000 (to t* 2.6-3.3)       Re 5000 (to 3.7), Re 9000 (to the peak)
  measured                1-node speed (H200), 8-GPU speed,   4-GPU speed, nodes to t* 0.79, memory/node
                          A100 4-node overlap, IC, padding
  modelled                16-GPU speed, steps                 16-GPU speed, nodes after t* 0.8, steps
Spectral is the safer bet on cost and accuracy; B's advantages (memory adaptivity, cheaper per node early) do not pay at
this Re because the post-reconnection cascade fills the box (1.5-3e8 nodes vs 3.6e9 spectral points, but at 3.2 kB per
node vs ~300 B per spectral point the memory is similar), and its scaling past one node is the open item.

G. Go/no-go checks (cheapest first)
===================================
G1. L3: Re 20000 with spectral at 1536 (keta ~5 by the same scaling) and 2048, t* 0 -> 2.6, 4 H200 each: ~1.2 h and ~3.3 h
    (e), ~18 H200-h. Paper: Omega/Omega0 peak 5.51 at t* 1.899, |w| ratio 202 at 1.918. Go if 1536 and 2048 agree to 0.3%
    (Omega) / 3% (|w| peak) and with the paper as at 5000/9000.
G2. Re 40000 ladder near the peak: 1536 to t* 1.4 (3850 steps, 4 H200, ~0.7 h), k. dump; padded to 2048 (4 H200), 2560
    and 3072 (16 H200) for t* 1.4 -> 1.95 (3800 / 4700 / 5700 steps at C 2): ~5 + ~25 + ~45 H200-h (e). Go for 3072 if
    2560 vs 3072 differ < 0.5% in the Omega peak and < 3% in the |w| peak; 2560 is enough if 2048 vs 2560 already do; if
    2560 vs 3072 still differ > 3% in |w|, the plan moves to 4096 on 32 H200 (~900 GPU-h, 29 h).
    This also measures U at the peak on the final grid (the step count) and the 16-GPU s/step.
G3. (B) the vm3 checks: L1 Re 5000 2pi/2048 window on one H200 (|w|max within 3% of 1024?), L2 Re 9000 2pi/2048 window,
    mode Q scaling 4/8/16/32 (go if 16 >= 2.5x over 4), node count at t* 1.0 (y40b).
G4. Before the production start: lab OK for ~600 GB netscratch; the 16-GPU s/step of G2 within 30% of the model.
