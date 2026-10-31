YAO40K: pseudo-spectral (boxsymmp) preparation for Yao & Hussain 2020a (JFM 883 A51) at Re 40000
================================================================================================

Branch yao-spec (worktree ~/dns-yspec, from main 06632cf). FASRC dir netscratch .../slitvinov/yspec (src/, run/, bench/).
Case notes and vm3 numbers: ~/dns/bench/ref/yao/NOTES.txt. Paper curves (fig 23 vector paths): ~/dns/flsr2/yao/.

1. IC, symmetry, orientation (2026-10-31)
-----------------------------------------
yao.py (their section 2): tubes along z, collision plane x = 0, symmetry plane z = 0, propagation along y.
Tube 2 = tube 1 mirrored in x = 0 with the opposite sign:
  omega(-x,y,z) = (omega_x, -omega_y, -omega_z)(x,y,z)  ->  u(-x,y,z) = (-u, v, w)    (boxsymmp Sx)
  centre curves even in s, t(-s) = (-t_x, -t_y, 1)     ->  u(x,y,-z) = (u, v, -w)    (boxsymmp Sz)
so the IC lies in boxsymmp's symmetric subspace with its orientation (x = collision normal, z = tube axis) and the
q. layout (x = i dx >= 0, z = k dz, y = -pi + j dy) is yao.py's grid. Box -X 2 -Y 2 -Z 2, nu = Gamma0/Re.
boxsymmp -I yao builds the IC on the GPU (same nearest-point Newton as yao.py, then u = i k x omega / k^2, which
also projects), so no IC file is needed at any grid. Step 0 against yao.py (full box, numpy):
  n     E (yao.py / -I yao)           Omega                          max|w|                  Gz (Stokes) / Gamma0 profile
  64    0.36631255 / 0.366312553459   2.60260894 / 2.602608943350    20.002722 / 20.0027224   7.680991 / 7.680385
  128   0.36631912 / 0.366319123416   2.60263700 / 2.602637004548    20.000325 / 20.0003247   7.680376
New out columns: Gz = circulation of omega_z through z = 0, x < 0 (Stokes: int v(0,y,0) dy - int v(-pi,y,0) dy,
spectrally exact in y), dt.

2. Code changes (mpi/boxsymmp.cu) and their checks (128^3, Re 5000, 100 steps, A100)
---------------------------------------------------------------------------------
  -I yao       IC on the GPU (above).
  -c C         adaptive dt = C / (Courant rate) at every diagnostic step (10 steps), growth <= 1.1x per 10 steps,
               <= -s; output times (-d -e -L -T, -t) are hit exactly; q.STEP.t / k.STEP.t hold t and dt, restarts
               reproduce the dt sequence.
  -S           all-to-all in float (fold kernels write float2, ext kernels read it).
  synth5       RK stages exchange u, v, w, omega_x, omega_y and form omega_z = d_x v - d_y u in the mixed layout
               after the exchange (x, y spectral there): 5 instead of 6 exchanges per synthesis, 8 instead of 9 per
               stage. BOXSYM_SYN6=1 restores the old path.
  nl_a/nl_b    order u, v, w, omega_y, omega_x; N_x = v omega_z - w omega_y is formed (one more field, b[6]) and its
               R2C, x-stage and exchange start while omega_x is still in flight, so the exchange stream has no gap
               between synthesis and analysis. BOXSYM_NL=0: synth5, cross, analysis (bitwise the same result).
  ring         with -S the two double exchange buffers hold four float buffers: z-stages run up to 4 fields ahead.
  self         the fold kernels write the rank's own section directly into the receive buffer (no local copy;
               bitwise the same; BOXSYM_SELFCOPY=1 restores the copy). 1 GPU: -6% s/step.
  fused diag   the diagnostics every 10 steps use the first RK stage's synthesis (no extra 6 exchanges);
               BOXSYM_NOFUSE=1 restores the own synthesis.
  mask         the filter mask is separable: three 1D tables instead of pow(x, 36) per mode and stage.
  -k, -r k.STEP [-G nx,ny,nz]   spectral checkpoints (exact restart; zero-padded restart on a larger grid).
  BOXSYM_A2AONLY=K, BOXSYM_NOCOMM=1   timing only: exchanges alone / compute alone.
Checks:
  old binary vs new with BOXSYM_SYN6=1: q dump bitwise equal.
  synth5 vs old: q max|diff|/max 3.7e-16, out columns <= 2e-15 (tail columns are rounding-level quantities).
  np 2 vs np 1: q bitwise equal.  -S np 2: repeatable bitwise; NCCL and MPI paths bitwise equal.
  -S vs double after 100 steps: E 1.8e-11, Omega 4.4e-11, wmax 1.9e-9, Gz 8.6e-11 relative. The z tail
  fraction (energy beyond 0.9 kmax) rises from 7e-22 to 1.8e-19: the float rounding floor of the spectrum.
  adaptive restart (q.): same dt sequence and times; fields equal to rounding (the q. read re-analyses).
  k. restart: final k. dump bitwise equal to the continuous run. Padded 64 -> 128 restart runs.
  fused diagnostics vs own synthesis: q bitwise equal, wmax 3.5e-16.  self vs copy (np 1, double and -S): bitwise.
  nl_a/nl_b vs synth5 + cross + analysis: bitwise.  new (all of the above) vs old binary: q 3.7e-16, out <= 2.6e-16.
Kernel time shares (nsys, 256^3, 1 rank, A100 MIG 3g.20gb, 10 steps; all kernels are memory bound and run at
~85-95% of the slice bandwidth): z2z FFTs 29%, ext/fold (parity extension + transpose) 45%, RK update 6%,
C2R/R2C 8%, spread5/cross/vortm 11%. The x and z FFTs run on parity-extended lines (2x the needed length);
removing that needs real DCT-I/DST-I transforms, which cuFFT does not have (not done).

3. Validation, Re 5000 (A100-40GB, 1 GPU, 512^3 = their 512 x 257^2, Hou-Li filter applied once per step, -c adaptive)
---------------------------------------------------------------------------------------------------------------
Paper (fig 23 vector paths): Omega/Omega0 peak 1.5401 at t* 1.785; |w|inf ratio peak 15.28 at t* 1.700.
vm3 (NOTES.txt): 1.539 at 1.789; 14.28 at 1.705.
  run          C     Omega/Omega0 peak    |w| ratio peak     Gz/G0 = 0.5 / 0.2     steps to t 7, s/step
  512          1     1.5400 at 1.7899     15.08 at 1.773     t* 1.8184 / 2.1980    6362, 0.49 (A100-40, old build)
  512 -S       1     1.5400 at 1.7899     15.08 at 1.773     1.8184 / -
  512          1.5   1.5400 at 1.7897     15.21 at 1.770     1.8184 / 2.1981
  512          2     1.5401 at 1.7895     15.34 at 1.770     1.8184 / 2.1981
  384          1     1.5395 at 1.7866     14.43 at 1.712     1.8176 / 2.2014
Ours against the paper curves at t* 1.0/1.4/1.5/1.6/1.7/1.8/1.9/2.0/2.2/2.5/3.0:
  Omega/Omega0  0.9412/0.980/1.034/1.157/1.406/1.538/1.423/1.287/1.129/0.967/0.876
  paper         0.9392/0.977/1.031/1.157/1.414/1.536/1.418/1.285/1.126/0.957/0.869
  |w| ratio     1.085/1.894/3.148/5.34/14.35/10.35/7.587/6.697/5.501/3.512/2.192
  paper         1.086/1.877/3.142/5.39/15.12/10.39/7.594/7.29/5.342/3.427/2.129
  (Omega within 0.3% to t* 2.2 and 1% after; |w| ratio within ~1% except on the spikes at t* 1.7 and 2.0, where the
  trace oscillates within a few hundredths of t*.) Gz/G0 lies on the colour-traced fig 22a points (yao5k.png).
fp32 exchange (-S) vs double, whole run: Omega <= 2.5e-9, wmax <= 1.2e-7, Gz <= 7e-9, E <= 3e-11 relative.
Time step: Omega changes by 5e-5 (C 1.5) and 8e-5 (C 2) against C 1, Gz by 3-4e-5, but the grid max |w| by up to
5-6% in t* 1.5-2.0: without -F the filter is applied once per step, so its strength per unit time follows dt. The dt
study is being repeated with -F 0.00125 (filter as a rate). Resolution: 384 -> 512 changes Omega by 0.2% and the
grid max |w| by up to 19% (t* 1.5-2.0): max |w| is the quantity that sets the grid.

4. Validation, Re 9000 (4 H200, one node, -S, C 1, no -F)
---------------------------------------------------------
Paper: Omega/Omega0 peak 2.2876 at t* 1.829; |w|inf ratio peak 33.98 at t* 1.838 (fig 23 vector paths).
  run     Omega/Omega0 peak    |w| ratio peak     Gz/G0 = 0.5 / 0.2     steps to t 5.6 (t* 2.61), s/step (4 H200)
  1024    2.2707 at 1.8348     37.87 at 1.842     1.8306 / 2.2338       11524, 0.172
  768     2.2672 at 1.8337     37.33 at 1.837     1.8298 / 2.2287       8584, 0.089
Against the paper at t* 1.0/1.4/1.5/1.6/1.7/1.8/1.9/2.0/2.2/2.5 (1024):
  Omega/Omega0  0.9953/1.063/1.151/1.413/1.987/2.188/2.263/2.015/1.476/1.215
  paper         0.9938/1.061/1.150/1.420/2.004/2.211/2.279/2.022/1.486/1.255
  |w| ratio     1.113/2.143/4.341/8.49/22.9/21.0/23.7/23.8/13.0/8.49
  paper         1.113/2.135/4.361/8.61/21.6/19.2/24.1/18.5/13.5/10.7
Omega peak 0.7% below the paper, peak time +0.006; Omega within 1% to t* 2.2, 3% low at t* 2.5 (768 is lower still,
so this is on our side of the resolution trend). Gz/G0 on the fig 22a trace (yao9k.png). The |w| ratio follows the
paper's spiky trace to t* 2.2; its largest spike is 37.9 here and 34.0 there (a single-sample event).
768 vs 1024: Omega 0.4% (t* 1.5-2), 1.1% after; Gz 0.2%; peak values 0.15% (Omega), 1.4% (|w|).

5. Resolution and time step (Re 5000, 512 unless noted)
------------------------------------------------------
Pointwise max differences relative to the reference's max (t* <1.5 / 1.5-2 / >2):
  384 vs 768       Omega 1e-8 / 2.5e-3 / 3.3e-3     wmax 2e-3 / 0.19 / 0.06     Gz 1e-5 / 2.9e-3 / 1e-3
  512 vs 768       Omega 7e-10 / 7.2e-4 / 4.4e-4    wmax 1e-3 / 0.085 / 0.02    Gz 3e-7 / 6.2e-4 / 1.7e-4
  peaks 384/512/768: Omega 1.5395 / 1.5400 / 1.5404, |w| 14.43 / 15.08 / 15.36.
Time step, against C 0.5 (no -F; filter applied once per step):
  C 1    Omega 1e-9 / 9.7e-5 / 3.6e-5    wmax 2e-4 / 0.090 / 0.016    Gz 7e-8 / 5.5e-5 / 2.3e-5
  C 1.5  Omega 5e-9 / 1.5e-4 / 4.8e-5    wmax 3e-4 / 0.096 / 0.015
  C 2    Omega 2e-8 / 1.8e-4 / 6.1e-5    wmax 8e-4 / 0.095 / 0.016    Gz 1e-7 / 9.7e-5 / 4.1e-5
With -F 0.00125 (filter as a rate), against C 0.5 -F, to t* 2.1:
  C 1 -F  Omega 1e-9 / 1.2e-5 / 1.6e-6    wmax 2e-4 / 0.085 / 0.017    Gz 6e-9 / 6.0e-6 / 2e-7
  C 2 -F  Omega 2e-8 / 4.4e-5 / 5.5e-6    wmax 8e-4 / 0.085 / 0.014    Gz 4e-8 / 2.0e-5 / 5e-7
  C 0.5 without -F vs with: wmax 0.8%. So -F removes most of the dt dependence of Omega and Gz (1e-4 -> 1e-5), but
  not that of the grid max |w| in t* 1.5-2 (8.5% even between C 1 and C 0.5): the spikes shift in time.
Reading: Omega, E and Gz are converged in dt to 1e-4 for any C <= 2 and need ~512 at Re 5000 / ~768-1024 at Re 9000
for 0.1-0.5%. The grid max |w| during the reconnection spikes (t* 1.6-2.0) moves by 5-10% with C and with the grid:
it is a pointwise maximum of thin sheets that roll up, sensitive to any perturbation; its peak value is stable to
~2% (15.08 / 15.21 / 15.34 / 15.26 for C 1 / 1.5 / 2 / 0.5). Use C 1 with -F (dt-independent filter) for reference
runs; C 2 costs half the steps for the same integral quantities.
Cost of 384 -> 512 -> 768: the Yao grid ladder (n ~ Re) is consistent with Omega converged to 0.1-0.4%.

6. Spectral checkpoints and padded restart (k.STEP)
-------------------------------------------------
-k writes k.STEP (unnormalized double complex [3][nx/2+1][nz/2+1][ny/2+1], one pwrite per rank and component; with
k.STEP.t = t dt). -r k.STEP restarts exactly (final dump bitwise equal to the continuous run, 128^3); -r k.STEP -G
n0,n0,n0 zero-pads from a smaller grid (its Nyquist modes dropped) and sets dt anew. seg.sh does this per run dir.
Check (Re 5000, -S): 512 from t = 0 to t* 1.4, k. dump, padded to 768 and run to t* 2.6, against 768 from t = 0:
  Omega 1.8e-6, Gz 5.4e-7, E 3.3e-10 (t* 1.5-2); wmax 3% on the spikes. Dump/read: 1.6 / 0.7 s at 512 (2 GB file).
So the early, smooth phase can run on a coarse grid at no measurable cost in Omega or Gz.

7. Timings, one node (4 H200, seas_gpu holygpu8a16104: AMD 9374F, NVLink NV6, 4x NDR400), cube n^3 quarter box,
   -I yao, Re 40000 nu, fixed dt 1e-4, 20 steps (two diagnostic steps included), s/step
------------------------------------------------------------------------------------------------------------------
"old" = boxsymmp_b2 with BOXSYM_SYN6=1 NOFUSE=1 SELFCOPY=1 NL=0: the old code path (6 + 3 exchanges per stage, own
diagnostic synthesis, local copy in the exchange; arithmetic bitwise equal to the old binary, which cannot start from
-I yao). b2 = all changes of section 2; b3 = b2 + z row offsets once per tile row (extz 2.81 -> 2.18 ms, foldz 1.37
-> 1.05 ms per call at 1024^3 on one H200).
  n      old      new b2   new b3   -S b2    -S b3    compute only b2 / b3   exchanges only: old (double) / -S
  1024   0.2027   0.1810   0.1715   0.1672   0.1560   0.1602 / 0.1509        0.0664 / 0.0286
  1536   0.7809   0.7143   0.6841   0.6693   0.6325   0.6427 / 0.6067        0.2185 / 0.1057
  2048   1.6201   1.4441   1.3654   1.3303   1.2362   1.2798 / 1.2001        0.5134 / 0.2472
  old -> -S b3: 1.30x (1024), 1.23x (1536), 1.31x (2048). Rank memory (cudaMemGetInfo after setup): 10.7 / 34.2 /
  74.7 GB old, 11.2 / 36.0 / 79.0 GB new (+ b[6] for nl).
Inside a node the step is compute bound (compute alone = 97% of the -S step; the exchanges alone run at ~350-380 GB/s
per rank over NVLink and are hidden). nl_a/nl_b and the chunked boundaries (BOXSYM_HC=2, 4) change nothing here
(newSnl0 0.1654 / 0.6648 / 1.3181, hc2 = hc4 = newS within 0.1%): they are for the network-bound case.
Kernel shares at 1024^3 on one H200 (b3, nsys, 10 steps): z2z FFT 28%, ext/fold 47%, update 6%, C2R/R2C 8%,
cross0/cross12/spreadw/vortm 8%. The ext/fold kernels run at 60-75% of HBM bandwidth, the FFTs at ~85%.

8. Re 40000: resolution needed
------------------------------
Run r40k1024 (1024^3, -S -F 0.000625, C 1, 4 H200, t* 0 -> 3, 15581 steps, 0.154 s/step, 40 min): deliberately
under-resolved, for the velocity history and kmax eta. Omega/Omega0 peak 9.74 at t* 1.788 (paper 11.81 at 1.882),
|w| ratio peak 124 (paper 569); |w| already 7% below the paper at t* 1.0 and Omega 7% low at t* 1.4. So 1024 is not
adequate even before the collision.
Criterion from our own ladders, with keta = kc eta, kc = (n-1)/3 per direction (the code's column; box-mean eta at the
enstrophy peak; the Hou-Li filter keeps ~0.45 n, so this kc is conservative for our runs):
  Re 5000   384  keta_min 2.96: Omega -0.3%, |w| peak -6% (vs 768)
            512  3.96:            Omega 0.07%, |w| peak -1.8%
            768  5.94
  Re 9000   768  4.02:            Omega peak -0.15%, |w| peak -1.4% (vs 1024)
            1024 5.37
  Re 40000  1024 1.77 (with its own, too low, enstrophy; with the paper's peak Omega 11.8: 1.69)
keta ~ 4 gave Omega to ~0.1-0.2% and the |w| peak to ~2% at both Re 5000 and 9000. Scaling keta with n at fixed
dissipation (Omega/Omega0 = 11.8): keta = 1.69 n / 1024 -> keta 4 needs n = 2420, keta 5 needs 3030.
Chosen grid: 3072^3 (keta ~5.1 at the peak if the true peak Omega is the paper's; 4.8 if it is 15% higher), 2560^3
as the minimum (keta 4.2). Yao's 4096 x 2049^2 with 2/3 truncation has keta 6.8 in this measure. Not done (not needed
for a reference now that the production run is SEM/vm3): the 2048/2560/3072 ladder near the peak that would confirm
this from the run itself.
Time step: C 1 with -F (filter rate F ~ 0.00125 * 512 / n, about one full filter per step at the start). C 2 halves
the steps and moves Omega and Gz by ~4e-5 (Re 5000), but the reconnection-time grid max |w| by ~6-9% pointwise.
Velocity scale U = max(|u|+|v|+|w|) from r40k1024: 4.8 (t* 0-1), 5.6 (1.3), 7.8 (1.58), 12.5 (1.74), 9-11 (1.8-3.0).
dt = C / (kc U): at 3072, C 1: 2.0e-4 early, 7e-5 at the peak.

9. Cost projection, Re 40000 to t* 3 (steps from U(t) of r40k1024; scripts yao/proj.py)
---------------------------------------------------------------------------------------
Measured on one H200 node (section 7): 1536 0.633, 2048 1.236 s/step (-S, b3). Multi-node s/step below is a MODEL until
the queued slots report (see 10): compute = one-node compute x (n/2048)^3 (log factor) x 4/P; network = 40 float
exchanges of X1 Z1 Q x 8 B per step, the off-node share (P-4)/P at ~40 GB/s per GPU (4x NDR400 per node: 50 GB/s
nominal); s/step ~ the larger of the two plus ~10%.
  plan                                   grid, GPUs     s/step (model)   steps     wall h   GPU-h
  3072 from t = 0                        3072, 16       1.5              46700     19.5     311
  1536 to t* 1.0, padded to 3072         1536 4 / 3072 16  0.63 / 1.5    5360 + 35975   15.9   244
  same, 2560 after t* 1.0 (minimum)      2560, 16       0.95             30000     8.9      130
  same, 4096 after t* 1.0 (Yao's grid)   4096, 32       2.1              48000     28.9     900
  3072 with C 2                          3072, 16       1.5              2680 + 17987   8.0   122
Memory per rank (b3): 3072 on 8 H200 ~125 GB (fits, tight), on 16 ~63 GB; 4096 needs >= 24 H200 (16 would need ~150
GB). Checkpoints (k.): 174 GB at 3072, 413 GB at 4096 per dump (netscratch has ~1.2 TB free for the whole lab).
Caveats: U at 3072 near the peak is probably higher than at 1024 (sharper sheets), i.e. more steps; the multi-node
s/step is modelled. No spectral production run is planned (Re 40000 production: SEM or vm3); run/prod40k.sh +
seg.sh (phased, padded, resumable) stay as a reference tool.

10. Timings across nodes
------------------------
A100 (seas_gpu, 4 nodes x 4 A100-80GB, ONE 100 Gb/s HDR port per node = ~3 GB/s per GPU off node), 2048^3, 16 ranks,
s/step (b2):
  old 13.08   new 11.33   -S 5.79   -S without nl 5.90   -S HC 2 5.74   -S HC 4 5.72   -S + NCCL_NCHANNELS_PER_NET_PEER=4,
  P2P_NET_CHUNKSIZE=512k 6.43 (worse)   exchanges alone: old 12.43, -S 5.53 (3.9 GB/s per rank)   compute alone 0.73
  -> network bound: old -> -S 2.26x (45 double -> 40 float exchanges: 2.25x fewer bytes). The overlap is near
  complete: -S step = exchanges alone + 5% (compute alone is 13%). nl_a/nl_b gives 1.8%, the chunked boundaries
  another 1.0-1.3%. These nodes are unusable for this code at scale (5.8 s/step at 2048 vs 1.24 on one H200 node).
H200 2/4/8 nodes (4x NDR400 per node, 16x the A100 nodes' bandwidth per GPU): queued (yspec-b2 8 GPUs 2048 3072,
yspec-b4 16 GPUs 2048 3072 4096, yspec-b8 32 GPUs 3072 4096; estimated start 15:00 today to 06:30 tomorrow, all with
--exclude=holygpu8a12204). They run bslot.sh (bmat.sh with b2, bmat3.sh with b3), append to netscratch
yspec/bench.txt and cancel themselves. The first attempt on 2 nodes landed on holygpu8a12204 (cudaSetDevice busy).
Model used in section 9 for H200: step = max(compute, exchange) x 1.05-1.1, exchange = off-node bytes / 40 GB/s per
GPU, as calibrated on the A100 overlap above.
