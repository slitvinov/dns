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
