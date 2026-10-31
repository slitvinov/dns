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
  synth5       RK stages synthesize u from 3 exchanges plus d_z u, d_z v (2), and form omega in the mixed layout
               after the exchange (x, y spectral there): 5 instead of 6 exchanges per synthesis, 8 instead of 9 per
               stage. BOXSYM_SYN6=1 restores the old path.
  -k, -r k.STEP [-G nx,ny,nz]   spectral checkpoints (exact restart; zero-padded restart on a larger grid).
  BOXSYM_A2AONLY=K, BOXSYM_NOCOMM=1   timing only: exchanges alone / compute alone.
Checks:
  old binary vs new with BOXSYM_SYN6=1: q dump bitwise equal.
  synth5 vs old: q max|diff|/max 3.7e-16, out columns <= 2e-15 (tail columns are rounding-level quantities).
  np 2 vs np 1: q bitwise equal.  -S np 2: repeatable bitwise; NCCL and MPI paths bitwise equal.
  -S vs double after 100 steps: E 1.8e-11, Omega 4.4e-11, wmax 1.9e-9, Gz 8.6e-11 relative. The z tail
  fraction (energy beyond 0.9 kmax) rises from 7e-22 to 1.8e-19: the float rounding floor of the spectrum.
  adaptive restart (q.): same dt sequence and times; fields equal to rounding (the q. read re-analyses).
