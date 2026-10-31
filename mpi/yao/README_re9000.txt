Yao & Hussain 2020a anti-parallel tubes, Re = Gamma0/nu = 9000 (nu 8.533761e-4), boxsymmp spectral run (quarter box by
the two mirror symmetries, 1024^3, Hou-Li filter, adaptive dt, Courant 1, -S), validated against the paper (Omega peak
2.271 at t* 1.835 vs 2.288 at 1.829). Spectral agent, 2026-10-31.
Files re9000_tstar<T>_{wabs,wz,lam2}.npy, T = t* = t Gamma0 / (2 pi b^2) = 0.4657717 t, at t* 0, 0.25, ..., 2.5.
Each: numpy float16, shape (512, 512, 512), C order [i][j][k] = (x, y, z): the full 2 pi periodic box, unfolded from the
quarter by both mirror symmetries and spectrally truncated to 512^3 (modes |k| < 256 per direction kept).
  x_i = -pi + i 2 pi/512, y_j = -pi + j 2 pi/512, z_k = -pi + k 2 pi/512 (index 256 = coordinate 0).
  x: collision normal (tube 1 at x = -0.81 with omega_z > 0, tube 2 at x = +0.81 with omega_z < 0; x = 0 collision plane)
  y: propagation direction of the pair;  z: tube axis (z = 0 symmetry plane).
wabs = |omega|; wz = omega_z (signed); lam2 = lambda2 / 1000 (lambda2 = middle eigenvalue of S^2 + W^2, negative in cores).
|omega| at t = 0 is 20 at the tube centres. Origin of the data: FASRC netscratch slitvinov/yspec/run/r9kvid.
