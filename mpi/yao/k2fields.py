# k2fields.py k.STEP n m outprefix [q.STEP]: full 2 pi box fields from a boxsymmp k. dump (quarter box, grid n^3), spectrally
# truncated to m^3 (modes |p|, |q|, |r| < m/2 kept): writes outprefix_wabs.npy and outprefix_wz.npy, float16, shape (m, m, m),
# index order [i][j][k] = (x, y, z), x_i = -pi + i 2pi/m (same for y, z): x collision normal, y propagation, z tube axis.
# With q.STEP (and m = n) it also checks the velocity against the physical dump.
import sys, numpy as np
kf, n, m, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
X1 = Z1 = Q = n // 2 + 1; M = m // 2
K = np.memmap(kf, '<c16', 'r', shape=(3, X1, Z1, Q))
ox, oz = (1, 0, 0), (0, 0, 1)                  # u odd in x, w odd in z
s = (m / n) ** 3                                # unnormalized n-grid coefficients -> numpy irfftn on m^3
U = []
for c in range(3):
    a = np.array(K[c, :M, :M, :M]) * s          # [p][r][q], p, r = 0 .. M-1, q = 0 .. M-1
    A = np.zeros((m, m, M + 1), complex)        # [px][pz][qy], full x and z, half y
    sx = -1 if ox[c] else 1; sz = -1 if oz[c] else 1
    A[:M, :M, :M] = a
    A[m - M + 1:, :M, :M] = sx * a[1:][::-1]            # -p
    A[:M, m - M + 1:, :M] = sz * a[:, 1:][:, ::-1]      # -r
    A[m - M + 1:, m - M + 1:, :M] = sx * sz * a[1:, 1:][::-1, ::-1]
    U.append(A)
kx = np.fft.fftfreq(m, 1 / m)[:, None, None]; kz = np.fft.fftfreq(m, 1 / m)[None, :, None]; ky = np.arange(M + 1)[None, None, :]
inv = lambda A: np.roll(np.roll(np.fft.irfftn(A, s=(m, m, m), axes=(0, 1, 2)), M, 0), M, 1)   # [x][z][y], x, z from -pi
if len(sys.argv) > 5:
    q = np.memmap(sys.argv[5], '<f8', 'r', shape=(3, X1, Z1, n))   # [c][i][k][j], x = i h >= 0, z = k h, y = -pi + j h
    for c in range(3):
        f = inv(U[c])                                          # f[x][z][y]
        print("u%d max|recon - q| / max|q| on the quarter: %.2e" % (c, np.abs(f[M:, M:, :] - q[c][:M, :M, :]).max() / np.abs(q[c]).max()))
wx = 1j * (ky * U[2] - kz * U[1]); wy = 1j * (kz * U[0] - kx * U[2]); wz = 1j * (kx * U[1] - ky * U[0])
del U
fz = inv(wz); del wz
fa = fz ** 2
fa += inv(wx) ** 2; del wx
fa += inv(wy) ** 2; del wy
fa = np.sqrt(fa)
print("%s: max|w| %.6g, max w_z %.6g, min w_z %.6g, <|w|^2>/2 %.8g" % (kf, fa.max(), fz.max(), fz.min(), 0.5 * (fa ** 2).mean()))
np.save(out + "_wabs.npy", fa.transpose(0, 2, 1).astype(np.float16))   # [x][z][y] -> [x][y][z]
np.save(out + "_wz.npy", fz.transpose(0, 2, 1).astype(np.float16))
