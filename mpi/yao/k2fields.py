# k2fields.py k.STEP n m outprefix [q.STEP]: full 2 pi box fields from a boxsymmp k. dump (quarter box, grid n^3), spectrally
# truncated to m^3 (modes |p|, |q|, |r| < m/2 kept): writes outprefix_wabs.npy, outprefix_wz.npy and outprefix_lam2.npy
# (lambda2 of S^2 + W^2 divided by 1000, negative in vortex cores), float16, shape (m, m, m),
# index order [i][j][k] = (x, y, z), x_i = -pi + i 2pi/m (same for y, z): x collision normal, y propagation, z tube axis.
# k.STEP may also be a .slab file: K[:, :m/2, :m/2, :m/2] of the dump only, same n.
# With q.STEP (and m = n) it also checks the velocity against the physical dump.
import sys, numpy as np
kf, n, m, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
X1 = Z1 = Q = n // 2 + 1; M = m // 2
# a .slab file holds only the needed low modes [3][M][M][M] (K[:, :M, :M, :M] of the dump, cut where numpy is missing)
K = np.memmap(kf, '<c16', 'r', shape=(3, M, M, M) if kf.endswith('.slab') else (3, X1, Z1, Q))
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
kk = (kx, ky, kz)
G = [[None] * 3 for _ in range(3)]                       # G[a][b] = d_b u_a, physical [x][z][y]
for a in range(3):
    for b in range(3):
        G[a][b] = inv(1j * kk[b] * U[a])
del U
ox = G[2][1] - G[1][2]; oy = G[0][2] - G[2][0]; fz = G[1][0] - G[0][1]
fa = np.sqrt(ox ** 2 + oy ** 2 + fz ** 2); del ox, oy
l2 = np.empty(fa.shape)
for i in range(m):                                       # lambda2 plane by plane
    g_ = np.stack([np.stack([G[a][b][i] for b in range(3)], -1) for a in range(3)], -2)   # [z][y][a][b]
    S = 0.5 * (g_ + np.swapaxes(g_, -1, -2)); W = 0.5 * (g_ - np.swapaxes(g_, -1, -2))
    l2[i] = np.linalg.eigvalsh(S @ S + W @ W)[..., 1]
del G
print("%s: max|w| %.6g, max w_z %.6g, min w_z %.6g, <|w|^2>/2 %.8g, min lambda2 %.6g" % (kf, fa.max(), fz.max(), fz.min(), 0.5 * (fa ** 2).mean(), l2.min()))
np.save(out + "_wabs.npy", fa.transpose(0, 2, 1).astype(np.float16))   # [x][z][y] -> [x][y][z]
np.save(out + "_wz.npy", fz.transpose(0, 2, 1).astype(np.float16))
np.save(out + "_lam2.npy", (l2 / 1000).transpose(0, 2, 1).astype(np.float16))
