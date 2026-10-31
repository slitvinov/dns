# davg.py d.STEP_X.npy out.npy: z block mean (factor 4, weights 1/2 at the ends, mirrored at z = 0 and pi) of a boxsymmp
# -V output (shape (Z1, Xc, Yc), already block-averaged in x, y): result (Zc, Xc, Yc) float16, order [K][I][J],
# z = K 4 dz, x = I 4 dx, y = -pi + J 4 dy: the whole symmetry quarter at 2 pi / (n / 4)
import sys, numpy as np
a = np.load(sys.argv[1], mmap_mode='r'); Z1 = a.shape[0]; Zc = (Z1 - 1) // 4 + 1; n2 = 2 * (Z1 - 1)
w = {-2: .5, -1: 1, 0: 1, 1: 1, 2: .5}
o = np.zeros((Zc,) + a.shape[1:], np.float32)
for K in range(Zc):
    for b, wb in w.items():
        k = K * 4 + b; k = -k if k < 0 else (n2 - k if k > Z1 - 1 else k)
        o[K] += wb * a[k].astype(np.float32)
np.save(sys.argv[2], (o / 4).astype(np.float16))
