import numpy as np
import math

L = 2 * math.pi
nv = 4
n = 32
dx = L / n
a = np.memmap("tgv.raw", dtype=float, mode="w+", shape=(nv, n, n, n))
U, V, W, P = a
i = np.arange(n)
s = np.sin(dx * i)
c = np.cos(dx * i)
np.einsum('i,j,k', s, c, c, out=U)
np.einsum('i,j,k', -c, s, c, out=V)
P.fill(0)
W.fill(0)
