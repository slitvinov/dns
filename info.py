import numpy as np
import math
import sys

nv = 4
a = np.memmap(sys.argv[1], dtype=float)
n = len(a) // nv
n = round(n**(1 / 3))
assert nv * n * n * n == len(a)

U, V, W, P = np.reshape(a, (nv, n, n, n))
print(np.min(U), np.max(U))
print(np.min(V), np.max(V))
print(np.min(W), np.max(W))
print(np.min(P), np.max(P))
