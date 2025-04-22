import numpy as np
import math
import sys

a = np.memmap(sys.argv[1], dtype=float)
n = len(a) // 3
n = round(n**(1/3))
assert 3 * n * n * n == len(a)

U, V, W = np.reshape(a, (3, n, n, n))
print(np.min(U), np.max(U))
print(np.min(V), np.max(V))
print(np.min(W), np.max(W))
#plt.imshow(U[0])
#plt.show()
