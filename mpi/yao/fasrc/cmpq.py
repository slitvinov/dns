# cmpq.py a.q b.q nx ny nz: max |a - b| / max |a| per component and in all of a q. file
import sys, numpy as np
nx, ny, nz = map(int, sys.argv[3:6])
sh = (3, nx // 2 + 1, nz // 2 + 1, ny)
a = np.memmap(sys.argv[1], "<f8", "r", shape=sh); b = np.memmap(sys.argv[2], "<f8", "r", shape=sh)
m = np.abs(a).max()
print("q max|a-b|/max|a|: " + " ".join("%.3e" % (np.abs(a[c] - b[c]).max() / m) for c in range(3)))
