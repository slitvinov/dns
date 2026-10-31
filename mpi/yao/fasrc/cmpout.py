# cmpout.py a.out b.out: max relative difference per column (common rows by step)
import sys
def rd(f):
    L = [l.split() for l in open(f) if l.strip()]
    h = L[0]; return h, {int(r[0]): [float(x) for x in r[1:]] for r in L[1:]}
ha, a = rd(sys.argv[1]); hb, b = rd(sys.argv[2])
ks = sorted(set(a) & set(b))
print("rows %d (steps %d..%d)" % (len(ks), ks[0], ks[-1]))
for i, n in enumerate(ha[1:]):
    if n not in hb: continue
    j = hb.index(n) - 1
    d = max(abs(a[k][i] - b[k][j]) / max(abs(a[k][i]), 1e-300) for k in ks)
    print("%-8s %.3e" % (n, d))
