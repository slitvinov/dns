# ladder.py out...: peaks of Omega/Omega0 and |w|max/|w|0 with the t = 0 values of the IC (Omega0 2.6026367, |w|0 20.0),
# for runs that start late (padded rungs); pairwise peak differences against the last file
import sys, numpy as np
TS = 0.4657717102700246; O0, W0 = 2.6026367, 20.0
def rd(f):
    L = [l.split() for l in open(f) if l.strip()]; h = L[0]; d = np.array([[float(x) for x in r] for r in L[1:]])
    return {n: d[:, i] for i, n in enumerate(h)}
def peak(t, y):
    i = int(np.argmax(y))
    if 0 < i < len(y) - 1:
        a, b, c = y[i - 1], y[i], y[i + 1]; dd = a - 2 * b + c; x = 0.5 * (a - c) / dd if dd else 0
        return b - 0.25 * (a - c) * x, t[i] + x * (t[i + 1] - t[i])
    return y[i], t[i]
R = []
for f in sys.argv[1:]:
    c = rd(f); ts = c['t'] * TS
    po, to = peak(ts, c['Omega'] / O0); pw, tw = peak(ts, c['wmax'] / W0)
    R.append((f, po, to, pw, tw, ts[-1], c['keta'].min()))
    print("%-28s t* to %.3f  Omega/Omega0 peak %.4f at %.4f   |w| ratio peak %.1f at %.4f   keta_min %.2f" % (f, ts[-1], po, to, pw, tw, c['keta'].min()))
ref = R[-1]
for r in R[:-1]:
    print("%s vs %s: Omega peak %+.2f%%, |w| peak %+.1f%%" % (r[0], ref[0], 100 * (r[1] / ref[1] - 1), 100 * (r[3] / ref[3] - 1)))
# extrapolation of the peaks when exactly three rungs on an equally spaced n ladder are given (e.g. 2048, 2560, 3072):
# Aitken delta^2 (exact for geometric = spectral convergence in n) and the algebraic fit phi = phi_inf + A n^-p;
# error bar = the larger of |phi_inf(Aitken) - phi_inf(algebraic)| and |phi(finest) - phi_inf(Aitken)| / 2
if len(R) == 3:
    import re
    ns = [int(re.findall(r'(\d{3,4})', r[0])[-1]) for r in R]
    for name, j in (('Omega/Omega0 peak', 1), ('|w| ratio peak', 3)):
        f1, f2, f3 = (r[j] for r in R)
        d1, d2 = f2 - f1, f3 - f2
        ait = f3 - d2 * d2 / (d2 - d1) if d2 != d1 else f3
        s = "%s: %s -> " % (name, " / ".join("%.4g" % f for f in (f1, f2, f3)))
        alg = None
        if d1 * d2 > 0 and abs(d2) < abs(d1):
            from scipy.optimize import brentq
            g = lambda p: (ns[1] ** -p - ns[2] ** -p) / (ns[0] ** -p - ns[1] ** -p) - d2 / d1
            try:
                p = brentq(g, 0.05, 40); A = d1 / (ns[1] ** -p - ns[0] ** -p); alg = f3 - A * ns[2] ** -p
                s += "algebraic p %.2f: %.4g; " % (p, alg)
            except ValueError:
                pass
        err = max(abs(f3 - ait) / 2, abs(ait - alg) if alg is not None else abs(f3 - ait))
        print(s + "Aitken %.4g +- %.2g (finest grid %.2f%% from it)" % (ait, err, 100 * (f3 / ait - 1)))
