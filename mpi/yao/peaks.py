# peaks.py out [out ...]: Yao case summary of boxsymmp out files: Omega/Omega0 and |w|max ratio peaks (parabola through
# the three samples around the largest), Gamma_z/Gamma0 crossing times, dE/dt vs -2 nu Omega; t* = 0.465772 t.
# Paper (fig 23 vector paths, ~/dns/flsr2/yao) peaks for Re 5000/9000/20000/40000 printed with -p.
import sys, numpy as np
TS = 7.680385244 / (2 * np.pi * 1.62 ** 2)
G0 = 7.680385244
def peak(t, y):
    i = int(np.argmax(y))
    if 0 < i < len(y) - 1:
        a, b, c = y[i - 1], y[i], y[i + 1]; d = a - 2 * b + c
        x = 0.5 * (a - c) / d if d != 0 else 0
        h = t[i + 1] - t[i]
        return b - 0.25 * (a - c) * x, t[i] + x * h
    return y[i], t[i]
def cross(t, y, v):
    j = np.where((y[:-1] > v) & (y[1:] <= v))[0]
    if len(j) == 0: return np.nan
    j = j[0]; return t[j] + (v - y[j]) / (y[j + 1] - y[j]) * (t[j + 1] - t[j])
args = sys.argv[1:]
if args and args[0] == '-p':
    args = args[1:]
    for re in (5000, 9000, 20000, 40000):
        a = np.loadtxt('/Users/lisergey/dns/flsr2/yao/yao23a_Re%d.txt' % re); c = np.loadtxt('/Users/lisergey/dns/flsr2/yao/yao23c_Re%d.txt' % re)
        ia = np.argmax(a[:, 1]); ic = np.argmax(c[:, 1])
        print("paper Re %5d: Omega peak %.4f at t* %.3f, |w| ratio peak %.2f at t* %.3f" % (re, a[ia, 1], a[ia, 0], c[ic, 1], c[ic, 0]))
for f in args:
    L = [l.split() for l in open(f) if l.strip()]
    h = L[0]; d = np.array([[float(x) for x in r] for r in L[1:]])
    col = {n: d[:, i] for i, n in enumerate(h)}
    t, E, Om, w = col['t'], col['E'], col['Omega'], col['wmax']
    ts = TS * t
    po, to = peak(ts, Om / Om[0]); pw, tw = peak(ts, w / w[0])
    s = "%s: t* to %.3f, Omega peak %.4f at t* %.4f, |w| ratio peak %.3f at t* %.4f" % (f, ts[-1], po, to, pw, tw)
    if 'Gz' in col:
        g = col['Gz'] / G0
        s += ", Gz/G0 0.5 at t* %.4f, 0.2 at %.4f, min %.3f" % (cross(ts, g, 0.5), cross(ts, g, 0.2), g.min())
    print(s)
