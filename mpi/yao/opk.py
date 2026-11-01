# opk.py out [out ...]: Omega/Omega0 and |w|max/20 peaks of boxsymmp out files, Omega0 = 2.6026367 fixed (the IC value; right
# also for runs restarted mid-way, unlike peaks.py which divides by the first row); parabola through the three top samples
import sys, numpy as np
TS = 0.4657717102700246
def pk(t, y):
    i = int(np.argmax(y)); a, b, c = y[i - 1:i + 2]; d = a - 2 * b + c; x = .5 * (a - c) / d
    return b - .25 * (a - c) * x, t[i] + x * (t[i + 1] - t[i])
for f in sys.argv[1:]:
    a = np.loadtxt(f, skiprows=1); t = a[:, 1] * TS
    o, to = pk(t, a[:, 3] / 2.6026367); w, tw = pk(t, a[:, 5] / 20)
    print("%s: Omega/Omega0 peak %.5f at t* %.4f, |w|max ratio peak %.3f at t* %.4f, last t* %.3f" % (f, o, to, w, tw, t[-1]))
