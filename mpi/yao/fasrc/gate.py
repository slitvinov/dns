# gate.py: the gates of the Yao Re 40000 fine phase; writes ../gate.ok or ../gate.fail (with the numbers) in yspec/
#  R1 Re 9000 coarse-to-fine (r9kpad out.0 + out.1 vs r9k1024): Omega peak < 0.1%, |w| peak < 2%
#  R2 Re 9000 -S (r9k768 vs r9k768noS): Omega peak < 0.2%, |w| peak < 2%
#  G1 Re 20000 (r20k1536 vs r20k2048p = 1536 to t* 1.3 padded to 2048, to t* 2.0): Omega peak < 0.3%, |w| peak < 3%, and
#     the 2048 Omega peak within 0.5% of the paper (5.5146, fig 23a exact path)
#  A file yspec/hold keeps the gate closed: the numbers go to gate.held instead (hold removed by hand after review)
import sys, os, numpy as np
Y = "/n/netscratch/koumoutsakos_lab/Everyone/slitvinov/yspec"; TS = 0.4657717102700246; O0, W0 = 2.6026367, 20.0
def rd(*fs):
    rows = {}
    for f in fs:
        L = [l.split() for l in open(f) if l.strip()]
        for r in L[1:]:
            rows[int(r[0])] = [float(x) for x in r]
    d = np.array([rows[k] for k in sorted(rows)]); return d[:, 1] * TS, d[:, 3] / O0, d[:, 5] / W0
def peak(t, y):
    i = int(np.argmax(y))
    if 0 < i < len(y) - 1:
        a, b, c = y[i - 1], y[i], y[i + 1]; dd = a - 2 * b + c; x = 0.5 * (a - c) / dd if dd else 0
        return b - 0.25 * (a - c) * x
    return y[i]
def cmp(name, A, B, to, tw, tmax=None):
    ta, oa, wa = A; tb, ob, wb = B
    if tmax is None: tmax = min(ta[-1], tb[-1])
    ma, mb = ta <= tmax, tb <= tmax
    do = abs(peak(ta[ma], oa[ma]) / peak(tb[mb], ob[mb]) - 1); dw = abs(peak(ta[ma], wa[ma]) / peak(tb[mb], wb[mb]) - 1)
    ok = do < to and dw < tw
    return ok, "%s: Omega peak %.4g vs %.4g (%.3f%%, limit %.1f%%), |w| peak %.4g vs %.4g (%.2f%%, limit %.0f%%): %s" % (
        name, peak(ta[ma], oa[ma]), peak(tb[mb], ob[mb]), 100 * do, 100 * to, peak(ta[ma], wa[ma]), peak(tb[mb], wb[mb]), 100 * dw, 100 * tw, "pass" if ok else "FAIL")
R = Y + "/run/"
res = [cmp("R1 Re 9000 512->1024 padded at t* 1.3 vs 1024", rd(R + "r9kpad/out.0", R + "r9kpad/out.1"), rd(R + "r9k1024/out"), 0.001, 0.02),
       cmp("R2 Re 9000 768 -S vs double", rd(R + "r9k768/out"), rd(R + "r9k768noS/out"), 0.002, 0.02),
       cmp("G1 Re 20000 1536 vs 2048", rd(R + "r20k1536/out"), rd(R + "r20k2048p/out.0", R + "r20k2048p/out.1"), 0.003, 0.03)]
t2, o2, w2 = rd(R + "r20k2048p/out.0", R + "r20k2048p/out.1"); t1, o1, w1 = rd(R + "r20k1536/out")
m1 = t1 <= t2[-1]; po2, po1 = peak(t2, o2), peak(t1[m1], o1[m1])
okp = abs(po2 / 5.5146 - 1) < 0.005
res.append((okp, "G1 vs paper: Omega/Omega0 peak 1536 %.4f -> 2048 %.4f, paper 5.5146 (2048 %+.2f%%, limit 0.5%%); |w| peak 1536 %.1f -> 2048 %.1f, paper 202: %s" % (
    po1, po2, 100 * (po2 / 5.5146 - 1), peak(t1[m1], w1[m1]), peak(t2, w2), "pass" if okp else "FAIL")))
ok = all(r[0] for r in res)
held = os.path.exists(Y + "/hold")
open(Y + ("/gate.held" if held else "/gate.ok" if ok else "/gate.fail"), "w").write("\n".join(r[1] for r in res) + ("\nall pass\n" if ok else "\nNOT all pass\n"))
print("\n".join(r[1] for r in res)); print("gate:", "held" if held else "ok" if ok else "FAIL")
