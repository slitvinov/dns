# plot.py out.png RE label=file [label=file ...]: Omega/Omega0, |w|max ratio, Gamma_z/Gamma0 vs t* for boxsymmp out files
# (several out.STEP of one run may be given as label=f1,f2,...), with the paper curves of fig 23 (Re RE) and fig 22a.
import sys, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
TS = 7.680385244 / (2 * np.pi * 1.62 ** 2); G0 = 7.680385244
png, re = sys.argv[1], int(sys.argv[2])
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
P = '/Users/lisergey/dns/flsr2/yao/yao23%s_Re%d.txt'
for i, k in enumerate('ac'):
    a = np.loadtxt(P % (k, re)); a = a[np.argsort(a[:, 0])]
    ax[i].plot(a[:, 0], a[:, 1], 'k.', ms=2, label='Yao 2020 Re %d' % re)
try:
    a = np.loadtxt('/Users/lisergey/dns/bench/ref/yao/dig/yao_fig22a_Re%d.txt' % re, comments='#')
    ax[2].plot(a[:, 0], a[:, 1], 'k.', ms=3, label='Yao fig 22a (colour trace)')
except Exception:
    pass
for arg in sys.argv[3:]:
    lab, fs = arg.split('=')
    rows = {}
    for f in fs.split(','):
        L = [l.split() for l in open(f) if l.strip()]
        h = L[0]
        for r in L[1:]:
            rows[int(r[0])] = [float(x) for x in r[1:]]
    st = sorted(rows); d = np.array([rows[s] for s in st]); col = {n: d[:, i - 1] for i, n in enumerate(h) if i > 0}
    ts = TS * col['t']
    ax[0].plot(ts, col['Omega'] / col['Omega'][0], lw=1, label=lab)
    ax[1].semilogy(ts, col['wmax'] / col['wmax'][0], lw=1, label=lab)
    if 'Gz' in col:
        ax[2].plot(ts, col['Gz'] / G0, lw=1, label=lab)
ax[1].set_yscale('log')
for a, y in zip(ax, ['Omega/Omega0', '|w|max/|w|max,0', 'Gamma_z/Gamma0']):
    a.set_xlabel('t*'); a.set_ylabel(y); a.set_xlim(0, 3.7); a.legend(fontsize=7)
fig.tight_layout(); fig.savefig(png, dpi=110)
