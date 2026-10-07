import glob
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def spectrum(path, dk):
    t = float(open(path).readline().split()[-1])
    b, e = np.loadtxt(path).T
    k = np.arange(1, int(b[-1]) + 1)
    E = np.array([e[(b >= q - dk / 2) & (b < q + dk / 2)].sum() / dk for q in k])
    return t, k, E


def fit(k, E, lo, hi):
    s = (k >= lo) & (k <= hi) & (E > 1e-27)
    X = np.c_[np.ones(s.sum()), -np.log(k[s]), -2 * k[s]]
    a, n, d = np.linalg.lstsq(X, np.log(E[s]), rcond=None)[0]
    return d, n, np.exp(a)


table1 = {
    0.5: (1.107, 1.107, 4.31, 4.31),
    1.0: (0.453, 0.451, 5.02, 5.14),
    1.5: (0.193, 0.192, 4.79, 4.86),
    2.0: (0.080, 0.080, 4.48, 4.50),
    2.5: (0.020, 0.034, 4.71, 4.18),
    3.0: (-0.022, 0.005, 5.63, 4.59),
    3.5: (-0.007, -0.002, 4.51, 4.56),
}
dk = float(sys.argv[1]) if len(sys.argv) > 1 else 2
fits = {}
fig, ax = plt.subplots(1, 2, figsize=(10, 6))
for kmax, d, hi, a in (42, "inv042", 36, None), (84, "inv084", 75, ax):
    for path in sorted(glob.glob(d + "/e.*"))[1:]:
        t, k, E = spectrum(path, dk)
        lo, top = (4, 22) if abs(t - 0.5) < 1e-6 else (10, hi)
        fits[kmax, round(t, 2)] = fit(k, E, lo, top)
        if a is not None:
            for q, w in (0, 1), (1, dk):
                tt, kk, EE = spectrum(path, w)
                a[q].plot(kk[EE > 0], np.log10(EE[EE > 0]), ".", ms=3)
for q, w in (0, 1), (1, dk):
    ax[q].set_xlim(0, 85)
    ax[q].set_ylim(-30, 0)
    ax[q].set_xlabel("wavenumber k")
    ax[q].set_ylabel("log E(k)")
    ax[q].set_title("inviscid, kmax = 84, dk = %g" % w)
fig.savefig("img/spectrum.svg")
print("Table 1: delta(t), n(t) fitted, dk = %g: ours (paper)" % dk)
print("%4s %16s %16s %14s %14s" % ("t", "delta 42", "delta 84", "n 42", "n 84"))
for t, (d42, d84, n42, n84) in table1.items():
    a42, b84 = fits.get((42, t)), fits.get((84, t))
    if a42 and b84:
        print("%4.1f %7.3f (%6.3f) %7.3f (%6.3f) %5.2f (%5.2f) %5.2f (%5.2f)"
              % (t, a42[0], d42, b84[0], d84, a42[1], n42, b84[1], n84))
plt.figure(figsize=(6, 6))
for kmax, m in (84, "o"), (42, "+"):
    T = sorted(t for q, t in fits if q == kmax and fits[q, t][0] > 0)
    plt.plot(T, [np.log(fits[kmax, t][0]) for t in T], m, mfc="none",
             label="kmax = %d" % kmax)
t = np.linspace(0, 4)
plt.plot(t, np.log(2.6 * np.exp(-t / 0.57)), "k-", label="(3.6): 2.6 exp(-t/0.57)")
plt.axhline(np.log(np.pi / 84), color="gray")
plt.xlabel("t")
plt.ylabel("ln delta(t)")
plt.legend()
plt.savefig("img/delta.svg")
fig, ax = plt.subplots(1, 2, figsize=(11, 5.5))
for r in "0200", "0400", "0800", "1600", "inf":
    d = np.loadtxt("inv084/out" if r == "inf" else "tg0256/" + r)
    nu = 0 if r == "inf" else 1 / int(r)
    t, O, P = d[:, 1], d[:, 3], d[:, 12]
    S = np.sqrt(135 / 98) * (np.gradient(O, t) + 2 * nu * P) / O**1.5
    for a, y in (ax[0], d[:, 4]), (ax[1], S):
        a.plot(t, y, "k" if r == "inf" else "-", label="R = " + r.lstrip("0"))
for a, title in (ax[0], "(5.6): from dvx/dx"), (ax[1], "(5.8): isotropic, from spectra"):
    a.set_xlim(0, 10)
    a.set_ylim(0, 1.2)
    a.set_xlabel("t")
    a.set_ylabel("S3(0)")
    a.set_title(title)
    a.legend()
fig.savefig("img/skewness.svg")
table8 = {"0200": (7, 0.45, 6.8, 9.9, 86, 1.4e3, 18, 773),
          "0400": (7, 0.61, 6.7, 13.9, 99, 2.4e3, 17.26, 823),
          "0800": (7, 0.47, 8.8, 17.6, 250, 1.4e4, None, None),
          "1600": (9, 0.65, 10.0, 23.1, 273, 1.9e4, 15.6, 660)}
print("Table 8: ours (paper); columns S3 S4 S5 S6 S8 Sb4 Sb6")
for r, (t, *p) in table8.items():
    d = np.loadtxt("tg0256/" + r)
    i = np.argmin(abs(d[:, 1] - t))
    ours = d[i, [4, 5, 6, 7, 9, 10, 11]]
    print("R = %4d t = %d " % (int(r), t) + " ".join(
        "%.3g (%s)" % (o, "-" if q is None else "%g" % q) for o, q in zip(ours, p)))


def side(name, page, draw):
    fig = plt.figure(figsize=(13, 7))
    a = fig.add_subplot(1, 2, 1)
    a.imshow(plt.imread("refs/pages/fig_%s.png" % page))
    a.axis("off")
    a.set_title("Brachet et al. (1983)")
    draw(fig)
    fig.savefig("img/%s.png" % name, dpi=110)


marks = "x1+s^odv"


def fig3(fig):
    a = fig.add_subplot(1, 2, 2)
    for i, path in enumerate(sorted(glob.glob("inv084/e.*"))[1:8]):
        t, k, E = spectrum(path, 1)
        s = E > 1e-30
        a.plot(k[s], np.log10(E[s]), marks[i], mfc="none", ms=4, color="k",
               label="t = %.1f" % t)
    a.set_xlim(0, 80)
    a.set_ylim(-30, 0)
    a.set_xlabel("wavenumber k")
    a.set_ylabel("log E(k)")
    a.set_title("tg -M 128, inviscid, dk = 1")
    a.legend(fontsize=8)


def fig4(fig):
    for q in 0, 1:
        a = fig.add_subplot(2, 2, 2 + 2 * q)
        for i, path in enumerate(sorted(glob.glob("inv084/e.*"))[1:8]):
            t, k, E = spectrum(path, 2)
            s = (E > 1e-30) & (k <= 80)
            x = k if q == 0 else np.log10(k)
            a.plot(x[s], np.log10(E[s]), marks[i], mfc="none", ms=3, color="k")
            lo, hi = (4, 22) if abs(t - 0.5) < 1e-6 else (10, 75)
            d, n, A = fit(k, E, lo, hi)
            f = (k >= lo) & (k <= hi)
            a.plot(x[f], np.log10(A * k[f] ** -n * np.exp(-2 * d * k[f])), "r-",
                   lw=1)
        a.set_xlim((0, 80) if q == 0 else (0, 2))
        a.set_ylim(-30, 0)
        a.set_xlabel("k" if q == 0 else "log k")
        a.set_ylabel("log E(k)")
    fig.axes[1].set_title("tg -M 128, inviscid, dk = 2, fits (3.5) in red")


def fig10(fig):
    for i, (r, t) in enumerate((("1600", 9), ("3000", 9), ("3000", 5))):
        a = fig.add_subplot(3, 2, 2 + 2 * i)
        tt, k, E = spectrum("tgspec/%s/e.%08d" % (r, round(t / 0.0025)), 1)
        s = E > 0
        a.semilogy(k[s], E[s], "k.", ms=3, label="linear-log")
        a.semilogy(50 * np.log10(k[s]), E[s], "b.", ms=3, label="log-log")
        d, n, A = fit(k, E, 13, 83)
        f = (k > 13) & (k < 83)
        g = A * k[f] ** -n * np.exp(-2 * d * k[f])
        a.semilogy(k[f], g, "r-", lw=1)
        a.semilogy(50 * np.log10(k[f]), g, "r-", lw=1)
        a.set_xlim(0, 100)
        a.set_ylim(1e-10, 1e-1)
        a.set_title("R = %s, t = %d: n = %.2f, delta = %.3f" % (r, t, n, d),
                    fontsize=9)
        if i == 0:
            a.legend(fontsize=7)
    a.set_xlabel("k (black) or 50 log k (blue)")
    fig.subplots_adjust(hspace=0.45)


side("fig3", "09", fig3)
side("fig4", "10", fig4)
side("fig10", "24", fig10)
