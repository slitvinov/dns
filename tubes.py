import sys
import math
import numpy as np

try:
    import scipy.fft as fft

    kw = {"workers": -1}
except ImportError:
    fft = np.fft
    kw = {}

if len(sys.argv) != 5:
    sys.stderr.write("usage: tubes.py nx ny nz u.raw\n")
    sys.exit(1)
nx, ny, nz = map(int, sys.argv[1:4])
out = sys.argv[4]
Lx, Ly, Lz = 6 * math.pi, 4 * math.pi, 2 * math.pi
A, rc, w0 = 0.2, 0.666, 26.093
K = 0.5 * math.exp(2) * math.log(2)
x = -Lx / 2 + Lx / nx * np.arange(nx)
y = -Ly / 2 + Ly / ny * np.arange(ny)
z = Lz / nz * np.arange(nz)


def profile(r):
    s = r / rc
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        f = w0 * (1 - np.exp(-K / s * np.exp(1 / (s - 1))))
    return np.where(s < 1, np.where(s > 0, f, w0), 0)


r = np.linspace(0, rc, 200001)
G = 2 * math.pi * np.trapezoid(profile(r) * r, r)
sys.stderr.write("tubes.py: circulation %.6f, Re %.1f (nu = 0.001)\n" % (G, G / 0.001))

wh = []
for d in range(3):
    o = np.zeros((nx, ny, nz))
    for xc, al, sg in (-0.866, math.pi / 3, 1), (0.866, 2 * math.pi / 3, -1):
        cx = xc + A * math.cos(al) * (1 + np.cos(z))
        cy = A * math.sin(al) * (1 + np.cos(z))
        f = sg * profile(np.sqrt((x[:, None, None] - cx) ** 2 + (y[None, :, None] - cy) ** 2))
        t = (-A * math.cos(al) * np.sin(z), -A * math.sin(al) * np.sin(z), np.ones(nz))[d]
        o += f * t
    wh.append(fft.rfftn(o, **kw))
    del o
kx = 2 * math.pi / Lx * np.fft.fftfreq(nx, 1 / nx)[:, None, None]
ky = 2 * math.pi / Ly * np.fft.fftfreq(ny, 1 / ny)[None, :, None]
kz = 2 * math.pi / Lz * np.arange(nz // 2 + 1)[None, None, :]
kk = kx ** 2 + ky ** 2 + kz ** 2
kk[0, 0, 0] = 1
with open(out, "wb") as f:
    for d in range(3):
        a, b = ((1, 2), (2, 0), (0, 1))[d]
        k = (kx, ky, kz)
        uh = 1j * (k[a] * wh[b] - k[b] * wh[a]) / kk
        uh[0, 0, 0] = 0
        u = fft.irfftn(uh, s=(nx, ny, nz), **kw)
        sys.stderr.write("tubes.py: u%d max %.4f\n" % (d, np.abs(u).max()))
        u.tofile(f)
