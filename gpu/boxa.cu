#include <cufft.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { nthread = 256, nblock = 1024 };
static const double pi = 3.141592653589793238;
static const double A[] = {0, -567301805773.0 / 1357537059087.0,
                           -2404267990393.0 / 2016746695238.0,
                           -3550918686646.0 / 2091501179385.0,
                           -1275806237668.0 / 842570457699.0};
static const double B[] = {1432997174477.0 / 9575080441755.0,
                           5161836677717.0 / 13612068292357.0,
                           1720146321549.0 / 2090206949498.0,
                           3134564353537.0 / 4481467310338.0,
                           2277821191437.0 / 14882151754819.0};
/* boxa: box.cu + the Kelvin-filtered (alpha) model, math/Q23_alpha_exact.md.
   v_t = P(u x curl v) + nu Lap v + F,  u = G v (one Fourier multiplier). Its curl is q_t + u.grad q = q.grad u + nu Lap q
   with q = curl v: vm3's raw-w model (Q17). -G 0 Navier-Stokes, 1 Helmholtz (LANS-alpha, -D alpha), 2 Gaussian-Laguerre
   order 2p = 20 (vm3 core, -D delta), 3 radial table (-T file: k delta, G; -D delta; e.g. the compact core).
   -K k0 -e E0: Taylor-Green forcing at k0 with energy control (amplitude keeps E_alpha at E0, relaxation time -c tau).
   -i zero: start from rest plus a small random low-k field (-R amplitude). -S dt_s: shell spectra lines "#s".
   Subfilter closures (math/Q26_closure.md), combinable:
   -L 1      Leray-alpha: v_t + u.grad v = -grad p + ..., N = -P div(u v) (the (grad u)^T v term dropped; invariant 1/2 |v|^2)
   -A n      ADM / van Cittert on the advecting velocity: u = G_n v, G_n = 1 - (1 - G)^(n+1) (default 0: u = G v)
   -B chi    relaxation (plain sink) chi (1 - G)^(n+1) v (Stolz-Adams; n from -A)
   -W C -P p strain-scaled band sink: v_t -= C H^(1/2) [s(x) H^(1/2) v], s = sqrt(2/3 S:S) >= max eig S (traceless S of u),
             H = (1 - G)^p (default p 0.5) */
struct G {
  long nx, ny, nz, nzf, nf, nr;
  double ax, ay, az;
  int f;
  double fp;
  int gk;        /* filter kind */
  double gd;     /* alpha or delta */
  const double *gt; /* table of G on a uniform grid in k delta, step gdk, gn entries (kind 3) */
  double gdk;
  int gn;
  double fa;     /* forcing amplitude of this step */
  int an;        /* -A: ADM order */
  double chi;    /* -B: relaxation rate */
  double wc, wp; /* -W, -P: strain-scaled band sink */
};
__device__ static double filt(G g, double kk);
/* the advecting velocity's filter: G_n = 1 - (1 - G)^(n+1) */
__device__ static double gu(G g, double kk) {
  double G0 = filt(g, kk);
  return g.an > 0 ? 1 - pow(1 - G0, g.an + 1) : G0;
}
__device__ static double hband(G g, double kk) { return pow(fmax(0.0, 1 - filt(g, kk)), 0.5 * g.wp); } /* H^(1/2) */
__device__ static double filt(G g, double kk) {
  if (g.gk == 1)
    return 1 / (1 + g.gd * g.gd * kk);
  if (g.gk == 2) {
    double x = kk * g.gd * g.gd / 4, t = 1, r = 1;
    for (int m = 1; m < 10; m++)
      t *= x / m, r += t;
    return exp(-x) * r;
  }
  if (g.gk == 3) {
    double x = sqrt(kk) * g.gd / g.gdk;
    int i = (int)x;
    if (i >= g.gn - 1)
      return g.gt[g.gn - 1];
    return g.gt[i] + (x - i) * (g.gt[i + 1] - g.gt[i]);
  }
  return 1;
}
static void cuda(cudaError_t e, const char *what) {
  if (e != cudaSuccess) {
    fprintf(stderr, "boxa: error: %s: %s\n", what, cudaGetErrorString(e));
    exit(1);
  }
}
static void cufft(cufftResult r, const char *what) {
  if (r != CUFFT_SUCCESS) {
    fprintf(stderr, "boxa: error: %s: cufft error %d\n", what, (int)r);
    exit(1);
  }
}
static void launched(const char *what) { cuda(cudaGetLastError(), what); }
static unsigned grid(long n) { return (n + nthread - 1) / nthread; }
__device__ static long wave(long i, long n) { return i <= n / 2 ? i : i - n; }
__device__ static void mode(G g, long l, long *i, long *j, long *k, double *kx,
                            double *ky, double *kz) {
  *i = l / (g.ny * g.nzf);
  *j = l / g.nzf % g.ny;
  *k = l % g.nzf;
  *kx = g.ax * wave(*i, g.nx);
  *ky = g.ay * wave(*j, g.ny);
  *kz = g.az * *k;
}
__device__ static double mask(G g, long i, long j, long k) {
  if (g.f) {
    double x = 2.0 * labs(wave(i, g.nx)) / g.nx, y = 2.0 * labs(wave(j, g.ny)) / g.ny,
           z = 2.0 * k / g.nz;
    return exp(-36 * (pow(x, 36) + pow(y, 36) + pow(z, 36)));
  }
  return 3 * labs(wave(i, g.nx)) < g.nx && 3 * labs(wave(j, g.ny)) < g.ny &&
         3 * k < g.nz;
}
__global__ static void spread(G g, int fl, const double2 *u, const double2 *v,
                              const double2 *w, double2 *b0, double2 *b1,
                              double2 *b2, double2 *b3, double2 *b4,
                              double2 *b5) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double2 a = u[l], b = v[l], c = w[l];
  double r = fl ? gu(g, kx * kx + ky * ky + kz * kz) : 1;
  b0[l] = make_double2(r * a.x, r * a.y);
  b1[l] = make_double2(r * b.x, r * b.y);
  b2[l] = make_double2(r * c.x, r * c.y);
  b3[l] = make_double2(-(ky * c.y - kz * b.y), ky * c.x - kz * b.x);
  b4[l] = make_double2(-(kz * a.y - kx * c.y), kz * a.x - kx * c.x);
  b5[l] = make_double2(-(kx * b.y - ky * a.y), kx * b.x - ky * a.x);
}
__global__ static void cross(G g, double s, double *b0, double *b1, double *b2,
                             const double *b3, const double *b4,
                             const double *b5) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= g.nr || l % (2 * g.nzf) >= g.nz)
    return;
  double u = b0[l] * s, v = b1[l] * s, w = b2[l] * s;
  double p = b3[l] * s, q = b4[l] * s, r = b5[l] * s;
  b0[l] = v * r - w * q;
  b1[l] = w * p - u * r;
  b2[l] = u * q - v * p;
}
/* Leray: b0..2 = u^ = G_n v^, b3..5 = v^ */
__global__ static void spreadL(G g, const double2 *u, const double2 *v, const double2 *w, double2 **b) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double r = gu(g, kx * kx + ky * ky + kz * kz);
  double2 a = u[l], c = v[l], e = w[l];
  b[0][l] = make_double2(r * a.x, r * a.y), b[1][l] = make_double2(r * c.x, r * c.y), b[2][l] = make_double2(r * e.x, r * e.y);
  b[3][l] = a, b[4][l] = c, b[5][l] = e;
}
/* c[3 j + i] = u_j v_i (physical, scaled) */
__global__ static void prod9(G g, double s, double **b, double **c) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= g.nr || l % (2 * g.nzf) >= g.nz)
    return;
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++)
      c[3 * j + i][l] = b[j][l] * s * b[3 + i][l] * s;
}
/* f_i = -i sum_j k_j (u_j v_i)^ */
__global__ static void ldiv(G g, double2 **c, double2 *f0, double2 *f1, double2 *f2) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double K[3] = {kx, ky, kz};
  double2 *f[3] = {f0, f1, f2};
  for (int q = 0; q < 3; q++) {
    double re = 0, im = 0;
    for (int d = 0; d < 3; d++)
      re += K[d] * c[3 * d + q][l].y, im -= K[d] * c[3 * d + q][l].x;
    f[q][l] = make_double2(re, im);
  }
}
/* strain of u: c0..c5 = S00 S11 S22 S01 S02 S12 (Fourier) */
__global__ static void sgrad(G g, const double2 *u, const double2 *v, const double2 *w, double2 **c) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double r = gu(g, kx * kx + ky * ky + kz * kz);
  double2 a = make_double2(r * u[l].x, r * u[l].y), b = make_double2(r * v[l].x, r * v[l].y),
          e = make_double2(r * w[l].x, r * w[l].y);
#define IK(K, z) make_double2(-(K) * (z).y, (K) * (z).x)
  c[0][l] = IK(kx, a), c[1][l] = IK(ky, b), c[2][l] = IK(kz, e);
  c[3][l] = make_double2(-0.5 * (ky * a.y + kx * b.y), 0.5 * (ky * a.x + kx * b.x));
  c[4][l] = make_double2(-0.5 * (kz * a.y + kx * e.y), 0.5 * (kz * a.x + kx * e.x));
  c[5][l] = make_double2(-0.5 * (kz * b.y + ky * e.y), 0.5 * (kz * b.x + ky * e.x));
#undef IK
}
/* c6 = C s(x), s = sqrt(2/3 S:S) (physical) */
__global__ static void smag(G g, double sc, double **c) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= g.nr || l % (2 * g.nzf) >= g.nz)
    return;
  double a = c[0][l] * sc, b = c[1][l] * sc, e = c[2][l] * sc, p = c[3][l] * sc, q = c[4][l] * sc, r = c[5][l] * sc;
  c[6][l] = g.wc * sqrt(2.0 / 3.0 * (a * a + b * b + e * e + 2 * (p * p + q * q + r * r)));
}
/* c0..2 = H^(1/2) v^ */
__global__ static void hhalf(G g, const double2 *u, const double2 *v, const double2 *w, double2 **c) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double h = hband(g, kx * kx + ky * ky + kz * kz);
  c[0][l] = make_double2(h * u[l].x, h * u[l].y), c[1][l] = make_double2(h * v[l].x, h * v[l].y),
  c[2][l] = make_double2(h * w[l].x, h * w[l].y);
}
/* c0..2 *= c6 (physical, scaled) */
__global__ static void smul(G g, double sc, double **c) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= g.nr || l % (2 * g.nzf) >= g.nz)
    return;
  for (int q = 0; q < 3; q++)
    c[q][l] *= c[6][l] * sc;
}
__global__ static void update(G g, int last, double nu, double dt, double a, double b,
                              const double2 *f0, const double2 *f1,
                              const double2 *f2, double2 *u, double2 *v,
                              double2 *w, double2 *du, double2 *dv,
                              double2 *dw, const double2 *F0, const double2 *F1, const double2 *F2,
                              const double2 *W0, const double2 *W1, const double2 *W2) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double kk = kx * kx + ky * ky + kz * kz, m = mask(g, i, j, k), visc = nu * kk;
  if (g.chi > 0)
    visc += g.chi * pow(fmax(0.0, 1 - filt(g, kk)), g.an + 1);
  double2 x = f0[l], y = f1[l], z = f2[l], p;
  if (W0 != NULL) {
    double hh = hband(g, kk);
    x.x -= hh * W0[l].x, x.y -= hh * W0[l].y;
    y.x -= hh * W1[l].x, y.y -= hh * W1[l].y;
    z.x -= hh * W2[l].x, z.y -= hh * W2[l].y;
  }
  if (F0 != NULL) {
    x.x += g.fa * F0[l].x, x.y += g.fa * F0[l].y;
    y.x += g.fa * F1[l].x, y.y += g.fa * F1[l].y;
    z.x += g.fa * F2[l].x, z.y += g.fa * F2[l].y;
  }
  x = make_double2(m * x.x, m * x.y);
  y = make_double2(m * y.x, m * y.y);
  z = make_double2(m * z.x, m * z.y);
  p = kk > 0 ? make_double2((kx * x.x + ky * y.x + kz * z.x) / kk,
                            (kx * x.y + ky * y.y + kz * z.y) / kk)
             : make_double2(0, 0);
  double2 U = u[l], V = v[l], W = w[l];
  double2 fx = make_double2(x.x - kx * p.x - visc * U.x,
                            x.y - kx * p.y - visc * U.y);
  double2 fy = make_double2(y.x - ky * p.x - visc * V.x,
                            y.y - ky * p.y - visc * V.y);
  double2 fz = make_double2(z.x - kz * p.x - visc * W.x,
                            z.y - kz * p.y - visc * W.y);
  double2 dx = du[l], dy = dv[l], dz = dw[l];
  dx = make_double2(a * dx.x + dt * fx.x, a * dx.y + dt * fx.y);
  dy = make_double2(a * dy.x + dt * fy.x, a * dy.y + dt * fy.y);
  dz = make_double2(a * dz.x + dt * fz.x, a * dz.y + dt * fz.y);
  du[l] = dx;
  dv[l] = dy;
  dw[l] = dz;
  double r = last && g.f ? (g.fp == 1 ? m : pow(m, g.fp)) : 1;
  u[l] = make_double2(r * (U.x + b * dx.x), r * (U.y + b * dx.y));
  v[l] = make_double2(r * (V.x + b * dy.x), r * (V.y + b * dy.y));
  w[l] = make_double2(r * (W.x + b * dz.x), r * (W.y + b * dz.y));
}
__device__ static void blockreduce(int K, int op, const double *v, double *s,
                                   double *out) {
  int i = threadIdx.x;
  for (int q = 0; q < K; q++)
    s[q * nthread + i] = v[q];
  __syncthreads();
  for (int h = nthread / 2; h > 0; h /= 2) {
    if (i < h)
      for (int q = 0; q < K; q++)
        s[q * nthread + i] = op ? fmax(s[q * nthread + i], s[q * nthread + i + h])
                                : s[q * nthread + i] + s[q * nthread + i + h];
    __syncthreads();
  }
  if (i == 0)
    for (int q = 0; q < K; q++)
      out[blockIdx.x * K + q] = s[q * nthread];
}
/* alpha norms (sums of h |v^|^2 weights): 0 E_v, 1 E_alpha (G), 2 E_u (G^2), 3 Omega_raw (k^2), 4 Omega_alpha (G k^2),
   5 Omega_u (G^2 k^2), 6 <F . u> (G Re F^* v^), 7 max-k tail of q (|k_i| > 0.9 k_c) */
__global__ static void anorms(G g, const double2 *u, const double2 *v, const double2 *w, const double2 *F0,
                              const double2 *F1, const double2 *F2, double *out) {
  __shared__ double s[8 * nthread];
  double e[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  double fx = g.f ? 0.4 * g.nx : (g.nx - 1) / 3, fy = g.f ? 0.4 * g.ny : (g.ny - 1) / 3,
         fz = g.f ? 0.4 * g.nz : (g.nz - 1) / 3;
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.nf; l += (long)gridDim.x * blockDim.x) {
    long i, j, k;
    double kx, ky, kz;
    mode(g, l, &i, &j, &k, &kx, &ky, &kz);
    double kk = kx * kx + ky * ky + kz * kz, h = k == 0 || 2 * k == g.nz ? 0.5 : 1, G = filt(g, kk);
    double2 a = u[l], b = v[l], c = w[l];
    double q = h * (a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y + c.x * c.x + c.y * c.y);
    e[0] += q, e[1] += G * q, e[2] += G * G * q, e[3] += kk * q, e[4] += G * kk * q, e[5] += G * G * kk * q;
    if (F0 != NULL)
      e[6] += 2 * h * G * (F0[l].x * a.x + F0[l].y * a.y + F1[l].x * b.x + F1[l].y * b.y + F2[l].x * c.x + F2[l].y * c.y);
    if (labs(wave(i, g.nx)) > 0.9 * fx || labs(wave(j, g.ny)) > 0.9 * fy || k > 0.9 * fz)
      e[7] += kk * q;
  }
  blockreduce(8, 0, e, s, out);
}
/* shell spectra (integer shells of |k| in units of the smallest wavenumber): sp[3 n + 0..2] += |v^|^2, G |v^|^2, G^2 |v^|^2 */
__global__ static void kspec(G g, const double2 *u, const double2 *v, const double2 *w, int ns, double *sp) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double kk = kx * kx + ky * ky + kz * kz, h = k == 0 || 2 * k == g.nz ? 0.5 : 1, G = filt(g, kk);
  int n = (int)lround(sqrt(kk) / fmin(g.ax, fmin(g.ay, g.az)));
  if (n >= ns)
    return;
  double2 a = u[l], b = v[l], c = w[l];
  double q = 2 * h * (a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y + c.x * c.x + c.y * c.y);
  atomicAdd(sp + 3 * n, q);
  atomicAdd(sp + 3 * n + 1, G * q);
  atomicAdd(sp + 3 * n + 2, G * G * q);
}
/* random low-k solenoidal field (|k| <= 4), amplitude a (spectral, before scaling) */
__device__ static double hrand(unsigned long x) {
  x ^= x >> 33, x *= 0xff51afd7ed558ccdUL, x ^= x >> 33, x *= 0xc4ceb9fe1a85ec53UL, x ^= x >> 33;
  return (x >> 11) * (1.0 / 9007199254740992.0) - 0.5;
}
__global__ static void knoise(G g, double a, double2 *u, double2 *v, double2 *w) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double kk = kx * kx + ky * ky + kz * kz;
  if (kk == 0 || kk > 16.5 || k == 0) {
    u[l] = v[l] = w[l] = make_double2(0, 0);
    return;
  }
  double r[6];
  for (int q = 0; q < 6; q++)
    r[q] = a * hrand(6 * l + q + 12345);
  double p0 = (kx * r[0] + ky * r[2] + kz * r[4]) / kk, p1 = (kx * r[1] + ky * r[3] + kz * r[5]) / kk;
  u[l] = make_double2(r[0] - kx * p0, r[1] - kx * p1);
  v[l] = make_double2(r[2] - ky * p0, r[3] - ky * p1);
  w[l] = make_double2(r[4] - kz * p0, r[5] - kz * p1);
}
__global__ static void norms(G g, const double2 *u, const double2 *v,
                             const double2 *w, double *out) {
  __shared__ double s[6 * nthread];
  double e[6] = {0, 0, 0, 0, 0, 0};
  double fx = g.f ? 0.4 * g.nx : (g.nx - 1) / 3, fy = g.f ? 0.4 * g.ny : (g.ny - 1) / 3,
         fz = g.f ? 0.4 * g.nz : (g.nz - 1) / 3;
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.nf;
       l += (long)gridDim.x * blockDim.x) {
    long i, j, k;
    double kx, ky, kz;
    mode(g, l, &i, &j, &k, &kx, &ky, &kz);
    double kk = kx * kx + ky * ky + kz * kz,
           h = k == 0 || 2 * k == g.nz ? 0.5 : 1;
    double2 a = u[l], b = v[l], c = w[l];
    double q = h * (a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y + c.x * c.x +
                    c.y * c.y);
    e[0] += q;
    e[1] += kk * q;
    if (kk > 0)
      e[2] += q / sqrt(kk);
    if (labs(wave(i, g.nx)) > 0.9 * fx)
      e[3] += q;
    if (labs(wave(j, g.ny)) > 0.9 * fy)
      e[4] += q;
    if (k > 0.9 * fz)
      e[5] += q;
  }
  blockreduce(6, 0, e, s, out);
}
__global__ static void peaks(G g, double s, double cx, double cy, double cz,
                             const double *b0, const double *b1,
                             const double *b2, const double *b3,
                             const double *b4, const double *b5, double *out) {
  __shared__ double sh[2 * nthread];
  double m[2] = {0, 0};
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.nr;
       l += (long)gridDim.x * blockDim.x) {
    if (l % (2 * g.nzf) >= g.nz)
      continue;
    double c = s * (cx * fabs(b0[l]) + cy * fabs(b1[l]) + cz * fabs(b2[l]));
    double o = s * s * (b3[l] * b3[l] + b4[l] * b4[l] + b5[l] * b5[l]);
    m[0] = fmax(m[0], c);
    m[1] = fmax(m[1], o);
  }
  blockreduce(2, 1, m, sh, out);
}
static void collect(int K, int op, const double *dP, double *P, double *r) {
  cuda(cudaMemcpy(P, dP, nblock * K * sizeof(double), cudaMemcpyDeviceToHost),
       "cudaMemcpy");
  for (int q = 0; q < K; q++)
    r[q] = 0;
  for (int b = 0; b < nblock; b++)
    for (int q = 0; q < K; q++)
      r[q] = op ? fmax(r[q], P[b * K + q]) : r[q] + P[b * K + q];
}
static long steps(const char *opt, double x, double dt) {
  long n = lround(x / dt);
  if (fabs(n * dt - x) > 1e-9 * x) {
    fprintf(stderr, "boxa: error: %s %g is not a multiple of -s %g\n", opt, x,
            dt);
    exit(1);
  }
  return n;
}
int main(int argc, char **argv) {
  G g;
  cufftHandle fplan, bplan;
  FILE *file;
  char path[FILENAME_MAX], *ipath, *rpath, *end, *name;
  double nu, dt, T, t, dd, F, X, Y, Z, x, N, s, cx, cy, cz, kc, sums[6], mx[2],
      *host, *dP, *P;
  double2 *u[3], *du[3], *b[6], *Fh[3] = {NULL, NULL, NULL};
  char *tpath = NULL;
  double k0 = 0, E0 = 0.5, tau = 0.5, Ramp = 1e-4, dts = 0, fmx = 1, as[8], *dsp = NULL, *hsp = NULL;
  int ns = 0;
  long nS = 0;
  g.gk = 0, g.gd = 0, g.gt = NULL, g.gn = 0, g.gdk = 1, g.fa = 0;
  g.an = 0, g.chi = 0, g.wc = 0, g.wp = 0.5;
  int leray = 0;
  double2 *c[9] = {NULL}, *cw[3] = {NULL}, **dc = NULL;
  long tstep, start, nd, nx, ny, nz, nt;
  size_t wf, wb;
  void *work;
  (void)argc;
  ipath = rpath = NULL;
  nx = ny = nz = 0;
  X = Y = Z = 2;
  nu = dt = -1;
  T = dd = 0;
  g.f = 0;
  F = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    if (argv[0][1] == 'h') {
      fprintf(stderr,
              "Usage: box -x <nx> -y <ny> -z <nz> [-X <Lx/pi>] [-Y <Ly/pi>] "
              "[-Z <Lz/pi>] (-i <u.raw> | -r <dump u.STEP>) -n <viscosity> "
              "-t <end time> -s <time step> [-d <dump interval>] [-f (Hou-Li "
              "filter instead of 2/3)] [-F <dt0> (filter "
              "the solution as m^(dt/dt0) per step: same damping per unit time as -s dt0)]\n");
      exit(1);
    }
    if (argv[0][1] == 'f') {
      g.f = 1;
      continue;
    }
    if (argv[1] == NULL) {
      fprintf(stderr, "boxa: error: %s needs an argument\n", argv[0]);
      exit(1);
    }
    if (argv[0][1] == 'i') {
      ipath = *++argv;
      continue;
    }
    if (argv[0][1] == 'r') {
      rpath = *++argv;
      continue;
    }
    if (argv[0][1] == 'T') {
      tpath = *++argv;
      continue;
    }
    x = strtod(argv[1], &end);
    if (*end != '\0') {
      fprintf(stderr, "boxa: error: '%s' is not a number\n", argv[1]);
      exit(1);
    }
    switch (argv[0][1]) {
    case 'x':
      nx = x;
      break;
    case 'y':
      ny = x;
      break;
    case 'z':
      nz = x;
      break;
    case 'X':
      X = x;
      break;
    case 'Y':
      Y = x;
      break;
    case 'Z':
      Z = x;
      break;
    case 'n':
      nu = x;
      break;
    case 's':
      dt = x;
      break;
    case 't':
      T = x;
      break;
    case 'd':
      dd = x;
      break;
    case 'F':
      F = x;
      break;
    case 'G':
      g.gk = (int)x;
      break;
    case 'D':
      g.gd = x;
      break;
    case 'K':
      k0 = x;
      break;
    case 'e':
      E0 = x;
      break;
    case 'c':
      tau = x;
      break;
    case 'R':
      Ramp = x;
      break;
    case 'S':
      dts = x;
      break;
    case 'M':
      fmx = x;
      break;
    case 'L':
      leray = (int)x;
      break;
    case 'A':
      g.an = (int)x;
      break;
    case 'B':
      g.chi = x;
      break;
    case 'W':
      g.wc = x;
      break;
    case 'P':
      g.wp = x;
      break;
    default:
      fprintf(stderr, "boxa: error: unknown option '%s'\n", *argv);
      exit(1);
    }
    argv++;
  }
  if (nx < 4 || ny < 4 || nz < 4 || nz % 2 != 0 || nu < 0 || dt <= 0 ||
      T <= 0 || X <= 0 || Y <= 0 || Z <= 0 || (ipath == NULL) == (rpath == NULL)) {
    fprintf(stderr, "boxa: error: need -x -y -z (nz even), -n, -s, -t and one "
                    "of -i, -r\n");
    exit(1);
  }
  g.nx = nx;
  g.ny = ny;
  g.nz = nz;
  g.nzf = nz / 2 + 1;
  g.fp = F > 0 ? dt / F : 1;
  g.nf = nx * ny * g.nzf;
  g.nr = 2 * g.nf;
  g.ax = 2 / X;
  g.ay = 2 / Y;
  g.az = 2 / Z;
  N = (double)nx * ny * nz;
  s = 1 / N;
  cx = g.ax * ((nx - 1) / 3);
  cy = g.ay * ((ny - 1) / 3);
  cz = g.az * ((nz - 1) / 3);
  kc = fmin(cx, fmin(cy, cz));
  nd = dd > 0 ? steps("-d", dd, dt) : 0;
  nt = steps("-t", T, dt);
  for (int c = 0; c < 3; c++) {
    cuda(cudaMalloc(&u[c], g.nf * sizeof(double2)), "cudaMalloc");
    cuda(cudaMalloc(&du[c], g.nf * sizeof(double2)), "cudaMalloc");
    cuda(cudaMemset(du[c], 0, g.nf * sizeof(double2)), "cudaMemset");
  }
  for (int c = 0; c < 6; c++)
    cuda(cudaMalloc(&b[c], g.nf * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&dP, 8 * nblock * sizeof(double)), "cudaMalloc");
  double2 **db = NULL;
  if (leray || g.wc > 0) {
    for (int q = 0; q < 9; q++)
      cuda(cudaMalloc(&c[q], g.nf * sizeof(double2)), "cudaMalloc");
    if (g.wc > 0)
      for (int q = 0; q < 3; q++)
        cuda(cudaMalloc(&cw[q], g.nf * sizeof(double2)), "cudaMalloc");
    cuda(cudaMalloc(&dc, 9 * sizeof(double2 *)), "cudaMalloc");
    cuda(cudaMemcpy(dc, c, 9 * sizeof(double2 *), cudaMemcpyHostToDevice), "cudaMemcpy");
    cuda(cudaMalloc(&db, 6 * sizeof(double2 *)), "cudaMalloc");
    cuda(cudaMemcpy(db, b, 6 * sizeof(double2 *), cudaMemcpyHostToDevice), "cudaMemcpy");
  }
  fprintf(stderr, "boxa: closure leray %d adm n %d chi %g strain sink C %g p %g\n", leray, g.an, g.chi, g.wc, g.wp);
  if (g.gk == 3) {
    /* radial table: lines "k delta, G" on a uniform grid starting at 0 */
    FILE *tf = tpath ? fopen(tpath, "r") : NULL;
    double a0, a1, *ht = (double *)malloc(1000000 * sizeof(double)), k1 = 0;
    if (tf == NULL) {
      fprintf(stderr, "boxa: error: -G 3 needs -T table\n");
      exit(1);
    }
    while (g.gn < 1000000 && fscanf(tf, "%lf %lf", &a0, &a1) == 2) {
      if (g.gn == 1)
        k1 = a0;
      ht[g.gn++] = a1;
    }
    fclose(tf);
    g.gdk = k1;
    double *dt_;
    cuda(cudaMalloc(&dt_, g.gn * sizeof(double)), "cudaMalloc");
    cuda(cudaMemcpy(dt_, ht, g.gn * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    g.gt = dt_;
    free(ht);
    fprintf(stderr, "boxa: table %s: %d entries, step %g in k delta\n", tpath, g.gn, g.gdk);
  }
  if (dts > 0) {
    nS = steps("-S", dts, dt);
    ns = (int)(sqrt((double)(nx * nx + ny * ny + nz * nz)) / 2) + 2;
    cuda(cudaMalloc(&dsp, 3 * ns * sizeof(double)), "cudaMalloc");
    hsp = (double *)malloc(3 * ns * sizeof(double));
  }
  if ((P = (double *)malloc(8 * nblock * sizeof(double))) == NULL ||
      (host = (double *)malloc(nx * ny * nz * sizeof(double))) == NULL) {
    fprintf(stderr, "boxa: error: malloc failed\n");
    exit(1);
  }
  cufft(cufftCreate(&fplan), "cufftCreate");
  cufft(cufftCreate(&bplan), "cufftCreate");
  cufft(cufftSetAutoAllocation(fplan, 0), "cufftSetAutoAllocation");
  cufft(cufftSetAutoAllocation(bplan, 0), "cufftSetAutoAllocation");
  cufft(cufftMakePlan3d(fplan, nx, ny, nz, CUFFT_D2Z, &wf), "cufftMakePlan3d");
  cufft(cufftMakePlan3d(bplan, nx, ny, nz, CUFFT_Z2D, &wb), "cufftMakePlan3d");
  cuda(cudaMalloc(&work, wf > wb ? wf : wb), "cudaMalloc");
  cufft(cufftSetWorkArea(fplan, work), "cufftSetWorkArea");
  cufft(cufftSetWorkArea(bplan, work), "cufftSetWorkArea");
  tstep = 0;
  name = ipath;
  if (k0 > 0) {
    /* Taylor-Green forcing F = (sin k0 x cos k0 y cos k0 z, -cos k0 x sin k0 y cos k0 z, 0) */
    double Lx = X * pi, Ly = Y * pi, Lz = Z * pi;
    for (int c = 0; c < 3; c++) {
      for (long i = 0; i < nx; i++)
        for (long j = 0; j < ny; j++)
          for (long k = 0; k < nz; k++) {
            double x = i * Lx / nx, y = j * Ly / ny, z = k * Lz / nz, f;
            f = c == 0 ? sin(k0 * x) * cos(k0 * y) * cos(k0 * z) : c == 1 ? -cos(k0 * x) * sin(k0 * y) * cos(k0 * z) : 0;
            host[(i * ny + j) * nz + k] = f;
          }
      cuda(cudaMalloc(&Fh[c], g.nf * sizeof(double2)), "cudaMalloc");
      cuda(cudaMemcpy2D(b[0], 2 * g.nzf * sizeof(double), host, nz * sizeof(double), nz * sizeof(double), nx * ny,
                        cudaMemcpyHostToDevice), "cudaMemcpy2D");
      cufft(cufftExecD2Z(fplan, (double *)b[0], b[0]), "cufftExecD2Z");
      cuda(cudaMemcpy(Fh[c], b[0], g.nf * sizeof(double2), cudaMemcpyDeviceToDevice), "cudaMemcpy");
    }
  }
  if (rpath == NULL && strcmp(ipath, "tgv") == 0) {
    /* Taylor-Green vortex u = (sin x cos y cos z, -cos x sin y cos z, 0) */
    double Lx = X * pi, Ly = Y * pi, Lz = Z * pi;
    for (int c = 0; c < 3; c++) {
      for (long i = 0; i < nx; i++)
        for (long j = 0; j < ny; j++)
          for (long k = 0; k < nz; k++) {
            double x = i * Lx / nx, y = j * Ly / ny, z = k * Lz / nz;
            host[(i * ny + j) * nz + k] = c == 0 ? sin(x) * cos(y) * cos(z) : c == 1 ? -cos(x) * sin(y) * cos(z) : 0;
          }
      cuda(cudaMemcpy2D(b[0], 2 * g.nzf * sizeof(double), host, nz * sizeof(double), nz * sizeof(double), nx * ny,
                        cudaMemcpyHostToDevice), "cudaMemcpy2D");
      cufft(cufftExecD2Z(fplan, (double *)b[0], b[0]), "cufftExecD2Z");
      cuda(cudaMemcpy(u[c], b[0], g.nf * sizeof(double2), cudaMemcpyDeviceToDevice), "cudaMemcpy");
    }
    name = NULL;
  }
  if (rpath == NULL && strcmp(ipath, "zero") == 0) {
    knoise<<<grid(g.nf), nthread>>>(g, Ramp * N, u[0], u[1], u[2]);
    launched("knoise");
    name = NULL;
  }
  if (name != NULL) {
  if (rpath != NULL) {
    const char *p = strrchr(rpath, '.');
    if (p == NULL || (tstep = strtol(p + 1, &end, 10), *end != '\0' || end == p + 1)) {
      fprintf(stderr, "boxa: error: '%s' is not named u.STEP\n", rpath);
      exit(1);
    }
    name = rpath;
  }
  if ((file = fopen(name, "r")) == NULL) {
    fprintf(stderr, "boxa: error: fail to open '%s'\n", name);
    exit(1);
  }
  for (int c = 0; c < 3; c++) {
    if (fread(host, sizeof(double), nx * ny * nz, file) != (size_t)(nx * ny * nz)) {
      fprintf(stderr, "boxa: error: fail to read '%s' (wrong -x -y -z?)\n", name);
      exit(1);
    }
    cuda(cudaMemcpy2D(b[0], 2 * g.nzf * sizeof(double), host, nz * sizeof(double),
                      nz * sizeof(double), nx * ny, cudaMemcpyHostToDevice),
         "cudaMemcpy2D");
    cufft(cufftExecD2Z(fplan, (double *)b[0], b[0]), "cufftExecD2Z");
    cuda(cudaMemcpy(u[c], b[0], g.nf * sizeof(double2), cudaMemcpyDeviceToDevice),
         "cudaMemcpy");
  }
  if (fgetc(file) != EOF) {
    fprintf(stderr, "boxa: error: '%s' is too long (wrong -x -y -z?)\n", name);
    exit(1);
  }
  fclose(file);
  }
  start = tstep;
  t = tstep * dt;
  for (;;) {
    int diag = tstep % 10 == 0, dump = nd > 0 && tstep % nd == 0 && (rpath == NULL || tstep != start);
    if (diag || dump) {
      spread<<<grid(g.nf), nthread>>>(g, 0, u[0], u[1], u[2], b[0], b[1], b[2], b[3],
                                      b[4], b[5]);
      launched("spread");
      for (int c = 0; c < 6; c++)
        cufft(cufftExecZ2D(bplan, b[c], (double *)b[c]), "cufftExecZ2D");
    }
    if (nS > 0 && tstep % nS == 0) {
      cuda(cudaMemset(dsp, 0, 3 * ns * sizeof(double)), "cudaMemset");
      kspec<<<grid(g.nf), nthread>>>(g, u[0], u[1], u[2], ns, dsp);
      launched("kspec");
      cuda(cudaMemcpy(hsp, dsp, 3 * ns * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
      printf("#s %.10e", t);
      for (int n = 0; n < ns; n++)
        printf(" %.6e %.6e %.6e", 0.5 * hsp[3 * n] * s * s, 0.5 * hsp[3 * n + 1] * s * s, 0.5 * hsp[3 * n + 2] * s * s);
      printf("\n");
    }
    if (k0 > 0 || diag) {
      anorms<<<nblock, nthread>>>(g, u[0], u[1], u[2], Fh[0], Fh[1], Fh[2], dP);
      launched("anorms");
      collect(8, 0, dP, P, as);
      for (int q = 0; q < 8; q++)
        as[q] *= s * s;
      if (k0 > 0) {
        /* energy control: injection fa <F.u> = 2 nu Omega_alpha + (E0 - E_alpha) / tau */
        double D = 2 * nu * as[4], I = as[6], fa = I > 0 ? (D + (E0 - as[1]) / tau) / I : fmx;
        g.fa = fmax(0, fmin(fa, fmx));
      }
    }
    if (diag) {
      norms<<<nblock, nthread>>>(g, u[0], u[1], u[2], dP);
      launched("norms");
      collect(6, 0, dP, P, sums);
      peaks<<<nblock, nthread>>>(g, s, cx, cy, cz, (double *)b[0], (double *)b[1],
                                 (double *)b[2], (double *)b[3], (double *)b[4],
                                 (double *)b[5], dP);
      launched("peaks");
      collect(2, 1, dP, P, mx);
      double E = sums[0] * s * s, Om = sums[1] * s * s, u2 = 2 * E / 3;
      double q[] = {t, as[0], as[1], as[2], as[3], as[4], as[5], 2 * nu * as[4], g.fa, g.fa * as[6], dt * mx[0], sqrt(mx[1]),
                    as[7] / as[3], sqrt(2 * as[0]) / (nu > 0 ? nu : 1), 2 * pi * sqrt(as[0] / as[3]),
                    pi / (2 * u2) * sums[2] * s * s};
      (void)Om;
      if (tstep == start)
        printf("step t Ev Ealpha Eu Omraw Omalpha Omu eps_alpha fa inj C wmax tail Re lambda L\n");
      printf("% 10ld", tstep);
      for (int i = 0; i < (int)(sizeof q / sizeof *q); i++)
        printf(" % .16e", q[i]);
      printf("\n");
      if (fflush(stdout) != 0 || ferror(stdout)) {
        fprintf(stderr, "boxa: error: fail to write stdout\n");
        exit(1);
      }
    }
    if (dump) {
      sprintf(path, "u.%08ld", tstep);
      if ((file = fopen(path, "w")) == NULL) {
        fprintf(stderr, "boxa: error: fail to open '%s'\n", path);
        exit(1);
      }
      for (int c = 0; c < 3; c++) {
        cuda(cudaMemcpy2D(host, nz * sizeof(double), b[c],
                          2 * g.nzf * sizeof(double), nz * sizeof(double),
                          nx * ny, cudaMemcpyDeviceToHost),
             "cudaMemcpy2D");
        for (long l = 0; l < nx * ny * nz; l++)
          host[l] *= s;
        if (fwrite(host, sizeof(double), nx * ny * nz, file) !=
            (size_t)(nx * ny * nz)) {
          fprintf(stderr, "boxa: error: fail to write '%s'\n", path);
          exit(1);
        }
      }
      if (ferror(file) || fclose(file) != 0) {
        fprintf(stderr, "boxa: error: fail to write '%s'\n", path);
        exit(1);
      }
    }
    if (tstep >= nt)
      break;
    for (int r = 0; r < 5; r++) {
      double2 *W[3] = {NULL, NULL, NULL};
      if (g.wc > 0) {
        /* strain-scaled band sink term, Fourier, in c0..c2 (before c is reused by Leray) */
        sgrad<<<grid(g.nf), nthread>>>(g, u[0], u[1], u[2], dc);
        launched("sgrad");
        for (int q = 0; q < 6; q++)
          cufft(cufftExecZ2D(bplan, c[q], (double *)c[q]), "cufftExecZ2D");
        smag<<<grid(g.nr), nthread>>>(g, s, (double **)dc);
        launched("smag");
        hhalf<<<grid(g.nf), nthread>>>(g, u[0], u[1], u[2], dc);
        launched("hhalf");
        for (int q = 0; q < 3; q++)
          cufft(cufftExecZ2D(bplan, c[q], (double *)c[q]), "cufftExecZ2D");
        smul<<<grid(g.nr), nthread>>>(g, s, (double **)dc);
        launched("smul");
        for (int q = 0; q < 3; q++) {
          cufft(cufftExecD2Z(fplan, (double *)c[q], c[q]), "cufftExecD2Z");
          cuda(cudaMemcpy(cw[q], c[q], g.nf * sizeof(double2), cudaMemcpyDeviceToDevice), "cudaMemcpy");
          W[q] = cw[q];
        }
      }
      if (leray) {
        /* N = -P div(u v) */
        spreadL<<<grid(g.nf), nthread>>>(g, u[0], u[1], u[2], db);
        launched("spreadL");
        for (int q = 0; q < 6; q++)
          cufft(cufftExecZ2D(bplan, b[q], (double *)b[q]), "cufftExecZ2D");
        prod9<<<grid(g.nr), nthread>>>(g, s, (double **)db, (double **)dc);
        launched("prod9");
        for (int q = 0; q < 9; q++)
          cufft(cufftExecD2Z(fplan, (double *)c[q], c[q]), "cufftExecD2Z");
        ldiv<<<grid(g.nf), nthread>>>(g, dc, b[0], b[1], b[2]);
        launched("ldiv");
      } else {
        spread<<<grid(g.nf), nthread>>>(g, 1, u[0], u[1], u[2], b[0], b[1], b[2], b[3],
                                        b[4], b[5]);
        launched("spread");
        for (int c = 0; c < 6; c++)
          cufft(cufftExecZ2D(bplan, b[c], (double *)b[c]), "cufftExecZ2D");
        cross<<<grid(g.nr), nthread>>>(g, s, (double *)b[0], (double *)b[1],
                                       (double *)b[2], (double *)b[3],
                                       (double *)b[4], (double *)b[5]);
        launched("cross");
        for (int c = 0; c < 3; c++)
          cufft(cufftExecD2Z(fplan, (double *)b[c], b[c]), "cufftExecD2Z");
      }
      update<<<grid(g.nf), nthread>>>(g, r == 4, nu, dt, A[r], B[r], b[0], b[1], b[2],
                                      u[0], u[1], u[2], du[0], du[1], du[2], Fh[0], Fh[1], Fh[2], W[0], W[1], W[2]);
      launched("update");
    }
    tstep++;
    t = tstep * dt;
  }
  cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
}
