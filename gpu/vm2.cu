#include <cufft.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>
#include <thrust/inner_product.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <cub/cub.cuh>

enum { nthread = 256, B = 16, NB = B * B, PM = 40, DM = 14, nw = 256 };
static const double pi = 3.141592653589793238;
static const double K8[4][8] = {
    {1, 0, -49. / 36, 0, -959. / 144, 2569. / 144, -727. / 48, 623. / 144},
    {0, -3. / 4, 3. / 4, 13. / 48, 95. / 24, -429. / 40, 1091. / 120, -623. / 240},
    {0, 3. / 20, -3. / 40, -1. / 6, -191. / 144, 2581. / 720, -2183. / 720, 623. / 720},
    {0, -1. / 60, 1. / 180, 1. / 48, 7. / 36, -37. / 72, 13. / 30, -89. / 720}};
__constant__ double cK[4][8], cR[8], cE[8];
__constant__ int nR, nE;
struct T {
  long N;
  int p, Pt, nleaf, nbox, Dm, Fn, mode;
  double L, cd, nu;
  int *fine, *key, *start, *lbox, *lev, *bxx, *byy;
  long *perm;
  double *x, *y, *g, *h, *v, *o, *ox, *oy, *og, *oh, *ov, *oo, *cx, *cy, *dg, *odg;
  int *f, *of;
  double2 *M, *Lc, *Tm2m, *Tl2l, *Tm2l, *Tfar, *w, *ow;
  int *uo, *ul, *us, *vo, *vl, *vs, *wo, *wl, *ws, *xo, *xl, *xs, *ro, *rl, *rs, *co, *cl, *cs;
  int *bdep, *bch, *bpx, *bpy, *bpar;
};
static void ck(cudaError_t e, const char *what) {
  if (e != cudaSuccess) {
    fprintf(stderr, "vm2: error: %s: %s\n", what, cudaGetErrorString(e));
    exit(1);
  }
}
static void cufft(cufftResult r, const char *what) {
  if (r != CUFFT_SUCCESS) {
    fprintf(stderr, "vm2: error: %s: cufft error %d\n", what, (int)r);
    exit(1);
  }
}
static void launched(const char *what) { ck(cudaGetLastError(), what); }
static unsigned grid(long n) { return (n + nthread - 1) / nthread; }
static void *dmalloc(size_t n) {
  void *a;
  ck(cudaMalloc(&a, n > 0 ? n : 1), "cudaMalloc");
  return a;
}
static void *hmalloc(size_t n) {
  void *a = malloc(n > 0 ? n : 1);
  if (a == NULL) {
    fprintf(stderr, "vm2: error: malloc failed\n");
    exit(1);
  }
  return a;
}
static void up(void *d, const void *h, size_t n) {
  ck(cudaMemcpy(d, h, n, cudaMemcpyHostToDevice), "cudaMemcpy");
}
static void down(void *h, const void *d, size_t n) {
  ck(cudaMemcpy(h, d, n, cudaMemcpyDeviceToHost), "cudaMemcpy");
}
__host__ __device__ static double2 mk(double a, double b) {
  double2 r;
  r.x = a;
  r.y = b;
  return r;
}
__host__ __device__ static double2 add(double2 a, double2 b) {
  return mk(a.x + b.x, a.y + b.y);
}
__host__ __device__ static double2 sub(double2 a, double2 b) {
  return mk(a.x - b.x, a.y - b.y);
}
__host__ __device__ static double2 mul(double2 a, double2 b) {
  return mk(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
__host__ __device__ static double2 scl(double2 a, double s) {
  return mk(a.x * s, a.y * s);
}
__host__ __device__ static double2 inv(double2 a) {
  double d = a.x * a.x + a.y * a.y;
  return mk(a.x / d, -a.y / d);
}
__host__ __device__ static double2 lg(double2 a) {
  return mk(log(hypot(a.x, a.y)), atan2(a.y, a.x));
}
__device__ static double2 sub2(double2 a, double2 b) { return mk(a.x - b.x, a.y - b.y); }
__device__ static double kw(double y) {
  double a = fabs(y);
  int j = (int)a;
  if (j >= 4)
    return 0;
  double u = a - j, w = 0;
  for (int e = 7; e >= 0; e--)
    w = w * u + cK[j][e];
  return w;
}
__device__ static double poly(const double *c, int n, double x) {
  double s = 0;
  for (int e = n - 1; e >= 0; e--)
    s = s * x + c[e];
  return s;
}
__device__ static double kcore(double x, double ex) {
  double s = 0;
  for (int e = nR - 1; e >= 1; e--)
    s = s * x + cR[e];
  return -expm1(-x) - ex * x * s;
}
__device__ static double wrap(double x, double L) {
  x -= L * floor(x / L);
  return x < L ? x : 0;
}
__global__ static void leafof(T t) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i >= t.N)
    return;
  double x = wrap(t.ox[i], t.L), y = wrap(t.oy[i], t.L);
  t.ox[i] = x;
  t.oy[i] = y;
  int a = (int)(x / t.L * t.Fn), b = (int)(y / t.L * t.Fn);
  a = a < t.Fn ? a : t.Fn - 1;
  b = b < t.Fn ? b : t.Fn - 1;
  t.key[i] = t.fine[(long)b * t.Fn + a];
}
__global__ static void gather(T t) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i >= t.N)
    return;
  long j = t.perm[i];
  t.x[i] = t.ox[j];
  t.y[i] = t.oy[j];
  t.g[i] = t.og[j];
  t.o[i] = t.oo[j];
  t.h[i] = t.oh[j];
  t.v[i] = t.ov[j];
  t.f[i] = t.of[j];
}
__global__ static void p2m(T t) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int P = t.p + 1;
  if (j >= (long)t.nleaf * P)
    return;
  int c = j / P, k = j % P, b = t.lbox[c];
  double sc = t.L / (double)(1L << t.bdep[b]);
  double2 s = mk(0, 0);
  for (int i = t.start[c]; i < t.start[c + 1]; i++) {
    double2 q = mk(0, -t.g[i] / (2 * pi)), d = mk((t.x[i] - t.cx[b]) / sc, (t.y[i] - t.cy[b]) / sc),
            pw = mk(1, 0);
    for (int m = 0; m < k; m++)
      pw = mul(pw, d);
    s = add(s, mul(q, pw));
  }
  t.M[(long)b * P + k] = k == 0 ? s : scl(s, -1.0 / k);
}
__global__ static void m2m(T t, int lo, int hi, int d) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int P = t.p + 1;
  if (j >= (long)(hi - lo) * P)
    return;
  int b = lo + j / P, r = j % P;
  if (t.bch[b] < 0)
    return;
  double2 s = mk(0, 0);
  for (int q = 0; q < 4; q++) {
    const double2 *A = t.Tm2m + (long)q * t.Pt * t.Pt + r * t.Pt;
    const double2 *m = t.M + (long)(t.bch[b] + q) * P;
    (void)d;
    for (int k = 0; k < P; k++)
      s = add(s, mul(A[k], m[k]));
  }
  t.M[(long)b * P + r] = s;
}
__global__ static void far(T t) {
  int r = threadIdx.x, P = t.p + 1;
  if (r >= P)
    return;
  double2 s = mk(0, 0);
  for (int k = 0; k < P; k++)
    s = add(s, mul(t.Tfar[r * t.Pt + k], t.M[k]));
  t.Lc[r] = s;
}
__global__ static void m2l(T t) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int P = t.p + 1;
  if (j >= (long)t.nbox * P)
    return;
  int b = j / P, r = j % P, d = t.bdep[b], n = 1 << d;
  double2 s = mk(0, 0);
  for (int e = t.vo[b]; e < t.vo[b + 1]; e++) {
    int c = t.vl[e], sc = t.vs[e];
    int ox = t.bpx[c] + (sc / 3 - 1) * n - t.bpx[b],
        oy = t.bpy[c] + (sc % 3 - 1) * n - t.bpy[b];
    const double2 *A = t.Tm2l + (long)((ox + 3) * 7 + oy + 3) * t.Pt * t.Pt + r * t.Pt;
    const double2 *m = t.M + (long)c * P;
    for (int k = 0; k < P; k++)
      s = add(s, mul(A[k], m[k]));
  }
  if (b > 0)
    t.Lc[j] = s;
}
__global__ static void l2l(T t, int lo, int hi, int d) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int P = t.p + 1;
  if (j >= (long)(hi - lo) * P)
    return;
  int b = lo + j / P, r = j % P;
  int q = (t.bpx[b] & 1) * 2 + (t.bpy[b] & 1);
  const double2 *A = t.Tl2l + (long)q * t.Pt * t.Pt + r * t.Pt;
  const double2 *m = t.Lc + (long)t.bpar[b] * P;
  (void)d;
  double2 u = mk(0, 0);
  for (int k = 0; k < P; k++)
    u = add(u, mul(A[k], m[k]));
  t.Lc[(long)b * P + r] = add(t.Lc[(long)b * P + r], u);
}
__global__ static void l2p(T t) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i >= t.N)
    return;
  int c = t.key[i], b = t.lbox[c], P = t.p + 1;
  if (t.mode == 2) {
    t.w[i] = mk(0, 0);
    return;
  }
  const double2 *a = t.Lc + (long)b * P;
  double sb = t.L / (double)(1L << t.bdep[b]);
  double2 d = mk((t.x[i] - t.cx[b]) / sb, (t.y[i] - t.cy[b]) / sb), w = mk(0, 0);
  for (int l = t.p; l > 0; l--)
    w = add(mul(w, d), scl(a[l], l));
  w = scl(w, 1 / sb);
  for (int e = t.wo[c]; e < t.wo[c + 1]; e++) {
    int s = t.wl[e], sc = t.ws[e];
    const double2 *m = t.M + (long)s * P;
    double ss = t.L / (double)(1L << t.bdep[s]);
    double2 z = mk(t.x[i] - t.cx[s] - (sc / 3 - 1) * t.L, t.y[i] - t.cy[s] - (sc % 3 - 1) * t.L),
            iz = inv(z), u = scl(iz, ss), acc = mk(0, 0);
    for (int k = t.p; k > 0; k--)
      acc = mul(add(acc, scl(m[k], k)), u);
    w = add(w, mul(sub(m[0], acc), iz));
  }
  t.w[i] = w;
}
__global__ static void p2p(T t) {
  __shared__ double sx[nw], sy[nw], sg[nw], sv[nw];
  __shared__ int sf[nw];
  int c = blockIdx.x;
  if (c >= t.nleaf)
    return;
  int n0 = t.start[c], cnt = t.start[c + 1] - n0;
  for (int t0 = 0; t0 < cnt; t0 += blockDim.x) {
    int i = n0 + t0 + threadIdx.x, ok = t0 + threadIdx.x < cnt;
    double xi = ok ? t.x[i] : 0, yi = ok ? t.y[i] : 0, de = ok ? t.h[i] : 1, oi = ok ? t.o[i] : 0;
    int fi = ok ? t.f[i] : 0;
    double2 w = mk(0, 0);
    double dgi = 0, de2 = de * de, ce = t.nu / (de2 * de2);
    for (int pass = 0; pass < 4; pass++) {
      const int *po = pass == 0 ? t.uo : (pass == 1 ? t.xo : (pass == 2 ? t.ro : t.co)),
                *pl = pass == 0 ? t.ul : (pass == 1 ? t.xl : (pass == 2 ? t.rl : t.cl)),
                *ps = pass == 0 ? t.us : (pass == 1 ? t.xs : (pass == 2 ? t.rs : t.cs));
      for (int e = po[c]; e < po[c + 1]; e++) {
        int s = pl[e], sc = ps[e] % 9, pse = pass < 3 || ps[e] >= 9;
        double shx = (sc / 3 - 1) * t.L, shy = (sc % 3 - 1) * t.L;
        int m0 = t.start[s], mc = t.start[s + 1] - m0;
        for (int s0 = 0; s0 < mc; s0 += blockDim.x) {
          int js = m0 + s0 + threadIdx.x;
          __syncthreads();
          if (s0 + threadIdx.x < mc) {
            sx[threadIdx.x] = t.x[js] + shx;
            sy[threadIdx.x] = t.y[js] + shy;
            sg[threadIdx.x] = t.g[js];
            sv[threadIdx.x] = t.v[js];
            sf[threadIdx.x] = t.f[js];
          }
          __syncthreads();
          int m = mc - s0 < (int)blockDim.x ? mc - s0 : (int)blockDim.x;
          if (ok)
            for (int j = 0; j < m; j++) {
              int ff = fi && sf[j];
              if ((t.mode == 1 && ff) || (t.mode == 2 && !ff))
                continue;
              double dx = xi - sx[j], dy = yi - sy[j], r2 = dx * dx + dy * dy;
              if (r2 > 0) {
                double xx = r2 / de2, ex = exp(-xx);
                double k = nR > 0 ? (pass < 3 ? kcore(xx, ex) : -ex * poly(cR, nR, xx)) : (pass < 3);
                double2 q = mk(0, -sg[j] / (2 * pi));
                if (pass != 2)
                  w = add(w, scl(mul(q, mk(dx / r2, -dy / r2)), k));
                if (pse)
                  dgi += ce * (sg[j] - sv[j] * oi) * ex * poly(cE, nE, xx);
              }
            }
        }
      }
    }
    if (ok) {
      t.w[i] = add(t.mode == 2 ? mk(0, 0) : t.w[i], w);
      t.dg[i] = dgi;
    }
  }
}
__global__ static void direct(T t, long m, const double2 *lf, double *out) {
  __shared__ double rx[nw], ry[nw];
  long i = (long)blockIdx.x * (t.N / m);
  double xi = t.ox[i], yi = t.oy[i], hi = t.oh[i];
  double2 W = mk(0, 0);
  for (long j = threadIdx.x; j < t.N; j += blockDim.x)
    for (int sx = -1; sx <= 1; sx++)
      for (int sy = -1; sy <= 1; sy++) {
        double dx = xi - t.ox[j] - sx * t.L, dy = yi - t.oy[j] - sy * t.L,
               r2 = dx * dx + dy * dy;
        if (r2 > 0) {
          double xx = r2 / (hi * hi);
          double k = nR > 0 ? kcore(xx, exp(-xx)) : 1;
          W = add(W, scl(mul(mk(0, -t.og[j] / (2 * pi)), mk(dx / r2, -dy / r2)), k));
        }
      }
  rx[threadIdx.x] = W.x;
  ry[threadIdx.x] = W.y;
  __syncthreads();
  for (int h = blockDim.x / 2; h > 0; h /= 2) {
    if ((int)threadIdx.x < h) {
      rx[threadIdx.x] += rx[threadIdx.x + h];
      ry[threadIdx.x] += ry[threadIdx.x + h];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    int P = t.p + 1;
    double2 d = mk((xi - t.L / 2) / t.L, (yi - t.L / 2) / t.L), w = mk(0, 0);
    for (int l = t.p; l > 0; l--)
      w = add(mul(w, d), scl(lf[l], l));
    w = scl(w, 1 / t.L);
    out[4 * blockIdx.x] = rx[0] + w.x;
    out[4 * blockIdx.x + 1] = ry[0] + w.y;
    out[4 * blockIdx.x + 2] = t.ow[i].x;
    out[4 * blockIdx.x + 3] = t.ow[i].y;
    (void)P;
  }
}
__global__ static void fastpair(T t, int sub) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i >= t.N)
    return;
  double2 w = mk(0, 0);
  double dgi = 0;
  if (t.of[i]) {
    int b = i / NB;
    double xi = t.ox[i], yi = t.oy[i], de = t.oh[i], oi = t.oo[i], de2 = de * de;
    for (int e = t.uo[b]; e < t.uo[b + 1]; e++) {
      int c = t.ul[e];
      for (long j = (long)c * NB; j < (long)(c + 1) * NB; j++) {
        if (!t.of[j] || j == i)
          continue;
        double dx = xi - t.ox[j], dy = yi - t.oy[j];
        dx -= t.L * rint(dx / t.L);
        dy -= t.L * rint(dy / t.L);
        double r2 = dx * dx + dy * dy;
        if (r2 > 0) {
          double xx = r2 / de2, ex = exp(-xx);
          double k = nR > 0 ? kcore(xx, ex) : 1;
          w = add(w, scl(mul(mk(0, -t.og[j] / (2 * pi)), mk(dx / r2, -dy / r2)), k));
          dgi += t.nu / (de2 * de2) * (t.og[j] - t.ov[j] * oi) * ex * poly(cE, nE, xx);
        }
      }
    }
  }
  if (sub) {
    t.ow[i] = sub2(t.ow[i], w);
    t.odg[i] -= dgi;
  } else {
    t.ow[i] = w;
    t.odg[i] = dgi;
  }
}
__global__ static void iota(long N, long *a) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i < N)
    a[i] = i;
}
__global__ static void starts(long N, int nl, const int *k, int *st) {
  int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c > nl)
    return;
  long lo = 0, hi = N;
  while (lo < hi) {
    long m = (lo + hi) / 2;
    if (k[m] < c)
      lo = m + 1;
    else
      hi = m;
  }
  st[c] = (int)lo;
}
__global__ static void sums3(long N, const double *g, const double *x, const double *y, double *o) {
  __shared__ double a[1024], b[1024], c[1024];
  double sa = 0, sb = 0, sc = 0;
  for (long i = threadIdx.x; i < N; i += blockDim.x) {
    sa += g[i];
    sb += g[i] * x[i];
    sc += g[i] * y[i];
  }
  a[threadIdx.x] = sa;
  b[threadIdx.x] = sb;
  c[threadIdx.x] = sc;
  __syncthreads();
  for (int h = blockDim.x / 2; h > 0; h /= 2) {
    if ((int)threadIdx.x < h) {
      a[threadIdx.x] += a[threadIdx.x + h];
      b[threadIdx.x] += b[threadIdx.x + h];
      c[threadIdx.x] += c[threadIdx.x + h];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    o[0] = a[0];
    o[1] = b[0];
    o[2] = c[0];
  }
}
__global__ static void dmax1(long N, const double *x, double *o) {
  __shared__ double a[1024];
  double s = 0;
  for (long i = threadIdx.x; i < N; i += blockDim.x)
    s = fmax(s, x[i]);
  a[threadIdx.x] = s;
  __syncthreads();
  for (int h = blockDim.x / 2; h > 0; h /= 2) {
    if ((int)threadIdx.x < h)
      a[threadIdx.x] = fmax(a[threadIdx.x], a[threadIdx.x + h]);
    __syncthreads();
  }
  if (threadIdx.x == 0)
    o[0] = a[0];
}
__global__ static void fcount(T t, int *fc) {
  int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= t.nleaf)
    return;
  int s = 0;
  for (int i = t.start[c]; i < t.start[c + 1]; i++)
    s += t.f[i];
  fc[c] = s;
}
__global__ static void unsort(T t) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i >= t.N)
    return;
  long j = t.perm[i];
  t.ow[j] = t.w[i];
  t.odg[j] = t.dg[i];
}
__global__ static void zbar(T t, const double *sm) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i >= t.N)
    return;
  double gt = sm[0], Dx = sm[1], Dy = sm[2];
  double A = t.L * t.L;
  double2 a = mk(gt * t.ox[i] - Dx, -gt * t.oy[i] + Dy);
  t.ow[i] = add(t.ow[i], mul(mk(0, 1 / (2 * A)), a));
}
__global__ static void spread(int n, const double2 *w, double2 *a, double2 *b) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int nh = n / 2 + 1;
  if (l >= (long)n * nh)
    return;
  long i = l % nh, j = l / nh;
  double kx = i, ky = j <= n / 2 ? j : j - n, kk = kx * kx + ky * ky;
  double2 W = w[l];
  double2 P = kk > 0 ? make_double2(W.x / kk, W.y / kk) : make_double2(0, 0);
  a[l] = make_double2(-ky * P.y, ky * P.x);
  b[l] = make_double2(kx * P.y, -kx * P.x);
}
static double binom(int n, int k) {
  double c = 1;
  for (int i = 1; i <= k; i++)
    c = c * (n - k + i) / i;
  return c;
}
static double2 cpow(double2 a, int n) {
  double2 r = mk(1, 0);
  for (int i = 0; i < n; i++)
    r = mul(r, a);
  return r;
}
static void solve(int n, double *A, double *b) {
  for (int c = 0; c < n; c++) {
    int pv = c;
    for (int r = c + 1; r < n; r++)
      if (fabs(A[r * n + c]) > fabs(A[pv * n + c]))
        pv = r;
    for (int k = 0; k < n; k++) {
      double s = A[c * n + k];
      A[c * n + k] = A[pv * n + k];
      A[pv * n + k] = s;
    }
    double s = b[c];
    b[c] = b[pv];
    b[pv] = s;
    for (int r = c + 1; r < n; r++) {
      double f = A[r * n + c] / A[c * n + c];
      for (int k = c; k < n; k++)
        A[r * n + k] -= f * A[c * n + k];
      b[r] -= f * b[c];
    }
  }
  for (int c = n - 1; c >= 0; c--) {
    for (int k = c + 1; k < n; k++)
      b[c] -= A[c * n + k] * b[k];
    b[c] /= A[c * n + c];
  }
}
static double fact(int n) { return n <= 1 ? 1 : n * fact(n - 1); }
static void cores(int ord, double *R, int *nr, double *E, int *ne) {
  int m = ord / 2;
  double A[64], b[8], c[8];
  for (int k = 0; k < m; k++) {
    for (int j = 0; j < m; j++)
      A[k * m + j] = fact(j + k);
    b[k] = k == 0 ? 1 : 0;
  }
  solve(m, A, b);
  for (int j = 0; j < m; j++)
    c[j] = b[j];
  for (int i = 0; i < m; i++) {
    R[i] = 0;
    for (int j = i; j < m; j++)
      R[i] += c[j] * fact(j) / fact(i);
  }
  *nr = m;
  for (int k = 1; k <= m; k++) {
    for (int j = 0; j < m; j++)
      A[(k - 1) * m + j] = pi * fact(j + k);
    b[k - 1] = k == 1 ? 4 : 0;
  }
  solve(m, A, b);
  for (int j = 0; j < m; j++)
    E[j] = b[j];
  *ne = m;
}
static int pair_;
static double omega(double x, double y, double rho, double eps) {
  if (pair_) {
    double s = -2.0 / (4 * pi * pi);
    for (int k = 0; k < 2; k++)
      for (int a = -1; a <= 1; a++)
        for (int b = -1; b <= 1; b++) {
          double dx = x - pi - (k ? eps : -eps) + 2 * pi * a, dy = y - pi + 2 * pi * b;
          s += exp(-(dx * dx + dy * dy) / (rho * rho)) / (pi * rho * rho);
        }
    return s;
  }
  double s = eps * cos(x);
  for (int n = -3; n <= 3; n++) {
    double a = cosh((y - pi / 2 + 2 * pi * n) / rho),
           b = cosh((y - 3 * pi / 2 + 2 * pi * n) / rho);
    s -= (1 / (a * a) - 1 / (b * b)) / rho;
  }
  return s;
}
static long steps(const char *opt, double x, double dt) {
  long n = lround(x / dt);
  if (fabs(n * dt - x) > 1e-9 * x) {
    fprintf(stderr, "vm2: error: %s %g is not a multiple of -s %g\n", opt, x,
            dt);
    exit(1);
  }
  return n;
}
static double now(void) {
  struct timeval v;
  gettimeofday(&v, NULL);
  return v.tv_sec + 1e-6 * v.tv_usec;
}
__global__ static void remesh(long N, int W, int G, double h0, const int *lev, const double *lw,
                              const double *lx, const double *ly, double *o) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  int b = n / NB, i = n % NB % B, j = n % NB / B;
  double h = h0 / (1 << lev[b]), s = 0;
  long c0 = (long)b * W * W;
  int i0 = i - G + 1 > -G ? i - G + 1 : -G, i1 = i + G - 1 < B + G - 1 ? i + G - 1 : B + G - 1;
  int j0 = j - G + 1 > -G ? j - G + 1 : -G, j1 = j + G - 1 < B + G - 1 ? j + G - 1 : B + G - 1;
  for (int cj = j0; cj <= j1; cj++)
    for (int ci = i0; ci <= i1; ci++) {
      long k = c0 + (long)(cj + G) * W + ci + G;
      double dx = i - ci - lx[k] / h, dy = j - cj - ly[k] / h;
      if (fabs(dx) < 4 && fabs(dy) < 4) {
        double fx = dx - floor(dx), fy = dy - floor(dy), sx = 0, sy = 0;
        for (int q = -3; q <= 4; q++) {
          sx += kw(fx - q);
          sy += kw(fy - q);
        }
        s += lw[k] / (sx * sy) * kw(dx) * kw(dy);
      }
    }
  o[n] = s;
}
__global__ static void dsp(long N, double L, const double *x, const double *y, const double *x0,
                           const double *y0, double *dx, double *dy) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  double a = x[n] - x0[n], c = y[n] - y0[n];
  dx[n] = a - L * rint(a / L);
  dy[n] = c - L * rint(c / L);
}
__global__ static void gfix(long N, double c0, const double *sm, double *o) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n < N)
    o[n] += (c0 - sm[1]) / sm[0];
}
__global__ static void gam(long N, const double *o, const double *v, double *g) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n < N)
    g[n] = o[n] * v[n];
}
__global__ static void labown(long N, int W, int G, const double *v, double *lab) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  int b = n / NB, i = n % NB % B, j = n % NB / B;
  lab[(long)b * W * W + (j + G) * W + i + G] = v[n];
}
__global__ static void labfill(long nr, const long *gr, const int *gp, const int *gc,
                               const double *gw, const double *v, double *lab) {
  long r = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (r >= nr)
    return;
  double s = 0;
  for (int e = gp[r]; e < gp[r + 1]; e++)
    s += gw[e] * v[gc[e]];
  lab[gr[r]] = s;
}
__global__ static void grad(long N, int W, int G, double h0, const int *lev,
                            const double *lu, const double *lv, double *out) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  int b = n / NB, i = n % NB % B, j = n % NB / B;
  double h = h0 / (1 << lev[b]);
  long c = (long)b * W * W + (j + G) * W + i + G;
  double ux = (lu[c + 1] - lu[c - 1]) / (2 * h), uy = (lu[c + W] - lu[c - W]) / (2 * h);
  double vx = (lv[c + 1] - lv[c - 1]) / (2 * h), vy = (lv[c + W] - lv[c - W]) / (2 * h);
  out[n] = sqrt(ux * ux + uy * uy + vx * vx + vy * vy);
}
__global__ static void detail(long N, int W, int G, const double *lw, double *out) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  int b = n / NB, i = n % NB % B, j = n % NB / B, oi = i & 1, oj = j & 1;
  const double DW[8] = {-5. / 2048, 49. / 2048, -245. / 2048, 1225. / 2048,
                        1225. / 2048, -245. / 2048, 49. / 2048, -5. / 2048};
  long c0 = (long)b * W * W;
  if (!oi && !oj) {
    out[n] = 0;
    return;
  }
  double p = 0;
  for (int a = 0; a < (oi ? 8 : 1); a++)
    for (int e = 0; e < (oj ? 8 : 1); e++) {
      int ii = oi ? i - 7 + 2 * a : i, jj = oj ? j - 7 + 2 * e : j;
      p += (oi ? DW[a] : 1) * (oj ? DW[e] : 1) * lw[c0 + (jj + G) * W + ii + G];
    }
  out[n] = fabs(lw[c0 + (j + G) * W + i + G] - p);
}
__global__ static void strain(long N, int W, int G, double h0, const int *lev, const double *lu,
                              const double *lv, double *out) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  int b = n / NB, i = n % NB % B, j = n % NB / B;
  double h = h0 / (1 << lev[b]);
  long c = (long)b * W * W + (j + G) * W + i + G;
  double ux = (lu[c + 1] - lu[c - 1]) / (2 * h), uy = (lu[c + W] - lu[c - W]) / (2 * h);
  double vx = (lv[c + 1] - lv[c - 1]) / (2 * h), vy = (lv[c + W] - lv[c - W]) / (2 * h);
  out[n] = sqrt(0.25 * (ux - vy) * (ux - vy) + 0.25 * (uy + vx) * (uy + vx));
}
__global__ static void nodes(long N, double h0, const int *lev, const int *bx,
                             const int *by, double *x, double *y, double *h) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  int b = n / NB, i = n % NB % B, j = n % NB / B;
  double hb = h0 / (1 << lev[b]);
  x[n] = (bx[b] * B + i) * hb;
  y[n] = (by[b] * B + j) * hb;
  h[n] = hb;
}
__global__ static void comb(long N, double *o, const double *a, double c1,
                            const double *k1, double c2, const double *k2) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  o[n] = a[n] + c1 * k1[n] + (k2 ? c2 * k2[n] : 0);
}
__global__ static void split(long N, const double2 *w, const double *dg,
                             double *ku, double *kv, double *kg) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  ku[n] = w[n].x;
  kv[n] = -w[n].y;
  kg[n] = dg[n];
}
__global__ static void disp(long N, double L, const double *x, const double *y,
                            const double *x0, const double *y0, double h0, const int *lev,
                            double *o) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  double dx = x[n] - x0[n], dy = y[n] - y0[n];
  dx -= L * rint(dx / L);
  dy -= L * rint(dy / L);
  o[n] = fmax(fabs(dx), fabs(dy)) / (h0 / (1 << lev[n / NB]));
}
struct H {
  int nb0, Lm, D0, Dm, nleaf, nbox, *lev, *bx, *by, **map, **bmap, *bd, *px,
      *py, *ch, *par, *leaf, *lbox, lo[DM + 2], hi[DM + 2];
  int *uo, *ul, *us, *vo, *vl, *vs, *wo, *wl, *ws, *xo, *xl, *xs, *ro, *rl, *rs, *co, *cl, *cs;
  int *rel, *nbid;
};
static int findb(H *m, int d, int x, int y, int *sc) {
  int n = 1 << d, xx = ((x % n) + n) % n, yy = ((y % n) + n) % n;
  *sc = ((x - xx) / n + 1) * 3 + (y - yy) / n + 1;
  return m->bmap[d][(long)yy * n + xx];
}
static int adjb(H *m, int a, int b, int sc) {
  long F = 1L << m->Dm, sa = m->Dm - m->bd[a], sb = m->Dm - m->bd[b];
  long ax0 = (long)m->px[a] << sa, ax1 = ((long)m->px[a] + 1) << sa;
  long ay0 = (long)m->py[a] << sa, ay1 = ((long)m->py[a] + 1) << sa;
  long bx0 = ((long)m->px[b] << sb) + (sc / 3 - 1) * F,
       bx1 = (((long)m->px[b] + 1) << sb) + (sc / 3 - 1) * F;
  long by0 = ((long)m->py[b] << sb) + (sc % 3 - 1) * F,
       by1 = (((long)m->py[b] + 1) << sb) + (sc % 3 - 1) * F;
  return ax0 <= bx1 && bx0 <= ax1 && ay0 <= by1 && by0 <= ay1;
}
static void push(int **l, int **s, int *n, int *cap, int v, int sc) {
  if (*n >= *cap) {
    *cap = *cap ? 2 * *cap : 1024;
    *l = (int *)realloc(*l, *cap * sizeof(int));
    *s = (int *)realloc(*s, *cap * sizeof(int));
    if (*l == NULL || *s == NULL) {
      fprintf(stderr, "vm2: error: malloc failed\n");
      exit(1);
    }
  }
  (*l)[*n] = v;
  (*s)[*n] = sc;
  (*n)++;
}
static int *lmap(H *m, int l, long x, long y) {
  long n = (long)m->nb0 << l;
  if (l < 0 || l > m->Lm)
    return NULL;
  x = ((x % n) + n) % n;
  y = ((y % n) + n) % n;
  return &m->map[l][y * n + x];
}
static int nbuf = 4;
static void build(H *m, int nb0, int Lm, double w, const int *tin) {
  int D0 = 0;
  while ((1 << D0) < nb0)
    D0++;
  if ((1 << D0) != nb0) {
    fprintf(stderr, "vm2: error: -n / %d must be a power of two\n", B);
    exit(1);
  }
  m->nb0 = nb0;
  m->Lm = Lm;
  m->D0 = D0;
  m->Dm = D0 + Lm;
  if (m->Dm > DM) {
    fprintf(stderr, "vm2: error: tree too deep\n");
    exit(1);
  }
  long nf = (long)nb0 << Lm;
  int *tgt = (int *)hmalloc(nf * nf * sizeof(int)), *lv = (int *)hmalloc(nf * nf * sizeof(int));
  for (long j = 0; j < nf; j++)
    for (long i = 0; i < nf; i++) {
      double y0 = 2 * pi * j / nf, y1 = 2 * pi * (j + 1) / nf, dm = 1e300;
      for (int q = 0; q < 2; q++) {
        double c = q ? 3 * pi / 2 : pi / 2;
        dm = fmin(dm, c < y0 ? y0 - c : (c > y1 ? c - y1 : 0));
      }
      tgt[j * nf + i] = tin ? tin[j * nf + i] : (dm < w ? Lm : 0);
    }
  for (;;) {
    for (long j = 0; j < nf; j++)
      for (long i = 0; i < nf; i++) {
        int l = 0;
        for (; l < Lm; l++) {
          long s = 1L << (Lm - l), i0 = i / s * s, j0 = j / s * s, mx = 0;
          for (long jj = j0; jj < j0 + s; jj++)
            for (long ii = i0; ii < i0 + s; ii++)
              mx = mx > tgt[jj * nf + ii] ? mx : tgt[jj * nf + ii];
          if (mx <= l)
            break;
        }
        lv[j * nf + i] = l;
      }
    int ok = 1;
    for (long j = 0; j < nf; j++)
      for (long i = 0; i < nf; i++) {
        int l = lv[j * nf + i];
        long s = 1L << (Lm - l), i0 = i / s * s, j0 = j / s * s, mx = 0;
        if (i != i0 || j != j0)
          continue;
        long r = l > 1 ? s * nbuf : s;
        for (long jj = j0 - r; jj < j0 + s + r; jj++)
          for (long ii = i0 - r; ii < i0 + s + r; ii++) {
            long q = ((jj + nf) % nf) * nf + (ii + nf) % nf;
            if (jj >= j0 - s && jj < j0 + 2 * s && ii >= i0 - s && ii < i0 + 2 * s)
              mx = mx > lv[q] ? mx : lv[q];
            if (lv[q] < l - 1) {
              tgt[q] = l - 1;
              ok = 0;
            }
          }
        if (mx > l + 1) {
          for (long jj = j0; jj < j0 + s; jj++)
            for (long ii = i0; ii < i0 + s; ii++)
              tgt[jj * nf + ii] = mx - 1;
          ok = 0;
        }
      }
    if (ok)
      break;
  }
  int cap = 0;
  m->nleaf = 0;
  m->lev = m->bx = m->by = NULL;
  m->map = (int **)hmalloc((Lm + 1) * sizeof(int *));
  for (int l = 0; l <= Lm; l++) {
    long n = (long)nb0 << l;
    m->map[l] = (int *)hmalloc(n * n * sizeof(int));
    for (long k = 0; k < n * n; k++)
      m->map[l][k] = -1;
  }
  for (int l = 0; l <= Lm; l++) {
    long n = (long)nb0 << l, s = 1L << (Lm - l);
    for (long j = 0; j < n; j++)
      for (long i = 0; i < n; i++)
        if (lv[j * s * nf + i * s] == l) {
          if (m->nleaf >= cap) {
            cap = cap ? 2 * cap : 1024;
            m->lev = (int *)realloc(m->lev, cap * sizeof(int));
            m->bx = (int *)realloc(m->bx, cap * sizeof(int));
            m->by = (int *)realloc(m->by, cap * sizeof(int));
          }
          m->lev[m->nleaf] = l;
          m->bx[m->nleaf] = i;
          m->by[m->nleaf] = j;
          m->map[l][j * n + i] = m->nleaf++;
        }
  }
  free(tgt);
  free(lv);
  int Dm = m->Dm, nbx = 1, cb = 1024;
  m->bd = (int *)hmalloc(cb * sizeof(int));
  m->px = (int *)hmalloc(cb * sizeof(int));
  m->py = (int *)hmalloc(cb * sizeof(int));
  m->ch = (int *)hmalloc(cb * sizeof(int));
  m->par = (int *)hmalloc(cb * sizeof(int));
  m->leaf = (int *)hmalloc(cb * sizeof(int));
  m->bmap = (int **)hmalloc((Dm + 1) * sizeof(int *));
  for (int d = 0; d <= Dm; d++) {
    long n = 1L << d;
    m->bmap[d] = (int *)hmalloc(n * n * sizeof(int));
    for (long k = 0; k < n * n; k++)
      m->bmap[d][k] = -1;
  }
  m->bd[0] = m->px[0] = m->py[0] = 0;
  m->par[0] = -1;
  m->bmap[0][0] = 0;
  m->lo[0] = 0;
  m->hi[0] = 1;
  for (int d = 0; d <= Dm; d++) {
    for (int b = m->lo[d]; b < m->hi[d]; b++) {
      int lf = -1;
      if (d >= D0)
        lf = m->map[d - D0][(long)m->py[b] * ((long)nb0 << (d - D0)) + m->px[b]];
      m->leaf[b] = lf;
      if (lf >= 0) {
        m->ch[b] = -1;
        continue;
      }
      if (d == Dm) {
        fprintf(stderr, "vm2: error: tree build\n");
        exit(1);
      }
      if (nbx + 4 > cb) {
        cb *= 2;
        m->bd = (int *)realloc(m->bd, cb * sizeof(int));
        m->px = (int *)realloc(m->px, cb * sizeof(int));
        m->py = (int *)realloc(m->py, cb * sizeof(int));
        m->ch = (int *)realloc(m->ch, cb * sizeof(int));
        m->par = (int *)realloc(m->par, cb * sizeof(int));
        m->leaf = (int *)realloc(m->leaf, cb * sizeof(int));
      }
      m->ch[b] = nbx;
      for (int q = 0; q < 4; q++) {
        m->bd[nbx] = d + 1;
        m->px[nbx] = 2 * m->px[b] + (q >> 1);
        m->py[nbx] = 2 * m->py[b] + (q & 1);
        m->par[nbx] = b;
        m->bmap[d + 1][(long)m->py[nbx] * (1L << (d + 1)) + m->px[nbx]] = nbx;
        nbx++;
      }
    }
    m->lo[d + 1] = m->hi[d];
    m->hi[d + 1] = nbx;
  }
  m->nbox = nbx;
  m->lbox = (int *)hmalloc(m->nleaf * sizeof(int));
  for (int b = 0; b < nbx; b++)
    if (m->leaf[b] >= 0)
      m->lbox[m->leaf[b]] = b;
  int nu = 0, cu = 0, nv = 0, cv = 0, nwn = 0, cw = 0, nx = 0, cx = 0, nx0 = 0, cx0 = 0, nxb = 0, cxb = 0, dum;
  int *ul = NULL, *us = NULL, *vl = NULL, *vs = NULL, *wl = NULL, *ws = NULL, *xl = NULL, *xs = NULL;
  int *xl0 = NULL, *xs0 = NULL, *xb0 = NULL, *xbs = NULL;
  m->uo = (int *)hmalloc((m->nleaf + 1) * sizeof(int));
  m->wo = (int *)hmalloc((m->nleaf + 1) * sizeof(int));
  m->xo = (int *)hmalloc((m->nleaf + 1) * sizeof(int));
  m->vo = (int *)hmalloc((nbx + 1) * sizeof(int));
  static int stk[1 << 16], sts[1 << 16];
  for (int c = 0; c < m->nleaf; c++) {
    int b = m->lbox[c], d = m->bd[b];
    m->uo[c] = nu;
    m->wo[c] = nwn;
    for (int dd = 0; dd < d; dd++) {
      int ax = m->px[b] >> (d - dd), ay = m->py[b] >> (d - dd);
      for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++) {
          int sc, e = findb(m, dd, ax + dx, ay + dy, &sc), dup = 0;
          if (e < 0 || m->ch[e] >= 0 || !adjb(m, b, e, sc))
            continue;
          for (int k = m->uo[c]; k < nu; k++)
            dup |= ul[k] == m->leaf[e] && us[k] == sc;
          if (!dup)
            push(&ul, &us, &nu, &cu, m->leaf[e], sc);
        }
    }
    int ns = 0;
    for (int dx = -1; dx <= 1; dx++)
      for (int dy = -1; dy <= 1; dy++) {
        int sc, e = findb(m, d, m->px[b] + dx, m->py[b] + dy, &sc);
        if (e >= 0) {
          stk[ns] = e;
          sts[ns++] = sc;
        }
      }
    while (ns > 0) {
      int e = stk[--ns], sc = sts[ns];
      if (m->ch[e] < 0) {
        push(&ul, &us, &nu, &cu, m->leaf[e], sc);
        continue;
      }
      for (int q = 0; q < 4; q++) {
        int k = m->ch[e] + q;
        if (adjb(m, b, k, sc)) {
          if (ns >= (1 << 16)) {
            fprintf(stderr, "vm2: error: list stack\n");
            exit(1);
          }
          stk[ns] = k;
          sts[ns++] = sc;
        } else
          push(&wl, &ws, &nwn, &cw, k, sc);
      }
    }
  }
  m->uo[m->nleaf] = nu;
  m->wo[m->nleaf] = nwn;
  for (int b = 0; b < nbx; b++) {
    int d = m->bd[b], pb = m->par[b];
    m->vo[b] = nv;
    if (d == 0)
      continue;
    for (int dx = -1; dx <= 1; dx++)
      for (int dy = -1; dy <= 1; dy++) {
        int sc, e = findb(m, d - 1, m->px[pb] + dx, m->py[pb] + dy, &sc);
        if (e < 0 || m->ch[e] < 0)
          continue;
        for (int q = 0; q < 4; q++)
          if (!adjb(m, b, m->ch[e] + q, sc))
            push(&vl, &vs, &nv, &cv, m->ch[e] + q, sc);
      }
    for (int dd = 0; dd < d; dd++) {
      int ax = m->px[b] >> (d - dd), ay = m->py[b] >> (d - dd);
      for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++) {
          int sc, e = findb(m, dd, ax + dx, ay + dy, &sc), dup = 0;
          if (e < 0 || m->ch[e] >= 0)
            continue;
          if (!adjb(m, pb, e, sc) || adjb(m, b, e, sc))
            continue;
          for (int k = 0; k < nx0; k++)
            dup |= xb0[k] == b && xl0[k] == m->leaf[e] && xs0[k] == sc;
          if (!dup) {
            push(&xl0, &xs0, &nx0, &cx0, m->leaf[e], sc);
            push(&xb0, &xbs, &nxb, &cxb, b, 0);
          }
        }
    }
  }
  m->vo[nbx] = nv;
  for (int c = 0; c < m->nleaf; c++) {
    m->xo[c] = nx;
    int b = m->lbox[c], bd = m->bd[b];
    for (int k = 0; k < nx0; k++) {
      int a = xb0[k], d = m->bd[a];
      if (bd < d || (m->px[b] >> (bd - d)) != m->px[a] || (m->py[b] >> (bd - d)) != m->py[a])
        continue;
      push(&xl, &xs, &nx, &cx, xl0[k], xs0[k]);
    }
  }
  m->xo[m->nleaf] = nx;
  m->ro = (int *)hmalloc((m->nleaf + 1) * sizeof(int));
  m->rl = (int *)hmalloc((nx + 1) * sizeof(int));
  m->rs = (int *)hmalloc((nx + 1) * sizeof(int));
  for (int c = 0; c <= m->nleaf; c++)
    m->ro[c] = 0;
  for (int c = 0; c < m->nleaf; c++)
    for (int e = m->xo[c]; e < m->xo[c + 1]; e++)
      m->ro[xl[e] + 1]++;
  for (int c = 0; c < m->nleaf; c++)
    m->ro[c + 1] += m->ro[c];
  {
    int *fill = (int *)hmalloc((m->nleaf + 1) * sizeof(int));
    memcpy(fill, m->ro, (m->nleaf + 1) * sizeof(int));
    for (int c = 0; c < m->nleaf; c++)
      for (int e = m->xo[c]; e < m->xo[c + 1]; e++) {
        int k = fill[xl[e]]++;
        m->rl[k] = c;
        m->rs[k] = (2 - xs[e] / 3) * 3 + 2 - xs[e] % 3;
      }
    free(fill);
  }
  (void)dum;
  m->ul = ul;
  m->us = us;
  m->vl = vl;
  m->vs = vs;
  m->wl = wl;
  m->ws = ws;
  m->xl = xl;
  m->xs = xs;
  free(xl0);
  free(xs0);
  free(xb0);
  free(xbs);
  m->rel = (int *)hmalloc(m->nleaf * 9 * sizeof(int));
  m->nbid = (int *)hmalloc(m->nleaf * 9 * 4 * sizeof(int));
  for (int c = 0; c < m->nleaf; c++) {
    int l = m->lev[c];
    for (int k = 0; k < 9; k++) {
      int cxx = k / 3 - 1, cyy = k % 3 - 1, *p;
      long X = m->bx[c] + cxx, Y = m->by[c] + cyy;
      m->rel[c * 9 + k] = 0;
      for (int q = 0; q < 4; q++)
        m->nbid[(c * 9 + k) * 4 + q] = -1;
      if (k == 4)
        continue;
      if ((p = lmap(m, l, X, Y)) != NULL && *p >= 0) {
        m->rel[c * 9 + k] = 1;
        m->nbid[(c * 9 + k) * 4] = *p;
      } else if ((p = lmap(m, l - 1, (X % ((long)nb0 << l) + ((long)nb0 << l)) % ((long)nb0 << l) >> 1,
                           (Y % ((long)nb0 << l) + ((long)nb0 << l)) % ((long)nb0 << l) >> 1)) != NULL &&
                 *p >= 0) {
        m->rel[c * 9 + k] = 2;
        m->nbid[(c * 9 + k) * 4] = *p;
      } else {
        m->rel[c * 9 + k] = 3;
        for (int q = 0; q < 4; q++) {
          p = lmap(m, l + 1, 2 * X + (q & 1), 2 * Y + (q >> 1));
          if (p == NULL || *p < 0) {
            fprintf(stderr, "vm2: error: mesh is not 2:1 balanced\n");
            exit(1);
          }
          m->nbid[(c * 9 + k) * 4 + q] = *p;
        }
      }
    }
  }
}
struct Tab {
  int hd[8], same[9][2], fine[9][2], cself[2], csame[9][2], cfine[9][2],
      ccoarse[9][4][2], pred[9][2];
  int *ops;
  double *q, R[4][8], E[4][8], *S;
  int ns, Pt;
  double2 *F;
};
static void tabload(Tab *t, const char *path) {
  FILE *f;
  if ((f = fopen(path, "rb")) == NULL) {
    fprintf(stderr, "vm2: error: fail to open '%s' (run gen_vm2.py)\n", path);
    exit(1);
  }
  size_t n = offsetof(Tab, ops);
  if (fread(t, 1, n, f) != n || t->hd[0] != 0x564D3254 || t->hd[1] != B || t->hd[3] != 8) {
    fprintf(stderr, "vm2: error: bad table '%s'\n", path);
    exit(1);
  }
  t->ops = (int *)hmalloc((size_t)t->hd[7] * 14 * sizeof(int));
  t->q = (double *)hmalloc(t->hd[4] * sizeof(double));
  if (fread(t->ops, sizeof(int), (size_t)t->hd[7] * 14, f) != (size_t)t->hd[7] * 14 ||
      fread(t->q, sizeof(double), t->hd[4], f) != (size_t)t->hd[4] ||
      fread(t->R, sizeof(double), 32, f) != 32 || fread(t->E, sizeof(double), 32, f) != 32 ||
      fread(&t->ns, sizeof(int), 1, f) != 1 || t->ns < 2 * PM + 4 ||
      (t->S = (double *)hmalloc(t->ns * sizeof(double))) == NULL ||
      fread(t->S, sizeof(double), t->ns, f) != (size_t)t->ns || fread(&t->Pt, sizeof(int), 1, f) != 1 ||
      t->Pt < 2 || t->Pt > PM + 1 ||
      (t->F = (double2 *)hmalloc(58L * t->Pt * t->Pt * sizeof(double2))) == NULL ||
      fread(t->F, sizeof(double2), 58L * t->Pt * t->Pt, f) != (size_t)58 * t->Pt * t->Pt) {
    fprintf(stderr, "vm2: error: fail to read '%s'\n", path);
    exit(1);
  }
  fclose(f);
}
enum { SVM = 1024 };
struct SV {
  int n, i[SVM];
  double w[SVM];
};
static const double DDW[8] = {-5. / 2048, 49. / 2048, -245. / 2048, 1225. / 2048,
                              1225. / 2048, -245. / 2048, 49. / 2048, -5. / 2048};
static void sexec(Tab *T, const int *sec, SV **bufs, int *nb, int own) {
  for (int k = sec[0]; k < sec[0] + sec[1]; k++) {
    int *o = T->ops + 14L * k, ty = o[0], bd = o[1], dst = o[2], bs = o[3], src = o[4], cnt = o[5];
    if (bs >= 3) {
      int blk = bs == 3 ? own : nb[bs - 4];
      for (int e = 0; e < (ty == 0 ? cnt : 1); e++) {
        SV *d = &bufs[bd][dst + e];
        d->n = 1;
        d->i[0] = blk * NB + src + e;
        d->w[0] = 1;
      }
      continue;
    }
    SV *d = &bufs[bd][dst];
    if (ty == 4 || ty == 0) {
      for (int e = 0; e < cnt; e++)
        bufs[bd][dst + e] = bufs[bs][src + e];
      continue;
    }
    int st = ty == 2 ? 1 : o[6];
    d->n = 0;
    for (int q = 0; q < 8; q++) {
      SV *s = &bufs[bs][src + q * st];
      if (s->n < 0) {
        d->n = -1;
        break;
      }
      for (int e = 0; e < s->n; e++) {
        int f = 0;
        for (; f < d->n; f++)
          if (d->i[f] == s->i[e])
            break;
        if (f == d->n) {
          if (d->n >= SVM) {
            fprintf(stderr, "vm2: error: sparse row too long\n");
            exit(1);
          }
          d->i[f] = s->i[e];
          d->w[f] = 0;
          d->n++;
        }
        d->w[f] += DDW[q] * s->w[e];
      }
    }
  }
}
static const double RW[9] = {-5. / 4096, 49. / 4096, -245. / 4096, 1225. / 4096, 0.5,
                             1225. / 4096, -245. / 4096, 49. / 4096, -5. / 4096};
static const int RO[9] = {-7, -5, -3, -1, 0, 1, 3, 5, 7};
static int fdiv(int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }
static int cbnode(H *m, int b, int I, int J) {
  int cx = fdiv(I, B), cy = fdiv(J, B), k = (cx + 1) * 3 + cy + 1, r;
  if (cx < -1 || cx > 1 || cy < -1 || cy > 1) {
    fprintf(stderr, "vm2: error: restriction stencil leaves the neighbourhood\n");
    exit(1);
  }
  int li = I - cx * B, lj = J - cy * B;
  if (k == 4)
    return b * NB + lj * B + li;
  r = m->rel[b * 9 + k];
  if (r == 1)
    return m->nbid[(b * 9 + k) * 4] * NB + lj * B + li;
  if (r == 3) {
    int X = 2 * li, Y = 2 * lj, q = X / B + 2 * (Y / B);
    return m->nbid[(b * 9 + k) * 4 + q] * NB + (Y % B) * B + X % B;
  }
  fprintf(stderr, "vm2: error: restriction stencil meets a coarser block\n");
  exit(1);
}
static int fnode(H *m, int b, int X, int Y) {
  int cx = fdiv(X, 2 * B), cy = fdiv(Y, 2 * B), k = (cx + 1) * 3 + cy + 1;
  if (k != 4 && cx >= -1 && cx <= 1 && cy >= -1 && cy <= 1 && m->rel[b * 9 + k] == 3) {
    int Xl = X - cx * 2 * B, Yl = Y - cy * 2 * B, q = Xl / B + 2 * (Yl / B);
    return m->nbid[(b * 9 + k) * 4 + q] * NB + (Yl % B) * B + Xl % B;
  }
  return -1;
}
enum { RM = 4096 };
static int rn_, ri_[RM];
static double rw_[RM];
static void radd(int idx, double w) {
  for (int e = 0; e < rn_; e++)
    if (ri_[e] == idx) {
      rw_[e] += w;
      return;
    }
  if (rn_ >= RM) {
    fprintf(stderr, "vm2: error: restriction row too long\n");
    exit(1);
  }
  ri_[rn_] = idx;
  rw_[rn_++] = w;
}
static void ffield(H *m, int b, int X, int Y, double w) {
  int f = fnode(m, b, X, Y);
  if (f >= 0) {
    radd(f, w);
    return;
  }
  int ox = ((X % 2) + 2) % 2, oy = ((Y % 2) + 2) % 2;
  for (int a = 0; a < (ox ? 8 : 1); a++)
    for (int e = 0; e < (oy ? 8 : 1); e++) {
      int I = ox ? fdiv(X - 1, 2) - 3 + a : fdiv(X, 2), J = oy ? fdiv(Y - 1, 2) - 3 + e : fdiv(Y, 2);
      radd(cbnode(m, b, I, J), w * (ox ? DDW[a] : 1) * (oy ? DDW[e] : 1));
    }
}
struct CSR {
  long nr, nnz;
  long *row;
  int *p, *c;
  double *w;
};
static void smoothr(H *m, int b, int k, int G, int W, SV *lab) {
  int cx = k / 3 - 1, cy = k % 3 - 1;
  int i0 = cx < 0 ? -G : (cx == 0 ? 0 : B), i1 = cx < 0 ? 0 : (cx == 0 ? B : B + G);
  int j0 = cy < 0 ? -G : (cy == 0 ? 0 : B), j1 = cy < 0 ? 0 : (cy == 0 ? B : B + G);
  for (int j = j0; j < j1; j++)
    for (int i = i0; i < i1; i++) {
      rn_ = 0;
      for (int a = 0; a < 9; a++)
        for (int e = 0; e < 9; e++)
          ffield(m, b, 2 * i + RO[a], 2 * j + RO[e], RW[a] * RW[e]);
      SV *s = &lab[(j + G) * W + i + G];
      if (rn_ > SVM) {
        fprintf(stderr, "vm2: error: restriction row %d > %d\n", rn_, SVM);
        exit(1);
      }
      s->n = rn_;
      for (int e = 0; e < rn_; e++) {
        s->i[e] = ri_[e];
        s->w[e] = rw_[e];
      }
    }
}
static void compile(H *m, Tab *T, CSR *g, CSR *fo, int smooth) {
  int W = T->hd[4], G = T->hd[2], C = T->hd[6], WW = W * W;
  SV *lab = (SV *)hmalloc(WW * sizeof(SV)), *buf = (SV *)hmalloc(C * C * sizeof(SV)),
     *tmp = (SV *)hmalloc(W * C * sizeof(SV));
  SV *bufs[3] = {lab, buf, tmp};
  long cr = 1 << 20, cn = 1 << 22;
  g->row = (long *)hmalloc(cr * sizeof(long));
  g->p = (int *)hmalloc((cr + 1) * sizeof(int));
  g->c = (int *)hmalloc(cn * sizeof(int));
  g->w = (double *)hmalloc(cn * sizeof(double));
  g->nr = g->nnz = 0;
  for (int b = 0; b < m->nleaf; b++) {
    for (int k = 0; k < WW; k++)
      lab[k].n = -1;
    for (int k = 0; k < C * C; k++)
      buf[k].n = -1;
    for (int k = 0; k < W * C; k++)
      tmp[k].n = -1;
    int anyc = 0, par = (m->bx[b] & 1) + 2 * (m->by[b] & 1);
    for (int k = 0; k < 9; k++)
      anyc |= m->rel[b * 9 + k] == 2;
    if (anyc)
      sexec(T, T->cself, bufs, NULL, b);
    for (int k = 0; k < 9; k++) {
      int r = m->rel[b * 9 + k], *nb = m->nbid + (b * 9 + k) * 4;
      if (r == 1) {
        sexec(T, T->same[k], bufs, nb, b);
        if (anyc)
          sexec(T, T->csame[k], bufs, nb, b);
      } else if (r == 3) {
        sexec(T, T->fine[k], bufs, nb, b);
        if (anyc)
          sexec(T, T->cfine[k], bufs, nb, b);
      } else if (r == 2)
        sexec(T, T->ccoarse[k][par], bufs, nb, b);
    }
    for (int k = 0; k < 9; k++)
      if (m->rel[b * 9 + k] == 2)
        sexec(T, T->pred[k], bufs, NULL, b);
    if (smooth)
      for (int k = 0; k < 9; k++)
        if (m->rel[b * 9 + k] == 3)
          smoothr(m, b, k, G, W, lab);
    for (int j = -G; j < B + G; j++)
      for (int i = -G; i < B + G; i++) {
        if (i >= 0 && i < B && j >= 0 && j < B)
          continue;
        SV *s = &lab[(j + G) * W + i + G];
        if (s->n < 0) {
          fprintf(stderr, "vm2: error: ghost (%d %d) of block %d not filled\n", i, j, b);
          exit(1);
        }
        if (g->nr + 1 >= cr) {
          cr *= 2;
          g->row = (long *)realloc(g->row, cr * sizeof(long));
          g->p = (int *)realloc(g->p, (cr + 1) * sizeof(int));
        }
        while (g->nnz + s->n >= cn) {
          cn *= 2;
          g->c = (int *)realloc(g->c, cn * sizeof(int));
          g->w = (double *)realloc(g->w, cn * sizeof(double));
        }
        g->row[g->nr] = (long)b * WW + (j + G) * W + i + G;
        g->p[g->nr] = g->nnz;
        for (int e = 0; e < s->n; e++) {
          g->c[g->nnz] = s->i[e];
          g->w[g->nnz++] = s->w[e];
        }
        g->nr++;
      }
  }
  g->p[g->nr] = g->nnz;
  long N = (long)m->nleaf * NB;
  fo->nr = N;
  fo->nnz = g->nnz;
  fo->row = NULL;
  fo->p = (int *)hmalloc((N + 1) * sizeof(int));
  fo->c = (int *)hmalloc((g->nnz + 1) * sizeof(int));
  fo->w = (double *)hmalloc((g->nnz + 1) * sizeof(double));
  for (long n = 0; n <= N; n++)
    fo->p[n] = 0;
  for (long e = 0; e < g->nnz; e++)
    fo->p[g->c[e] + 1]++;
  for (long n = 0; n < N; n++)
    fo->p[n + 1] += fo->p[n];
  int *fill = (int *)hmalloc((N + 1) * sizeof(int));
  memcpy(fill, fo->p, (N + 1) * sizeof(int));
  for (long r = 0; r < g->nr; r++)
    for (int e = g->p[r]; e < g->p[r + 1]; e++) {
      int k = fill[g->c[e]]++;
      fo->c[k] = (int)g->row[r];
      fo->w[k] = g->w[e];
    }
  free(fill);
  free(lab);
  free(buf);
  free(tmp);
}
static double *dup(const void *h, size_t n) {
  double *d = (double *)dmalloc(n);
  up(d, h, n);
  return d;
}
static int *idup(const int *h, size_t n) {
  int *d = (int *)dmalloc(n * sizeof(int));
  up(d, h, n * sizeof(int));
  return d;
}
static void tables(const double2 *F, int Pt, T *t) {
  size_t q = (size_t)Pt * Pt * sizeof(double2);
  t->Pt = Pt;
  t->Tm2m = (double2 *)dmalloc(4 * q);
  t->Tl2l = (double2 *)dmalloc(4 * q);
  t->Tm2l = (double2 *)dmalloc(49 * q);
  t->Tfar = (double2 *)dmalloc(q);
  up(t->Tm2m, F, 4 * q);
  up(t->Tl2l, F + 4L * Pt * Pt, 4 * q);
  up(t->Tm2l, F + 8L * Pt * Pt, 49 * q);
  up(t->Tfar, F + 57L * Pt * Pt, q);
}
static void setup(H *m, T *t, long N, int p, double L, const double2 *F, int Pt) {
  t->N = N;
  t->p = p;
  t->L = L;
  t->nleaf = m->nleaf;
  t->nbox = m->nbox;
  t->Dm = m->Dm;
  t->Fn = 1 << m->Dm;
  int *fine = (int *)hmalloc((size_t)t->Fn * t->Fn * sizeof(int));
  for (int c = 0; c < m->nleaf; c++) {
    int b = m->lbox[c], s = 1 << (m->Dm - m->bd[b]);
    for (int j = 0; j < s; j++)
      for (int i = 0; i < s; i++)
        fine[(long)(m->py[b] * s + j) * t->Fn + m->px[b] * s + i] = c;
  }
  t->fine = idup(fine, (size_t)t->Fn * t->Fn);
  free(fine);
  double *cx = (double *)hmalloc(m->nbox * sizeof(double)), *cy = (double *)hmalloc(m->nbox * sizeof(double));
  for (int b = 0; b < m->nbox; b++) {
    double h = L / (double)(1L << m->bd[b]);
    cx[b] = (m->px[b] + 0.5) * h;
    cy[b] = (m->py[b] + 0.5) * h;
  }
  t->cx = dup(cx, m->nbox * sizeof(double));
  t->cy = dup(cy, m->nbox * sizeof(double));
  free(cx);
  free(cy);
  t->bdep = idup(m->bd, m->nbox);
  t->bch = idup(m->ch, m->nbox);
  t->bpx = idup(m->px, m->nbox);
  t->bpy = idup(m->py, m->nbox);
  t->bpar = idup(m->par, m->nbox);
  t->lbox = idup(m->lbox, m->nleaf);
  t->uo = idup(m->uo, m->nleaf + 1);
  t->ul = idup(m->ul, m->uo[m->nleaf]);
  t->us = idup(m->us, m->uo[m->nleaf]);
  t->wo = idup(m->wo, m->nleaf + 1);
  t->wl = idup(m->wl, m->wo[m->nleaf]);
  t->ws = idup(m->ws, m->wo[m->nleaf]);
  t->xo = idup(m->xo, m->nleaf + 1);
  t->xl = idup(m->xl, m->xo[m->nleaf]);
  t->xs = idup(m->xs, m->xo[m->nleaf]);
  t->ro = idup(m->ro, m->nleaf + 1);
  t->rl = idup(m->rl, m->ro[m->nleaf]);
  t->rs = idup(m->rs, m->ro[m->nleaf]);
  t->vo = idup(m->vo, m->nbox + 1);
  t->vl = idup(m->vl, m->vo[m->nbox]);
  t->vs = idup(m->vs, m->vo[m->nbox]);
  int P = p + 1;
  t->M = (double2 *)dmalloc((size_t)m->nbox * P * sizeof(double2));
  t->Lc = (double2 *)dmalloc((size_t)m->nbox * P * sizeof(double2));
  t->key = (int *)dmalloc(N * sizeof(int));
  t->start = (int *)dmalloc((m->nleaf + 1) * sizeof(int));
  t->perm = (long *)dmalloc(N * sizeof(long));
  t->x = (double *)dmalloc(N * sizeof(double));
  t->y = (double *)dmalloc(N * sizeof(double));
  t->g = (double *)dmalloc(N * sizeof(double));
  t->o = (double *)dmalloc(N * sizeof(double));
  t->oo = (double *)dmalloc(N * sizeof(double));
  t->h = (double *)dmalloc(N * sizeof(double));
  t->v = (double *)dmalloc(N * sizeof(double));
  t->f = (int *)dmalloc(N * sizeof(int));
  t->ox = (double *)dmalloc(N * sizeof(double));
  t->oy = (double *)dmalloc(N * sizeof(double));
  t->og = (double *)dmalloc(N * sizeof(double));
  t->dg = (double *)dmalloc(N * sizeof(double));
  t->odg = (double *)dmalloc(N * sizeof(double));
  t->w = (double2 *)dmalloc(N * sizeof(double2));
  t->ow = (double2 *)dmalloc(N * sizeof(double2));
  tables(F, Pt, t);
}
static long npair, nfpair;
static int graphing, ptime;
static double ptacc[8];
static cudaEvent_t pev[8];
static void pmark(int k) {
  if (ptime && !graphing)
    ck(cudaEventRecord(pev[k], cudaStreamPerThread), "cudaEventRecord");
}
static int nozbar;
static void fmm(T *t, H *m, int mode, const double *x, const double *y,
                const double *g, double *ku, double *kv, double *kg) {
  long N = t->N;
  int P = t->p + 1;
  static int *k2 = NULL, *cnt = NULL;
  static long kN = 0, kl = 0, *perm0 = NULL;
  static void *ctmp = NULL;
  static size_t ctn = 0;
  static double *dsm = NULL;
  static int kbits = 32;
  if (dsm == NULL)
    dsm = (double *)dmalloc(3 * sizeof(double));
  if (kN != N || kl != t->nleaf) {
    if (k2) {
      cudaFree(k2);
      cudaFree(perm0);
      cudaFree(ctmp);
      free(cnt);
    }
    k2 = (int *)dmalloc(N * sizeof(int));
    perm0 = (long *)dmalloc(N * sizeof(long));
    cnt = (int *)hmalloc((t->nleaf + 1) * sizeof(int));
    ctn = 0;
    ctmp = NULL;
    ck(cub::DeviceRadixSort::SortPairs(ctmp, ctn, t->key, k2, perm0, t->perm, (int)N, 0, 32), "SortPairs size");
    kbits = 1;
    while ((1L << kbits) <= t->nleaf)
      kbits++;
    ctmp = dmalloc(ctn);
    kN = N;
    kl = t->nleaf;
  }
  t->mode = 0;
  ck(cudaMemcpyAsync(t->ox, x, N * sizeof(double), cudaMemcpyDeviceToDevice, cudaStreamPerThread), "cudaMemcpyAsync");
  ck(cudaMemcpyAsync(t->oy, y, N * sizeof(double), cudaMemcpyDeviceToDevice, cudaStreamPerThread), "cudaMemcpyAsync");
  ck(cudaMemcpyAsync(t->oo, g, N * sizeof(double), cudaMemcpyDeviceToDevice, cudaStreamPerThread), "cudaMemcpyAsync");
  gam<<<grid(N), nthread>>>(N, g, t->ov, t->og);
  if (mode == 2) {
    fastpair<<<grid(N), nthread>>>(*t, 0);
    launched("fastpair");
    split<<<grid(N), nthread>>>(N, t->ow, t->odg, ku, kv, kg);
    launched("split");
    return;
  }
  pmark(0);
  leafof<<<grid(N), nthread>>>(*t);
  iota<<<grid(N), nthread>>>(N, perm0);
  ck(cub::DeviceRadixSort::SortPairs(ctmp, ctn, t->key, k2, perm0, t->perm, (int)N, 0, kbits, cudaStreamPerThread), "SortPairs");
  ck(cudaMemcpyAsync(t->key, k2, N * sizeof(int), cudaMemcpyDeviceToDevice, cudaStreamPerThread), "cudaMemcpyAsync");
  starts<<<grid(t->nleaf + 1), nthread>>>(N, t->nleaf, k2, t->start);
  gather<<<grid(N), nthread>>>(*t);
  launched("gather");
  if (!graphing) {
    down(cnt, t->start, (t->nleaf + 1) * sizeof(int));
    for (int c = 0; c < t->nleaf; c++) {
      long s2 = 0;
      for (int e = m->uo[c]; e < m->uo[c + 1]; e++)
        s2 += cnt[m->ul[e] + 1] - cnt[m->ul[e]];
      for (int e = m->xo[c]; e < m->xo[c + 1]; e++)
        s2 += cnt[m->xl[e] + 1] - cnt[m->xl[e]];
      for (int e = m->co[c]; e < m->co[c + 1]; e++)
        s2 += cnt[m->cl[e] + 1] - cnt[m->cl[e]];
      npair += s2 * (cnt[c + 1] - cnt[c]);
    }
  }
  pmark(1);
  if (mode != 2) {
    p2m<<<grid((long)t->nleaf * P), nthread>>>(*t);
    launched("p2m");
    for (int d = m->Dm - 1; d >= 0; d--) {
      if (m->hi[d] == m->lo[d])
        continue;
      m2m<<<grid((long)(m->hi[d] - m->lo[d]) * P), nthread>>>(*t, m->lo[d], m->hi[d], d);
      launched("m2m");
    }
    far<<<1, P>>>(*t);
    launched("far");
    m2l<<<grid((long)t->nbox * P), nthread>>>(*t);
    launched("m2l");
    for (int d = 1; d <= m->Dm; d++) {
      if (m->hi[d] == m->lo[d])
        continue;
      l2l<<<grid((long)(m->hi[d] - m->lo[d]) * P), nthread>>>(*t, m->lo[d], m->hi[d], d);
      launched("l2l");
    }
  }
  pmark(2);
  l2p<<<grid(N), nthread>>>(*t);
  launched("l2p");
  pmark(3);
  p2p<<<t->nleaf, nw>>>(*t);
  launched("p2p");
  pmark(4);
  unsort<<<grid(N), nthread>>>(*t);
  launched("unsort");
  if (mode != 2 && !nozbar) {
    sums3<<<1, 1024>>>(N, t->og, t->ox, t->oy, dsm);
    zbar<<<grid(N), nthread>>>(*t, dsm);
    launched("zbar");
  }
  if (mode == 1) {
    fastpair<<<grid(N), nthread>>>(*t, 1);
    launched("fastpair");
  }
  split<<<grid(N), nthread>>>(N, t->ow, t->odg, ku, kv, kg);
  launched("split");
  pmark(5);
  if (ptime && !graphing && mode != 2) {
    ck(cudaEventSynchronize(pev[5]), "sync");
    for (int k = 0; k < 5; k++) {
      float ms;
      ck(cudaEventElapsedTime(&ms, pev[k], pev[k + 1]), "elapsed");
      ptacc[k] += ms;
    }
  }
}
struct St {
  double *x, *y, *g;
};
static St salloc(long N) {
  St s;
  s.x = (double *)dmalloc(N * sizeof(double));
  s.y = (double *)dmalloc(N * sizeof(double));
  s.g = (double *)dmalloc(N * sizeof(double));
  return s;
}
static void scomb(long N, St o, St a, double c1, St k1, double c2, St *k2) {
  comb<<<grid(N), nthread>>>(N, o.x, a.x, c1, k1.x, c2, k2 ? k2->x : NULL);
  comb<<<grid(N), nthread>>>(N, o.y, a.y, c1, k1.y, c2, k2 ? k2->y : NULL);
  comb<<<grid(N), nthread>>>(N, o.g, a.g, c1, k1.g, c2, k2 ? k2->g : NULL);
  launched("comb");
}
static void rhs(T *t, H *m, int mode, St y, St k) {
  fmm(t, m, mode, y.x, y.y, y.g, k.x, k.y, k.g);
}
static long nfull, nslow, nfast;
static void heun(T *t, H *m, St y, double H0, int M, St *w) {
  long N = t->N;
  double h = H0 / M;
  for (int s = 0; s < M; s++) {
    rhs(t, m, 0, y, w[0]);
    scomb(N, w[3], y, h / 3, w[0], 0, NULL);
    rhs(t, m, 0, w[3], w[1]);
    scomb(N, w[3], y, 2 * h / 3, w[1], 0, NULL);
    rhs(t, m, 0, w[3], w[2]);
    scomb(N, y, y, h / 4, w[0], 3 * h / 4, &w[2]);
    nfull += 3;
  }
}
static const double G0[3][3] = {{1. / 3, 0, 0}, {-1. / 3, 2. / 3, 0}, {0, -2. / 3, 1}};
static const double G1[3][3] = {{0, 0, 0}, {0, 0, 0}, {1. / 2, 0, -1. / 2}};
static void force(T *t, H *m, long N, int i, St v, double th, double H0, St *S, St k) {
  rhs(t, m, 2, v, k);
  nfast++;
  scomb(N, k, k, 1. / 3 - 1, k, 0, NULL);
  for (int j = 0; j <= i; j++)
    scomb(N, k, k, G0[i][j] + G1[i][j] * th / H0, S[j], 0, NULL);
}
static void mri(T *t, H *m, St y, double H0, int mi, St *w, St *S) {
  long N = t->N;
  double h = H0 / mi;
  for (int i = 0; i < 3; i++) {
    rhs(t, m, 1, y, S[i]);
    nslow++;
    double th = 0;
    for (int s = 0; s < mi; s++) {
      force(t, m, N, i, y, th, H0, S, w[0]);
      scomb(N, w[3], y, h / 3, w[0], 0, NULL);
      force(t, m, N, i, w[3], th + h / 3, H0, S, w[1]);
      scomb(N, w[3], y, 2 * h / 3, w[1], 0, NULL);
      force(t, m, N, i, w[3], th + 2 * h / 3, H0, S, w[2]);
      scomb(N, y, y, h / 4, w[0], 3 * h / 4, &w[2]);
      th += h;
    }
  }
}
static H m;
static Tab tab;
static CSR gh, fo;
static T t;
static long N;
static int n0, Lm, G, W, nosm, p = 24;
static double L, h0, wband, cd, nu, ramp = 16, band = 16, *dsm2;
static int *dlev, *dbx, *dby, *dgp, *dgc, *hf;
static long *dgr;
static double *dgw, *lab, *lab2, *lab3, *hh, *hx, *hy, *hg, *hu, *hv, *hgr, rgsum;
static St s1, s0, wk[4], S[3], u, dd;
static void dfree(void *a) {
  if (a)
    ck(cudaFree(a), "cudaFree");
}
static void sfree(St s) {
  dfree(s.x);
  dfree(s.y);
  dfree(s.g);
}
static void meshfree(void) {
  void *d[] = {t.fine, t.cx, t.cy, t.bdep, t.bch, t.bpx, t.bpy, t.bpar, t.lbox, t.uo, t.ul, t.us,
               t.wo, t.wl, t.ws, t.xo, t.xl, t.xs, t.ro, t.rl, t.rs, t.vo, t.vl, t.vs, t.M, t.Lc,
               t.key, t.start, t.perm, t.x, t.y, t.g, t.h, t.v, t.f, t.ox, t.oy, t.og, t.dg,
               t.odg, t.w, t.ow, t.Tm2m, t.Tl2l, t.Tm2l, t.Tfar, t.oh, t.ov, t.of, dlev, dbx,
               dby, dgp, dgc, dgr, dgw, lab, lab2, lab3, t.o, t.oo, t.co, t.cl, t.cs};
  for (size_t k = 0; k < sizeof d / sizeof *d; k++)
    dfree(d[k]);
  sfree(s1);
  sfree(s0);
  sfree(u);
  sfree(dd);
  for (int q = 0; q < 4; q++)
    sfree(wk[q]);
  for (int q = 0; q < 3; q++)
    sfree(S[q]);
  free(hh);
  free(hx);
  free(hy);
  free(hg);
  free(hu);
  free(hv);
  free(hgr);
  free(hf);
}
static void hfree(H *a, CSR *g, CSR *f) {
  for (int l = 0; l <= a->Lm; l++)
    free(a->map[l]);
  for (int d = 0; d <= a->Dm; d++)
    free(a->bmap[d]);
  void *h[] = {a->map, a->bmap, a->lev, a->bx, a->by, a->bd, a->px, a->py, a->ch, a->par, a->leaf,
               a->lbox, a->uo, a->ul, a->us, a->vo, a->vl, a->vs, a->wo, a->wl, a->ws, a->xo,
               a->xl, a->xs, a->ro, a->rl, a->rs, a->co, a->cl, a->cs, a->rel, a->nbid, g->row, g->p, g->c, g->w,
               f->p, f->c, f->w};
  for (size_t k = 0; k < sizeof h / sizeof *h; k++)
    free(h[k]);
}
static double bern(double t, int d) {
  double s = 0, c = 1;
  t = t < 0 ? 0 : (t > 1 ? 1 : t);
  for (int j = 0; j <= d; j++) {
    if (2 * j > d)
      s += c * pow(t, j) * pow(1 - t, d - j);
    c = c * (d - j) / (j + 1);
  }
  return s;
}
static double gap(double a, double a0, double sz) {
  double d = a - (a0 + sz / 2);
  d -= L * rint(d / L);
  return fmax(0, fabs(d) - sz / 2);
}
static void pou(void) {
  int nl = m.nleaf, *cand = (int *)hmalloc(nl * sizeof(int));
  double *dmx = (double *)hmalloc(nl * sizeof(double)), Rm = ramp * h0 + (8 + band) * h0 / 2;
  for (int b = 0; b < nl; b++) {
    int l = m.lev[b], nc = 0;
    double hb = h0 / (1 << l), bs = B * hb;
    for (int c = 0; c < nl; c++) {
      if (m.lev[c] >= l)
        continue;
      double hc = h0 / (1 << m.lev[c]), cs = B * hc, ox = m.bx[c] * cs + cs / 2 - (m.bx[b] * bs + bs / 2),
             oy = m.by[c] * cs + cs / 2 - (m.by[b] * bs + bs / 2);
      ox -= L * rint(ox / L);
      oy -= L * rint(oy / L);
      if (fabs(ox) - (cs + bs) / 2 < Rm * hc / h0 && fabs(oy) - (cs + bs) / 2 < Rm * hc / h0)
        cand[nc++] = c;
    }
    dmx[b] = 0;
    for (int q = 0; q < NB; q++) {
      long n = (long)b * NB + q;
      long gi = (long)m.bx[b] * B + q % B, gj = (long)m.by[b] * B + q / B;
      double ph[DM + 2], sk[DM + 2], P = 0, de = h0;
      for (int k = 0; k <= l + 1; k++)
        ph[k] = sk[k] = k <= l;
      for (int e = 0; e < nc; e++) {
        int c = cand[e];
        double hc = h0 / (1 << m.lev[c]), cs = B * hc, gx = gap(hx[n], m.bx[c] * cs, cs),
               gy = gap(hy[n], m.by[c] * cs, cs);
        if (gx >= Rm * hc / h0 || gy >= Rm * hc / h0)
          continue;
        for (int k = m.lev[c] + 1; k <= l; k++) {
          double hk = h0 / (1 << k), wr = ramp * 2 * hk, o = wr + 8 * hk;
          ph[k] *= 1 - (1 - bern(gx / wr, 23)) * (1 - bern(gy / wr, 23));
          sk[k] *= 1 - (1 - bern((gx - o) / (band * hk), 15)) * (1 - bern((gy - o) / (band * hk), 15));
        }
      }
      for (int k = 0; k <= l; k++) {
        long s = 1L << (l - k);
        double hk = h0 / (1 << k);
        if (gi % s == 0 && gj % s == 0)
          P += hk * hk * ph[k] * (1 - ph[k + 1]);
        if (k > 0)
          de += (hk - 2 * hk) * sk[k];
      }
      hh[n] = P;
      hgr[n] = cd * de;
      dmx[b] = fmax(dmx[b], cd * de);
    }
  }
  int nc = 0, cc = 0, *cl = NULL, *cs = NULL;
  m.co = (int *)hmalloc((nl + 1) * sizeof(int));
  for (int c = 0; c < nl; c++) {
    m.co[c] = nc;
    double hc = h0 / (1 << m.lev[c]), sc0 = B * hc, R = 7 * dmx[c];
    for (int e = 0; e < nl; e++) {
      double he = h0 / (1 << m.lev[e]), se = B * he;
      for (int sc = 0; sc < 9; sc++) {
        double ox = m.bx[e] * se + se / 2 + (sc / 3 - 1) * L - (m.bx[c] * sc0 + sc0 / 2),
               oy = m.by[e] * se + se / 2 + (sc % 3 - 1) * L - (m.by[c] * sc0 + sc0 / 2);
        double gx = fmax(0, fabs(ox) - (se + sc0) / 2), gy = fmax(0, fabs(oy) - (se + sc0) / 2);
        if (hypot(gx, gy) >= R)
          continue;
        int in = 0, pse = 1;
        for (int k = m.uo[c]; k < m.uo[c + 1]; k++)
          in |= m.ul[k] == e && m.us[k] == sc;
        for (int k = m.xo[c]; k < m.xo[c + 1]; k++)
          in |= m.xl[k] == e && m.xs[k] == sc;
        for (int k = m.ro[c]; k < m.ro[c + 1]; k++)
          pse &= !(m.rl[k] == e && m.rs[k] == sc);
        if (!in)
          push(&cl, &cs, &nc, &cc, e, sc + 9 * pse);
      }
    }
  }
  m.co[nl] = nc;
  m.cl = cl ? cl : (int *)hmalloc(sizeof(int));
  m.cs = cs ? cs : (int *)hmalloc(sizeof(int));
  t.co = idup(m.co, nl + 1);
  t.cl = idup(m.cl, nc + 1);
  t.cs = idup(m.cs, nc + 1);
  free(cand);
  free(dmx);
}
static void meshup(const int *tgt) {
  build(&m, n0 / B, Lm, wband, tgt);
  compile(&m, &tab, &gh, &fo, !nosm);
  N = (long)m.nleaf * NB;
  memset(&t, 0, sizeof t);
  setup(&m, &t, N, p, L, tab.F, tab.Pt);
  t.cd = cd;
  t.nu = nu;
  dlev = idup(m.lev, m.nleaf);
  dbx = idup(m.bx, m.nleaf);
  dby = idup(m.by, m.nleaf);
  dgr = (long *)dmalloc((gh.nr + 1) * sizeof(long));
  up(dgr, gh.row, gh.nr * sizeof(long));
  dgp = idup(gh.p, gh.nr + 1);
  dgc = idup(gh.c, gh.nnz + 1);
  dgw = dup(gh.w, (gh.nnz + 1) * sizeof(double));
  lab = (double *)dmalloc((size_t)m.nleaf * W * W * sizeof(double));
  lab2 = (double *)dmalloc((size_t)m.nleaf * W * W * sizeof(double));
  lab3 = (double *)dmalloc((size_t)m.nleaf * W * W * sizeof(double));
  s1 = salloc(N);
  s0 = salloc(N);
  u = salloc(N);
  dd = salloc(N);
  for (int q = 0; q < 4; q++)
    wk[q] = salloc(N);
  for (int q = 0; q < 3; q++)
    S[q] = salloc(N);
  hh = (double *)hmalloc(N * sizeof(double));
  hx = (double *)hmalloc(N * sizeof(double));
  hy = (double *)hmalloc(N * sizeof(double));
  hg = (double *)hmalloc(N * sizeof(double));
  hu = (double *)hmalloc(N * sizeof(double));
  hv = (double *)hmalloc(N * sizeof(double));
  hgr = (double *)hmalloc(N * sizeof(double));
  hf = (int *)hmalloc(N * sizeof(int));
  nodes<<<grid(N), nthread>>>(N, h0, dlev, dbx, dby, s0.x, s0.y, t.h);
  launched("nodes");
  down(hx, s0.x, N * sizeof(double));
  down(hy, s0.y, N * sizeof(double));
  pou();
  t.ov = dup(hh, N * sizeof(double));
  t.oh = dup(hgr, N * sizeof(double));
  t.of = (int *)dmalloc(N * sizeof(int));
  ck(cudaMemset(t.of, 0, N * sizeof(int)), "cudaMemset");
}
static void blockmax(const double *d, double *o) {
  double *hb = (double *)hmalloc(N * sizeof(double));
  down(hb, d, N * sizeof(double));
  for (int b = 0; b < m.nleaf; b++) {
    o[b] = 0;
    for (int k = 0; k < NB; k++)
      o[b] = fmax(o[b], hb[(long)b * NB + k]);
  }
  free(hb);
}
static double tval(H *o, CSR *og, const int *lr, const double *wo, int ob, int I, int J) {
  int WW = W * W;
  if (I >= 0 && I < B && J >= 0 && J < B)
    return wo[(long)ob * NB + J * B + I];
  if (I < -G || I >= B + G || J < -G || J >= B + G) {
    fprintf(stderr, "vm2: error: transfer stencil outside the lab\n");
    exit(1);
  }
  int r = lr[(long)ob * WW + (J + G) * W + I + G];
  double s = 0;
  for (int k = og->p[r]; k < og->p[r + 1]; k++)
    s += og->w[k] * wo[og->c[k]];
  (void)o;
  return s;
}
static double pred(H *mo, CSR *go, const int *lr, const double *wo, int ob, int k, long X, long Y) {
  if (k == 0)
    return tval(mo, go, lr, wo, ob, (int)X, (int)Y);
  int ox = (int)(X & 1), oy = (int)(Y & 1);
  double v = 0;
  for (int a = 0; a < (ox ? 8 : 1); a++)
    for (int e = 0; e < (oy ? 8 : 1); e++) {
      long I = ox ? (X - 1) / 2 - 3 + a : X / 2, J = oy ? (Y - 1) / 2 - 3 + e : Y / 2;
      if (X < 0)
        I = ox ? -((-(X - 1) + 1) / 2) - 3 + a : -((-X + 1) / 2);
      if (Y < 0)
        J = oy ? -((-(Y - 1) + 1) / 2) - 3 + e : -((-Y + 1) / 2);
      v += (ox ? DDW[a] : 1) * (oy ? DDW[e] : 1) * pred(mo, go, lr, wo, ob, k - 1, I, J);
    }
  return v;
}
static void regrid(double er, double kap, double tau, int ic, double rho, double eps) {
  long N0 = N;
  int nl0 = m.nleaf, nb0 = n0 / B;
  double *db = (double *)hmalloc(nl0 * sizeof(double)), *sb = (double *)hmalloc(nl0 * sizeof(double));
  ck(cudaMemcpy(wk[0].x, s0.g, N * sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy");
  labown<<<grid(N), nthread>>>(N, W, G, wk[0].x, lab);
  labfill<<<grid(gh.nr), nthread>>>(gh.nr, dgr, dgp, dgc, dgw, wk[0].x, lab);
  detail<<<grid(N), nthread>>>(N, W, G, lab, wk[0].y);
  launched("detail");
  blockmax(wk[0].y, db);
  double wm = 0;
  {
    double *ho = (double *)hmalloc(N * sizeof(double));
    down(ho, wk[0].x, N * sizeof(double));
    for (long n = 0; n < N; n++)
      wm = fmax(wm, fabs(ho[n]));
    free(ho);
  }
  for (int b = 0; b < nl0; b++)
    sb[b] = 0;
  if (kap > 0) {
    rhs(&t, &m, 0, s0, u);
    labown<<<grid(N), nthread>>>(N, W, G, u.x, lab);
    labfill<<<grid(gh.nr), nthread>>>(gh.nr, dgr, dgp, dgc, dgw, u.x, lab);
    labown<<<grid(N), nthread>>>(N, W, G, u.y, lab2);
    labfill<<<grid(gh.nr), nthread>>>(gh.nr, dgr, dgp, dgc, dgw, u.y, lab2);
    strain<<<grid(N), nthread>>>(N, W, G, h0, dlev, lab, lab2, wk[0].g);
    launched("strain");
    blockmax(wk[0].g, sb);
  }
  long nf = (long)nb0 << Lm;
  int *tg = (int *)hmalloc(nf * nf * sizeof(int)), *nlv = (int *)hmalloc(nl0 * sizeof(int)), *fl = (int *)hmalloc(nl0 * sizeof(int));
  double ec = er / 256;
  for (int b = 0; b < nl0; b++) {
    double d = db[b] / wm * exp(kap * sb[b] * tau);
    int l = m.lev[b];
    fl[b] = d > er;
    nlv[b] = d > er && l < Lm ? l + 1 : (d < ec && l > 0 ? l - 1 : l);
  }
  for (int b = 0; b < nl0; b++) {
    if (nlv[b] >= m.lev[b])
      continue;
    int l = m.lev[b], px = m.bx[b] & ~1, py = m.by[b] & ~1, ok = 1;
    for (int q = 0; q < 4; q++) {
      int *c = lmap(&m, l, px + (q & 1), py + (q >> 1));
      ok &= c != NULL && *c >= 0 && nlv[*c] < l;
    }
    if (!ok)
      nlv[b] = l;
  }
  for (long k = 0; k < nf * nf; k++)
    tg[k] = 0;
  for (int pass = 0; pass < 2; pass++)
    for (int b = 0; b < nl0; b++) {
      int l = m.lev[b], s = 1 << (Lm - l), r = pass && fl[b] ? (1 << (Lm - nlv[b])) * nbuf : 0;
      for (long j = (long)m.by[b] * s - r; j < (long)(m.by[b] + 1) * s + r; j++)
        for (long i = (long)m.bx[b] * s - r; i < (long)(m.bx[b] + 1) * s + r; i++) {
          long q = ((j % nf + nf) % nf) * nf + (i % nf + nf) % nf;
          int v = pass ? (fl[b] ? nlv[b] : 0) : nlv[b];
          tg[q] = tg[q] > v ? tg[q] : v;
        }
    }
  H mo = m;
  CSR go = gh, foo = fo;
  double *wo = (double *)hmalloc(N0 * sizeof(double));
  down(wo, wk[0].x, N0 * sizeof(double));
  meshfree();
  meshup(tg);
  for (long n = 0; n < N; n++) {
    if (ic) {
      hg[n] = omega(hx[n], hy[n], rho, eps);
      continue;
    }
    hg[n] = 0;
  }
  if (!ic) {
    int WW = W * W, *lr = (int *)hmalloc((long)mo.nleaf * WW * sizeof(int));
    for (long k = 0; k < (long)mo.nleaf * WW; k++)
      lr[k] = -1;
    for (long r = 0; r < go.nr; r++)
      lr[go.row[r]] = r;
    double gn = 0, vn = 0;
    for (long n = 0; n < N; n++) {
      int b = n / NB, i = n % NB % B, j = n % NB / B, l = m.lev[b], *p;
      long gx = (long)m.bx[b] * B + i, gy = (long)m.by[b] * B + j;
      double v;
      int lo = -1, ob = -1;
      for (int q = 0; q <= mo.Lm && lo < 0; q++) {
        long X = q >= l ? gx << (q - l) : gx >> (l - q), Y = q >= l ? gy << (q - l) : gy >> (l - q);
        if ((p = lmap(&mo, q, X / B, Y / B)) != NULL && *p >= 0) {
          lo = q;
          ob = *p;
        }
      }
      if (lo < 0) {
        fprintf(stderr, "vm2: error: regrid: no old block covers a new node\n");
        exit(1);
      }
      if (lo >= l) {
        long X = gx << (lo - l), Y = gy << (lo - l);
        v = wo[(long)ob * NB + (Y % B) * B + X % B];
      } else {
        int k = l - lo;
        long X = gx - ((long)mo.bx[ob] * B << k), Y = gy - ((long)mo.by[ob] * B << k);
        v = pred(&mo, &go, lr, wo, ob, k, X, Y);
      }
      hg[n] = v;
    }
    free(lr);
    for (long n = 0; n < N; n++) {
      gn += hg[n] * hh[n];
      vn += hh[n];
    }
    double c = (rgsum - gn) / vn;
    for (long n = 0; n < N; n++)
      hg[n] += c;
  }
  up(s0.g, hg, N * sizeof(double));
  hfree(&mo, &go, &foo);
  free(wo);
  free(tg);
  free(nlv);
  free(fl);
  free(db);
  free(sb);
}
int main(int argc, char **argv) {
  FILE *file;
  char *end, *rpath = NULL, *tpath = NULL, *wpath = NULL;
  double x, Tend = -1, H0 = -1, rho = pi / 15, eps = 0.05, D = -1, fr = 2, er = 1e-6, kap = 0;
  int M = 1, lts = 0, ord = 8, pord = 8, chk = 0, rn = 0, verb = 0, tst = 0, rg = 0, gr = 0;
  nu = -1;
  cd = 2;
  L = 2 * pi;
  (void)argc;
  memset(&t, 0, sizeof t);
  while (*++argv != NULL && argv[0][0] == '-') {
    if (argv[0][1] == 'h') {
      fprintf(stderr,
              "Usage: vm2 -n <base nodes> -t <end time> -s <macro step> -u <viscosity> "
              "-d <output interval> [-L <levels> -w <band>] [-M <substeps>] [-l (MRI)] "
              "[-f <fast fraction>] [-p <terms>] [-o <core order>] [-O <pse order>] "
              "[-k <core width / h>] [-a <ramp, coarse cells>] [-A <core band, cells>] [-r <rho>] [-e <eps>] [-g <table>] "
              "[-i <reference> -I <reference N>] [-c (velocity check)]\n");
      exit(1);
    }
    if (argv[0][1] == 'l') {
      lts = 1;
      continue;
    }
    if (argv[0][1] == 'c') {
      chk = 1;
      continue;
    }
    if (argv[0][1] == 'v') {
      verb = 1;
      continue;
    }
    if (argv[0][1] == 'j') {
      nosm = 1;
      continue;
    }
    if (argv[0][1] == 'T') {
      tst = 1;
      continue;
    }
    if (argv[0][1] == 'P') {
      pair_ = 1;
      continue;
    }
    if (argv[0][1] == 'G' && argv[0][2] == 'r') {
      gr = 1;
      continue;
    }
    if (argv[1] == NULL) {
      fprintf(stderr, "vm2: error: %s needs an argument\n", argv[0]);
      exit(1);
    }
    if (argv[0][1] == 'i') {
      rpath = *++argv;
      continue;
    }
    if (argv[0][1] == 'g') {
      tpath = *++argv;
      continue;
    }
    if (argv[0][1] == 'W') {
      wpath = *++argv;
      continue;
    }
    x = strtod(argv[1], &end);
    if (*end != '\0') {
      fprintf(stderr, "vm2: error: '%s' is not a number\n", argv[1]);
      exit(1);
    }
    switch (argv[0][1]) {
    case 'n':
      n0 = x;
      break;
    case 'L':
      Lm = x;
      break;
    case 'w':
      wband = x;
      break;
    case 'R':
      rg = x;
      break;
    case 'E':
      er = x;
      break;
    case 'K':
      kap = x;
      break;
    case 't':
      Tend = x;
      break;
    case 's':
      H0 = x;
      break;
    case 'M':
      M = x;
      break;
    case 'f':
      fr = x;
      break;
    case 'u':
      nu = x;
      break;
    case 'd':
      D = x;
      break;
    case 'p':
      p = x;
      break;
    case 'o':
      ord = x;
      break;
    case 'O':
      pord = x;
      break;
    case 'k':
      cd = x;
      break;
    case 'a':
      ramp = x;
      break;
    case 'A':
      band = x;
      break;
    case 'r':
      rho = x;
      break;
    case 'e':
      eps = x;
      break;
    case 'I':
      rn = x;
      break;
    default:
      fprintf(stderr, "vm2: error: unknown option '%s'\n", *argv);
      exit(1);
    }
    argv++;
  }
  if (n0 < B || n0 % B || Tend <= 0 || H0 <= 0 || nu < 0 || D <= 0 || M < 1 || p < 4 ||
      p > PM || ord % 2 || ord > 16 || pord < 2 || pord % 2 || pord > 16 ||
      (lts && M % 3)) {
    fprintf(stderr, "vm2: error: need -n (multiple of %d) -t -s -u -d; -M >= 1 "
                    "(multiple of 3 with -l); 4 <= -p <= %d; even -o, -O <= 16\n",
            B, PM);
    exit(1);
  }
  long nst = steps("-t", Tend, H0), nd = steps("-d", D, H0);
  nbuf = (int)ceil((2 * ramp + 8 + band) / B);
  tabload(&tab, tpath ? tpath : "vm2_g8.bin");
  double R[8] = {0}, E[8] = {0}, Rr[8] = {0}, Er[8] = {0}, Ex[8] = {0};
  int nr = ord / 2, ne = pord / 2, nrr, ner;
  if (ord > 8 || pord > 8) {
    fprintf(stderr, "vm2: error: -o, -O <= 8\n");
    exit(1);
  }
  for (int k = 0; k < 8; k++) {
    R[k] = ord > 0 ? tab.R[ord / 2 - 1][k] : 0;
    E[k] = tab.E[pord / 2 - 1][k];
  }
  cores(ord > 0 ? ord : 2, Rr, &nrr, Ex, &ner);
  cores(pord, Ex, &nrr, Er, &ner);
  if (verb) {
    double dr = 0, de = 0;
    for (int k = 0; k < 8; k++) {
      dr = fmax(dr, ord > 0 ? fabs(R[k] - Rr[k]) : 0);
      de = fmax(de, fabs(E[k] - Er[k]) / fabs(E[0]));
    }
    fprintf(stderr, "vm2: core/pse tables (mpmath) vs double solve: %.3e %.3e\n", dr, de);
  }
  if (nr > 0 && R[0] != 1) {
    fprintf(stderr, "vm2: error: core table R[0] = %.17g, kcore needs exactly 1\n", R[0]);
    exit(1);
  }
  ck(cudaMemcpyToSymbol(cR, R, sizeof R), "cudaMemcpyToSymbol");
  ck(cudaMemcpyToSymbol(nR, &nr, sizeof nr), "cudaMemcpyToSymbol");
  ck(cudaMemcpyToSymbol(cE, E, sizeof E), "cudaMemcpyToSymbol");
  ck(cudaMemcpyToSymbol(nE, &ne, sizeof ne), "cudaMemcpyToSymbol");
  ck(cudaMemcpyToSymbol(cK, K8, sizeof K8), "cudaMemcpyToSymbol");
  G = tab.hd[2];
  W = tab.hd[4];
  ptime = getenv("VM2TIME") != NULL;
  for (int k = 0; k < 8; k++)
    ck(cudaEventCreate(&pev[k]), "cudaEventCreate");
  h0 = L / n0;
  if (p + 1 > tab.Pt) {
    fprintf(stderr, "vm2: error: -p %d exceeds the table (%d)\n", p, tab.Pt - 1);
    exit(1);
  }
  if (rg) {
    long nf = (long)(n0 / B) << Lm;
    int *tg = (int *)hmalloc(nf * nf * sizeof(int));
    for (long k = 0; k < nf * nf; k++)
      tg[k] = 0;
    meshup(tg);
    free(tg);
  } else
    meshup(NULL);
  for (long n = 0; n < N; n++)
    hg[n] = omega(hx[n], hy[n], rho, eps);
  up(s0.g, hg, N * sizeof(double));
  dsm2 = (double *)dmalloc(3 * sizeof(double));
  if (rg)
    for (int k = 0; k <= Lm; k++)
      regrid(er, kap, rg * H0, 1, rho, eps);
  down(hg, s0.g, N * sizeof(double));
  rgsum = 0;
  for (long n = 0; n < N; n++)
    rgsum += hg[n] * hh[n];
  if (tst) {
    int WW = W * W;
    double em = 0, fm = 0, ec = 0;
    for (int a = 0; a <= 7; a++)
      for (int e = 0; e <= 7; e++) {
        for (long r = 0; r < gh.nr; r++) {
          int b = gh.row[r] / WW, c = gh.row[r] % WW, gi = c % W - G, gj = c / W - G;
          double hb = h0 / (1 << m.lev[b]), Dd = L / hb, s2 = 0;
          for (int k = gh.p[r]; k < gh.p[r + 1]; k++) {
            int n = gh.c[k], nb = n / NB, li = n % NB % B, lj = n % NB / B;
            double hn = h0 / (1 << m.lev[nb]);
            double xx = ((m.bx[nb] * B + li) * hn - m.bx[b] * B * hb) / hb,
                   yy = ((m.by[nb] * B + lj) * hn - m.by[b] * B * hb) / hb;
            xx -= Dd * rint((xx - gi) / Dd);
            yy -= Dd * rint((yy - gj) / Dd);
            s2 += gh.w[k] * pow(xx / 16, a) * pow(yy / 16, e);
          }
          double f = pow(gi / 16.0, a) * pow(gj / 16.0, e);
          em = fmax(em, fabs(s2 - f));
          fm = fmax(fm, fabs(f));
        }
      }
    for (long r = 0; r < gh.nr; r++) {
      double s2 = 0;
      for (int k = gh.p[r]; k < gh.p[r + 1]; k++)
        s2 += gh.w[k];
      ec = fmax(ec, fabs(s2 - 1));
    }
    printf("ghost rows %ld nnz %ld: tensor degrees <= 7 max error %.3e (max |f| %.3e), max |row sum - 1| %.3e\n",
           gh.nr, gh.nnz, em, fm, ec);
    return 0;
  }
  {
    double a = 0, mn = 1e300, mx = 0, qa = 0, qb = 0;
    for (long n = 0; n < N; n++) {
      a += hh[n];
      double r = hh[n] * (n0 << m.lev[n / NB]) * (n0 << m.lev[n / NB]) / (L * L);
      mn = fmin(mn, r);
      mx = fmax(mx, r);
      qa += hh[n] * cos(hx[n] + 0.3) * cos(hx[n] + 0.3) * sin(hy[n]) * sin(hy[n]);
      qb += hh[n] * exp(sin(hx[n] + 0.1)) * cos(hy[n] - 0.2) * cos(hy[n] - 0.2);
    }
    fprintf(stderr, "vm2: sum of volumes / area - 1 = %.3e, V/h^2 in [%.6f, %.6f], quadrature errors %.3e %.3e\n",
            a / (L * L) - 1, mn, mx, qa / (pi * pi) - 1, qb / (2 * pi * pi * 1.2660658777520083) - 1);
  }
  fprintf(stderr, "vm2: blocks %d (levels %d) particles %ld boxes %d ghost rows %ld nnz %ld U %d V %d W %d X %d\n",
          m.nleaf, Lm + 1, N, m.nbox, gh.nr, gh.nnz, m.uo[m.nleaf], m.vo[m.nbox], m.wo[m.nleaf], m.xo[m.nleaf]);
  if (chk) {
    double t0 = now();
    rhs(&t, &m, 0, s0, u);
    ck(cudaDeviceSynchronize(), "sync");
    double t1 = now();
    down(hu, u.x, N * sizeof(double));
    down(hv, u.y, N * sizeof(double));
    int n = rpath ? rn : n0, nh = n / 2 + 1;
    if (!rpath && Lm != 0) {
      fprintf(stderr, "vm2: error: -c on an adaptive mesh needs -i -I (vorticity on a fine grid)\n");
      exit(1);
    }
    if (n % (n0 << Lm)) {
      fprintf(stderr, "vm2: error: -I must be a multiple of %d\n", n0 << Lm);
      exit(1);
    }
    double *dr, *gr = (double *)hmalloc((long)n * n * sizeof(double));
    double2 *c0, *c1, *c2;
    cufftHandle fw, bw;
    dr = (double *)dmalloc((long)n * n * sizeof(double));
    c0 = (double2 *)dmalloc((long)n * nh * sizeof(double2));
    c1 = (double2 *)dmalloc((long)n * nh * sizeof(double2));
    c2 = (double2 *)dmalloc((long)n * nh * sizeof(double2));
    cufft(cufftPlan2d(&fw, n, n, CUFFT_D2Z), "cufftPlan2d");
    cufft(cufftPlan2d(&bw, n, n, CUFFT_Z2D), "cufftPlan2d");
    if (rpath) {
      if ((file = fopen(rpath, "r")) == NULL || fread(gr, sizeof(double), (long)n * n, file) != (size_t)n * n) {
        fprintf(stderr, "vm2: error: fail to read '%s'\n", rpath);
        exit(1);
      }
      fclose(file);
    } else
      for (long q = 0; q < N; q++) {
        int b = q / NB, i = q % NB % B, j = q % NB / B;
        gr[(long)(m.by[b] * B + j) * n + m.bx[b] * B + i] = hg[q];
      }
    up(dr, gr, (long)n * n * sizeof(double));
    cufft(cufftExecD2Z(fw, dr, c0), "cufftExecD2Z");
    spread<<<grid((long)n * nh), nthread>>>(n, c0, c1, c2);
    launched("spread");
    double eu = 0, um = 0, eb = 0;
    for (int q = 0; q < 2; q++) {
      cufft(cufftExecZ2D(bw, q ? c2 : c1, dr), "cufftExecZ2D");
      down(gr, dr, (long)n * n * sizeof(double));
      for (long k = 0; k < N; k++) {
        int b = k / NB, i = k % NB % B, j = k % NB / B, s = n / (n0 << m.lev[b]);
        double f = gr[(long)(m.by[b] * B + j) * s * n + (long)(m.bx[b] * B + i) * s] / ((double)n * n),
               v = q ? hv[k] : hu[k];
        eu = fmax(eu, fabs(v - f));
        um = fmax(um, fabs(f));
        int bd = 0;
        for (int kk = 0; kk < 9; kk++)
          bd |= m.rel[b * 9 + kk] != 1 && kk != 4;
        if (bd)
          eb = fmax(eb, fabs(v - f));
      }
    }
    fprintf(stderr, "vm2: max velocity error in blocks next to a level change %.3e\n", eb / um);
    {
      nozbar = 1;
      rhs(&t, &m, 0, s0, u);
      nozbar = 0;
      long mm = 4096 < N ? 4096 : N;
      double *dout = (double *)dmalloc(4 * mm * sizeof(double)), *ho = (double *)hmalloc(4 * mm * sizeof(double));
      direct<<<mm, nw>>>(t, mm, t.Lc, dout);
      launched("direct");
      down(ho, dout, 4 * mm * sizeof(double));
      double e = 0, mx = 0;
      long wi = 0;
      for (long q = 0; q < mm; q++) {
        double d = hypot(ho[4 * q] - ho[4 * q + 2], ho[4 * q + 1] - ho[4 * q + 3]);
        if (d > e) {
          e = d;
          wi = q * (N / mm);
        }
        mx = fmax(mx, hypot(ho[4 * q], ho[4 * q + 1]));
      }
      fprintf(stderr, "vm2: FMM vs direct (3x3 images + root lattice local), %ld targets: max rel %.3e at particle %ld block %ld lev %d\n",
              mm, e / mx, wi, wi / NB, m.lev[wi / NB]);
    }
    printf("n p order k maxerr maxu rel time\n");
    printf("%d %d %d % .16e % .16e % .16e % .16e % .16e\n", n0, p, ord, cd, eu, um, eu / um, t1 - t0);
    return 0;
  }
  if (chk == 0 && lts && fr < -1) {
    for (long n = 0; n < N; n++)
      hf[n] = (n / NB) % 3 == 0;
    up(t.of, hf, N * sizeof(int));
    St a = salloc(N), b2 = salloc(N);
    double *pa[3] = {a.x, a.y, a.g}, *pb[3] = {b2.x, b2.y, b2.g}, *p0[3] = {s0.x, s0.y, s0.g};
    for (int q = 0; q < 3; q++) {
      ck(cudaMemcpy(pa[q], p0[q], N * sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy");
      ck(cudaMemcpy(pb[q], p0[q], N * sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy");
    }
    heun(&t, &m, a, H0, 1, wk);
    {
      St f0 = salloc(N), f1 = salloc(N), f2 = salloc(N);
      rhs(&t, &m, 0, a, f0);
      rhs(&t, &m, 1, a, f1);
      rhs(&t, &m, 2, a, f2);
      double *A0 = (double *)hmalloc(N * sizeof(double)), *A1 = (double *)hmalloc(N * sizeof(double)),
             *A2 = (double *)hmalloc(N * sizeof(double));
      down(A0, f0.x, N * sizeof(double));
      down(A1, f1.x, N * sizeof(double));
      down(A2, f2.x, N * sizeof(double));
      double e = 0, mf = 0;
      for (long n = 0; n < N; n++) {
        e = fmax(e, fabs(A0[n] - A1[n] - A2[n]));
        mf = fmax(mf, fabs(A2[n]));
      }
      fprintf(stderr, "vm2: moved state: max |full - slow - fast| %.3e, max |fast| %.3e\n", e, mf);
    }
    for (int q = 0; q < 3; q++)
      ck(cudaMemcpy(pa[q], p0[q], N * sizeof(double), cudaMemcpyDeviceToDevice), "cudaMemcpy");
    heun(&t, &m, a, H0, 1, wk);
    mri(&t, &m, b2, H0, M / 3, wk, S);
    double *A = (double *)hmalloc(N * sizeof(double)), *Bb = (double *)hmalloc(N * sizeof(double)),
           *C0 = (double *)hmalloc(N * sizeof(double));
    for (int q = 0; q < 3; q++) {
      down(A, pa[q], N * sizeof(double));
      down(Bb, pb[q], N * sizeof(double));
      down(C0, p0[q], N * sizeof(double));
      double ef = 0, es = 0, mv = 0;
      for (long n = 0; n < N; n++) {
        double d = fabs(A[n] - Bb[n]);
        if (hf[n])
          ef = fmax(ef, d);
        else
          es = fmax(es, d);
        mv = fmax(mv, fabs(A[n] - C0[n]));
      }
      fprintf(stderr, "vm2: one step H %g comp %d: max |heun - mri| fast %.3e slow %.3e, max increment %.3e\n",
              H0, q, ef, es, mv);
    }
    return 0;
  }
  if (chk == 0 && lts && fr < 0) {
    for (long n = 0; n < N; n++)
      hf[n] = (n / NB) % 3 == 0;
    up(t.of, hf, N * sizeof(int));
    St a = wk[0], b2 = wk[1], c = wk[2];
    rhs(&t, &m, 0, s0, a);
    rhs(&t, &m, 1, s0, b2);
    rhs(&t, &m, 2, s0, c);
    double *A = (double *)hmalloc(N * sizeof(double)), *Bb = (double *)hmalloc(N * sizeof(double)),
           *Cc = (double *)hmalloc(N * sizeof(double));
    double *pa[3] = {a.x, a.y, a.g}, *pb[3] = {b2.x, b2.y, b2.g}, *pc[3] = {c.x, c.y, c.g};
    for (int q = 0; q < 3; q++) {
      down(A, pa[q], N * sizeof(double));
      down(Bb, pb[q], N * sizeof(double));
      down(Cc, pc[q], N * sizeof(double));
      double e = 0, mx = 0, ef = 0;
      for (long n = 0; n < N; n++) {
        e = fmax(e, fabs(A[n] - Bb[n] - Cc[n]));
        mx = fmax(mx, fabs(A[n]));
        if (!hf[n])
          ef = fmax(ef, fabs(Cc[n]));
      }
      fprintf(stderr, "vm2: split check component %d: max |full - slow - fast| %.3e (max |full| %.3e), max |fast| on slow particles %.3e\n",
              q, e, mx, ef);
    }
    return 0;
  }
  if (verb) {
    rhs(&t, &m, 0, s0, u);
    int P = p + 1;
    double2 *hM = (double2 *)hmalloc((size_t)m.nbox * P * sizeof(double2)),
            *hL = (double2 *)hmalloc((size_t)m.nbox * P * sizeof(double2));
    down(hM, t.M, (size_t)m.nbox * P * sizeof(double2));
    down(hL, t.Lc, (size_t)m.nbox * P * sizeof(double2));
    for (int d = 0; d <= m.Dm; d++) {
      double a0 = 1e300, a1 = 0, b0 = 1e300, b1 = 0;
      for (int b = m.lo[d]; b < m.hi[d]; b++)
        for (int k = 0; k < P; k++) {
          double x1 = hypot(hM[(long)b * P + k].x, hM[(long)b * P + k].y),
                 x2 = hypot(hL[(long)b * P + k].x, hL[(long)b * P + k].y);
          if (x1 > 0)
            a0 = fmin(a0, x1), a1 = fmax(a1, x1);
          if (x2 > 0)
            b0 = fmin(b0, x2), b1 = fmax(b1, x2);
        }
      fprintf(stderr, "vm2: depth %d boxes %d |M_k| in [%.1e, %.1e] |L_k| in [%.1e, %.1e]\n", d,
              m.hi[d] - m.lo[d], a0, a1, b0, b1);
    }
    down(hg, s0.g, N * sizeof(double));
    for (int l = 0; l <= Lm; l++) {
      double g0 = 1e300, g1 = 0, v0 = 1e300, v1 = 0, o1 = 0;
      for (long n = 0; n < N; n++)
        if (m.lev[n / NB] == l) {
          g0 = fmin(g0, fabs(hg[n] * hh[n]) > 0 ? fabs(hg[n] * hh[n]) : g0);
          g1 = fmax(g1, fabs(hg[n] * hh[n]));
          v0 = fmin(v0, hh[n]);
          v1 = fmax(v1, hh[n]);
          o1 = fmax(o1, fabs(hg[n]));
        }
      fprintf(stderr, "vm2: level %d |Gamma| in [%.1e, %.1e] V in [%.3e, %.3e] max|omega| %.3e\n", l, g0, g1, v0, v1, o1);
    }
    free(hM);
    free(hL);
  }
  double wt = now();
  long cost0 = 0;
  long nreg = 0, gdf = 0, gds = 0, gda = 0;
  double nodet = 0, *dmx = (double *)dmalloc(sizeof(double));
  cudaGraphExec_t gexe = NULL;
  printf("step t circ circabs Z E wmax nfast nfull nslow nfast_evals pairs fastpairs wall nodes\n");
  int nfst = 0;
  long pf = 0;
  for (long st = 0;; st++) {
    if (st % nd == 0 || lts) {
      rhs(&t, &m, 0, s0, u);
      nfull++;
    }
    if (st % nd == 0) {
      down(hg, s0.g, N * sizeof(double));
      down(hu, u.x, N * sizeof(double));
      down(hv, u.y, N * sizeof(double));
      double c = 0, ca = 0, z = 0, e = 0, wm = 0;
      for (long n = 0; n < N; n++) {
        double a = hh[n], o = hg[n];
        c += o * a;
        ca += fabs(o * a);
        z += 0.5 * o * o * a;
        e += 0.5 * (hu[n] * hu[n] + hv[n] * hv[n]) * a;
        wm = fmax(wm, fabs(o));
      }
      printf("% 10ld % .16e % .16e % .16e % .16e % .16e % .16e %d %ld %ld %ld %ld %ld % .16e %ld\n", st,
             st * H0, c, ca, z, e, wm, nfst, nfull, nslow, nfast, npair, nfpair, now() - wt, N);
      if (fflush(stdout) != 0 || ferror(stdout)) {
        fprintf(stderr, "vm2: error: fail to write stdout\n");
        exit(1);
      }
    }
    if (st >= nst)
      break;
    if (rg && st > 0 && st % rg == 0) {
      regrid(er, kap, rg * H0, 0, rho, eps);
      nreg++;
      if (gexe != NULL) {
        ck(cudaGraphExecDestroy(gexe), "GraphExecDestroy");
        gexe = NULL;
      }
    }
    nodet += N;
    if (lts) {
      labown<<<grid(N), nthread>>>(N, W, G, u.x, lab);
      labfill<<<grid(gh.nr), nthread>>>(gh.nr, dgr, dgp, dgc, dgw, u.x, lab);
      labown<<<grid(N), nthread>>>(N, W, G, u.y, lab2);
      labfill<<<grid(gh.nr), nthread>>>(gh.nr, dgr, dgp, dgc, dgw, u.y, lab2);
      grad<<<grid(N), nthread>>>(N, W, G, h0, dlev, lab, lab2, S[0].x);
      launched("grad");
      down(hgr, S[0].x, N * sizeof(double));
      double gm = 0;
      for (long n = 0; n < N; n++)
        gm = fmax(gm, hgr[n]);
      int *fb = (int *)hmalloc(m.nleaf * sizeof(int)), *fb2 = (int *)hmalloc(m.nleaf * sizeof(int));
      for (int b = 0; b < m.nleaf; b++) {
        double q = 0;
        for (int k = 0; k < NB; k++)
          q = fmax(q, hgr[(long)b * NB + k]);
        fb[b] = q > fr * gm;
      }
      for (int b = 0; b < m.nleaf; b++) {
        fb2[b] = fb[b];
        for (int e = m.uo[b]; e < m.uo[b + 1]; e++)
          fb2[b] |= fb[m.ul[e]];
      }
      nfst = 0;
      for (long n = 0; n < N; n++) {
        hf[n] = fb2[n / NB];
        nfst += hf[n];
      }
      pf = 0;
      for (int b = 0; b < m.nleaf; b++)
        if (fb2[b])
          for (int e = m.uo[b]; e < m.uo[b + 1]; e++)
            pf += (long)NB * NB * fb2[m.ul[e]];
      up(t.of, hf, N * sizeof(int));
      free(fb);
      free(fb2);
    }
    long c0f = nfull, c0s = nslow, c0a = nfast, c0p = npair;
    if (gr && gexe == NULL) {
      ck(cudaStreamBeginCapture(cudaStreamPerThread, cudaStreamCaptureModeThreadLocal), "BeginCapture");
      graphing = 1;
    }
    if (!gr || graphing) {
      ck(cudaMemcpyAsync(s1.x, s0.x, N * sizeof(double), cudaMemcpyDeviceToDevice, cudaStreamPerThread), "cudaMemcpyAsync");
      ck(cudaMemcpyAsync(s1.y, s0.y, N * sizeof(double), cudaMemcpyDeviceToDevice, cudaStreamPerThread), "cudaMemcpyAsync");
      ck(cudaMemcpyAsync(s1.g, s0.g, N * sizeof(double), cudaMemcpyDeviceToDevice, cudaStreamPerThread), "cudaMemcpyAsync");
      if (lts) {
        long f0 = nfast;
        mri(&t, &m, s1, H0, M / 3, wk, S);
        nfpair += pf * (nfast - f0);
      } else
        heun(&t, &m, s1, H0, M, wk);
      disp<<<grid(N), nthread>>>(N, L, s1.x, s1.y, s0.x, s0.y, h0, dlev, S[0].x);
      dmax1<<<1, 1024>>>(N, S[0].x, dmx);
      dsp<<<grid(N), nthread>>>(N, L, s1.x, s1.y, s0.x, s0.y, dd.x, dd.y);
      labown<<<grid(N), nthread>>>(N, W, G, s1.g, lab);
      labfill<<<grid(gh.nr), nthread>>>(gh.nr, dgr, dgp, dgc, dgw, s1.g, lab);
      labown<<<grid(N), nthread>>>(N, W, G, dd.x, lab2);
      labfill<<<grid(gh.nr), nthread>>>(gh.nr, dgr, dgp, dgc, dgw, dd.x, lab2);
      labown<<<grid(N), nthread>>>(N, W, G, dd.y, lab3);
      labfill<<<grid(gh.nr), nthread>>>(gh.nr, dgr, dgp, dgc, dgw, dd.y, lab3);
      remesh<<<grid(N), nthread>>>(N, W, G, h0, dlev, lab, lab2, lab3, s0.g);
      sums3<<<1, 1024>>>(N, t.ov, s0.g, s0.g, dsm2);
      gfix<<<grid(N), nthread>>>(N, rgsum, dsm2, s0.g);
      launched("step");
    }
    if (graphing) {
      cudaGraph_t gph;
      ck(cudaStreamEndCapture(cudaStreamPerThread, &gph), "EndCapture");
      ck(cudaGraphInstantiate(&gexe, gph, NULL, NULL, 0), "GraphInstantiate");
      ck(cudaGraphDestroy(gph), "GraphDestroy");
      graphing = 0;
      gdf = nfull - c0f;
      gds = nslow - c0s;
      gda = nfast - c0a;
      nfull = c0f;
      nslow = c0s;
      nfast = c0a;
      npair = c0p;
    }
    if (gr) {
      ck(cudaGraphLaunch(gexe, cudaStreamPerThread), "GraphLaunch");
      nfull += gdf;
      nslow += gds;
      nfast += gda;
    }
    double dm;
    ck(cudaMemcpy(&dm, dmx, sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    if (dm > G - 4) {
      fprintf(stderr, "vm2: error: step %ld: displacement %g cells exceeds ghost width %d - 4\n", st, dm, G);
      exit(1);
    }
  }
  (void)cost0;
  if (ptime)
    fprintf(stderr, "vm2: time [s] sort+gather %.3f fmm-up+m2l+down %.3f l2p+m2p %.3f p2p %.3f unsort+zbar %.3f\n",
            ptacc[0] / 1e3, ptacc[1] / 1e3, ptacc[2] / 1e3, ptacc[3] / 1e3, ptacc[4] / 1e3);
  if (wpath != NULL) {
    down(hg, s0.g, N * sizeof(double));
    if ((file = fopen(wpath, "w")) == NULL || fwrite(hg, sizeof(double), N, file) != (size_t)N ||
        fclose(file) != 0) {
      fprintf(stderr, "vm2: error: fail to write '%s'\n", wpath);
      exit(1);
    }
  }
  if (rpath != NULL) {
    if (rn <= 0 || rn % (n0 << Lm)) {
      fprintf(stderr, "vm2: error: -I must be a multiple of %d\n", n0 << Lm);
      exit(1);
    }
    double *ref = (double *)hmalloc((long)rn * rn * sizeof(double));
    if ((file = fopen(rpath, "r")) == NULL) {
      fprintf(stderr, "vm2: error: fail to open '%s'\n", rpath);
      exit(1);
    }
    if (fread(ref, sizeof(double), (long)rn * rn, file) != (size_t)rn * rn) {
      fprintf(stderr, "vm2: error: fail to read '%s'\n", rpath);
      exit(1);
    }
    fclose(file);
    down(hg, s0.g, N * sizeof(double));
    double a2 = 0, b2 = 0, ai = 0, bi = 0, el[DM + 1] = {0}, eg = 0;
    for (long n = 0; n < N; n++) {
      int b = n / NB, i = n % NB % B, j = n % NB / B, s = rn / (n0 << m.lev[b]);
      double f = ref[(long)(m.by[b] * B + j) * s * rn + (long)(m.bx[b] * B + i) * s],
             o = hg[n], a = hh[n];
      a2 += (o - f) * (o - f) * a;
      b2 += f * f * a;
      ai = fmax(ai, fabs(o - f));
      bi = fmax(bi, fabs(f));
      el[m.lev[b]] = fmax(el[m.lev[b]], fabs(o - f));
      int bd = 0;
      for (int k = 0; k < 9; k++)
        bd |= m.rel[b * 9 + k] != 1 && k != 4;
      if (bd)
        eg = fmax(eg, fabs(o - f));
    }
    for (int l = 0; l <= Lm; l++) {
      int nbl = 0;
      for (int c = 0; c < m.nleaf; c++)
        nbl += m.lev[c] == l;
      fprintf(stderr, "vm2: level %d blocks %d Linf %.3e\n", l, nbl, el[l] / bi);
    }
    fprintf(stderr, "vm2: Linf in blocks next to a level change %.3e\n", eg / bi);
    printf("t L2 Linf nodes_mean regrids\n% .16e % .16e % .16e % .16e %ld\n", nst * H0, sqrt(a2 / b2), ai / bi,
           nst > 0 ? nodet / nst : (double)N, nreg);
  }
  return 0;
}
