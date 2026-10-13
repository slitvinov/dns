#include <cuda_profiler_api.h>
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
static unsigned grid2(long n, int nt) { return (n + nt - 1) / nt; }
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
__host__ __device__ static double omega(double x, double y, double rho, double eps, int pr) {
  const double q = 3.141592653589793238;
  if (pr) {
    double s = -2.0 / (4 * q * q);
    for (int k = 0; k < 2; k++)
      for (int a = -1; a <= 1; a++)
        for (int b = -1; b <= 1; b++) {
          double dx = x - q - (k ? eps : -eps) + 2 * q * a, dy = y - q + 2 * q * b;
          s += exp(-(dx * dx + dy * dy) / (rho * rho)) / (q * rho * rho);
        }
    return s;
  }
  double s = eps * cos(x);
  for (int n = -3; n <= 3; n++) {
    double a = cosh((y - q / 2 + 2 * q * n) / rho),
           b = cosh((y - 3 * q / 2 + 2 * q * n) / rho);
    s -= (1 / (a * a) - 1 / (b * b)) / rho;
  }
  return s;
}
__global__ static void icw(long N, int pr, double rho, double eps, const double *x, const double *y, double *o) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n < N)
    o[n] = omega(x[n], y[n], rho, eps, pr);
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
    s = fma(gw[e], v[gc[e]], s);
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
/* Exactly rounded double operations (no FMA contraction) for everything the
   host reference and the GPU regrid must produce bit-identically: node weights,
   core widths, the old -> new value transfer. Host code is built with
   -ffp-contract=off. */
#ifdef __CUDA_ARCH__
#define XM(a, b) __dmul_rn((a), (b))
#define XA(a, b) __dadd_rn((a), (b))
#define XS(a, b) __dsub_rn((a), (b))
#define XD(a, b) __ddiv_rn((a), (b))
#else
#define XM(a, b) ((double)(a) * (double)(b))
#define XA(a, b) ((double)(a) + (double)(b))
#define XS(a, b) ((double)(a) - (double)(b))
#define XD(a, b) ((double)(a) / (double)(b))
#endif
enum { NS = 45, NCODE = 4 * 6561, UL1 = 128, SK = 128, CL = 1024, PC = 512 };
struct Off {
  long o[DM + 3];
};
struct PP {
  double h0, L, ramp, band, cd, Rm;
  int nb0, Lm;
};
__host__ __device__ static unsigned spr(unsigned v) {
  v &= 0xffffu;
  v = (v | (v << 8)) & 0x00ff00ffu;
  v = (v | (v << 4)) & 0x0f0f0f0fu;
  v = (v | (v << 2)) & 0x33333333u;
  v = (v | (v << 1)) & 0x55555555u;
  return v;
}
__host__ __device__ static unsigned cmp1(unsigned v) {
  v &= 0x55555555u;
  v = (v | (v >> 1)) & 0x33333333u;
  v = (v | (v >> 2)) & 0x0f0f0f0fu;
  v = (v | (v >> 4)) & 0x00ff00ffu;
  v = (v | (v >> 8)) & 0x0000ffffu;
  return v;
}
/* Morton code, x bit above y bit: the order of the FMM children q = 2 xbit + ybit */
__host__ __device__ static long mort(long x, long y) { return ((long)spr((unsigned)x) << 1) | spr((unsigned)y); }
__host__ __device__ static long fdivl(long a, long b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }
/* leaf of level l at level-l block (X, Y) (periodic), -1 if that block is not a leaf */
__host__ __device__ static int lmapd(int l, long X, long Y, int nb0, int Lm, const int *fine, const int *lev) {
  if (l < 0 || l > Lm)
    return -1;
  long n = (long)nb0 << l, nf = (long)nb0 << Lm;
  X = ((X % n) + n) % n;
  Y = ((Y % n) + n) % n;
  int c = fine[(Y << (Lm - l)) * nf + (X << (Lm - l))];
  return lev[c] == l ? c : -1;
}
__host__ __device__ static int adjq(int Dm, const int *bd, const int *px, const int *py, int a, int b, int sc) {
  long F = 1L << Dm, sa = Dm - bd[a], sb = Dm - bd[b];
  long ax0 = (long)px[a] << sa, ax1 = ((long)px[a] + 1) << sa;
  long ay0 = (long)py[a] << sa, ay1 = ((long)py[a] + 1) << sa;
  long bx0 = ((long)px[b] << sb) + (sc / 3 - 1) * F, bx1 = (((long)px[b] + 1) << sb) + (sc / 3 - 1) * F;
  long by0 = ((long)py[b] << sb) + (sc % 3 - 1) * F, by1 = (((long)py[b] + 1) << sb) + (sc % 3 - 1) * F;
  return ax0 <= bx1 && bx0 <= ax1 && ay0 <= by1 && by0 <= ay1;
}
/* smoothstep: tail of the Bernstein sum, sum_{j > d/2} C(d,j) t^j (1-t)^(d-j) */
__host__ __device__ static double bern(double t, int d) {
  if (t <= 0)
    return 0;
  if (t >= 1)
    return 1;
  double s = 0, c = 1, u = XS(1.0, t);
  for (int j = 0; j <= d; j++) {
    if (2 * j > d) {
      double a = c;
      for (int k = 0; k < j; k++)
        a = XM(a, t);
      for (int k = 0; k < d - j; k++)
        a = XM(a, u);
      s = XA(s, a);
    }
    c = XD(XM(c, (double)(d - j)), (double)(j + 1));
  }
  return s;
}
__host__ __device__ static double gap_(double a, double a0, double sz, double L) {
  double d = XS(a, XA(a0, XD(sz, 2.0)));
  d = XS(d, XM(L, rint(XD(d, L))));
  return fmax(0.0, XS(fabs(d), XD(sz, 2.0)));
}
/* coarser block c is a PoU candidate of block b */
__host__ __device__ static int pcand(int lb, int bxb, int byb, int lc, int bxc, int byc, PP q) {
  double hb = XD(q.h0, (double)(1 << lb)), bs = XM((double)B, hb), hc = XD(q.h0, (double)(1 << lc)),
         cs = XM((double)B, hc);
  double ox = XS(XA(XM((double)bxc, cs), XD(cs, 2.0)), XA(XM((double)bxb, bs), XD(bs, 2.0)));
  double oy = XS(XA(XM((double)byc, cs), XD(cs, 2.0)), XA(XM((double)byb, bs), XD(bs, 2.0)));
  ox = XS(ox, XM(q.L, rint(XD(ox, q.L))));
  oy = XS(oy, XM(q.L, rint(XD(oy, q.L))));
  double lim = XD(XM(q.Rm, hc), q.h0), hw = XD(XA(cs, bs), 2.0);
  return XS(fabs(ox), hw) < lim && XS(fabs(oy), hw) < lim;
}
/* PoU weight P and core width delta of one node (BORDER.md section 5) */
__host__ __device__ static void pnode(int l, double x, double y, long gi, long gj, const int *cand, int nc,
                                      const int *lev, const int *bx, const int *by, PP q, double *Po,
                                      double *Do) {
  double ph[DM + 2], sk[DM + 2], P = 0, de = q.h0;
  for (int k = 0; k <= l + 1; k++)
    ph[k] = sk[k] = k <= l;
  for (int e = 0; e < nc; e++) {
    int c = cand[e], lc = lev[c];
    double hc = XD(q.h0, (double)(1 << lc)), cs = XM((double)B, hc);
    double gx = gap_(x, XM((double)bx[c], cs), cs, q.L), gy = gap_(y, XM((double)by[c], cs), cs, q.L);
    double lim = XD(XM(q.Rm, hc), q.h0);
    if (gx >= lim || gy >= lim)
      continue;
    for (int k = lc + 1; k <= l; k++) {
      double hk = XD(q.h0, (double)(1 << k)), wr = XM(XM(q.ramp, 2.0), hk), o = XA(wr, XM(8.0, hk)),
             bw = XM(q.band, hk);
      ph[k] = XM(ph[k], XS(1.0, XM(XS(1.0, bern(XD(gx, wr), 23)), XS(1.0, bern(XD(gy, wr), 23)))));
      sk[k] = XM(sk[k], XS(1.0, XM(XS(1.0, bern(XD(XS(gx, o), bw), 15)), XS(1.0, bern(XD(XS(gy, o), bw), 15)))));
    }
  }
  for (int k = 0; k <= l; k++) {
    long s = 1L << (l - k);
    double hk = XD(q.h0, (double)(1 << k));
    if (gi % s == 0 && gj % s == 0)
      P = XA(P, XM(XM(XM(hk, hk), ph[k]), XS(1.0, ph[k + 1])));
    if (k > 0)
      de = XA(de, XM(XS(hk, XM(2.0, hk)), sk[k]));
  }
  *Po = P;
  *Do = XM(q.cd, de);
}
/* leaf e (image sc) is within R of leaf c: C-list candidate */
__host__ __device__ static int cnear(int lc, int bxc, int byc, int le, int bxe, int bye, int sc, double R,
                                     double h0, double L) {
  double hc = XD(h0, (double)(1 << lc)), s0 = XM((double)B, hc), he = XD(h0, (double)(1 << le)),
         se = XM((double)B, he);
  double ox = XS(XA(XA(XM((double)bxe, se), XD(se, 2.0)), XM((double)(sc / 3 - 1), L)), XA(XM((double)bxc, s0), XD(s0, 2.0)));
  double oy = XS(XA(XA(XM((double)bye, se), XD(se, 2.0)), XM((double)(sc % 3 - 1), L)), XA(XM((double)byc, s0), XD(s0, 2.0)));
  double hw = XD(XA(se, s0), 2.0), gx = fmax(0.0, XS(fabs(ox), hw)), gy = fmax(0.0, XS(fabs(oy), hw));
  return sqrt(XA(XM(gx, gx), XM(gy, gy))) < R;
}
/* value at node (X, Y) of an old block refined k times: recursive DD8
   prediction from the old lab (own nodes + ghost rows) */
__host__ __device__ static double predv(const double *lb, int W, int G, int ob, int k, long X, long Y, int *err) {
  const double DW[8] = {-5. / 2048, 49. / 2048, -245. / 2048, 1225. / 2048,
                        1225. / 2048, -245. / 2048, 49. / 2048, -5. / 2048};
  if (k == 0) {
    if (X < -G || X >= B + G || Y < -G || Y >= B + G) {
      *err = 2;
      return 0;
    }
    return lb[(long)ob * W * W + (Y + G) * W + X + G];
  }
  int ox = (int)(X & 1), oy = (int)(Y & 1);
  double v = 0;
  for (int a = 0; a < (ox ? 8 : 1); a++)
    for (int e = 0; e < (oy ? 8 : 1); e++) {
      long I = ox ? (X - 1) / 2 - 3 + a : X / 2, J = oy ? (Y - 1) / 2 - 3 + e : Y / 2;
      if (X < 0)
        I = ox ? -((-(X - 1) + 1) / 2) - 3 + a : -((-X + 1) / 2);
      if (Y < 0)
        J = oy ? -((-(Y - 1) + 1) / 2) - 3 + e : -((-Y + 1) / 2);
      v = XA(v, XM(XM(ox ? DW[a] : 1.0, oy ? DW[e] : 1.0), predv(lb, W, G, ob, k - 1, I, J, err)));
    }
  return v;
}
struct H {
  int nb0, Lm, D0, Dm, nleaf, nbox, *lev, *bx, *by, **map, **bmap, *bd, *px, *py, *ch, *par, *leaf, *lbox,
      lo[DM + 3], hi[DM + 3];
  int *uo, *ul, *us, *vo, *vl, *vs, *wo, *wl, *ws, *xo, *xl, *xs, *ro, *rl, *rs, *co, *cl, *cs;
  int *rel, *nbid;
};
static int findb(H *m, int d, int x, int y, int *sc) {
  int n = 1 << d, xx = ((x % n) + n) % n, yy = ((y % n) + n) % n;
  *sc = ((x - xx) / n + 1) * 3 + (y - yy) / n + 1;
  return m->bmap[d][(long)yy * n + xx];
}
static int adjb(H *m, int a, int b, int sc) { return adjq(m->Dm, m->bd, m->px, m->py, a, b, sc); }
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
static int *bandtg(int nb0, int Lm, double w) {
  long nf = (long)nb0 << Lm;
  int *tgt = (int *)hmalloc(nf * nf * sizeof(int));
  for (long j = 0; j < nf; j++)
    for (long i = 0; i < nf; i++) {
      double y0 = 2 * pi * j / nf, y1 = 2 * pi * (j + 1) / nf, dm = 1e300;
      for (int q = 0; q < 2; q++) {
        double c = q ? 3 * pi / 2 : pi / 2;
        dm = fmin(dm, c < y0 ? y0 - c : (c > y1 ? c - y1 : 0));
      }
      tgt[j * nf + i] = dm < w ? Lm : 0;
    }
  return tgt;
}
/* Host reference mesh build. Fine-cell targets tin -> leaf levels lv by the
   fixpoint of (a) the buffer rule: a level-l leaf keeps every cell within
   nbuf level-l block widths at level >= l - 1, (b) 2:1 balance. Within one
   pass the leaves are scanned in row-major order of their corner and a later
   write replaces an earlier one (the GPU reproduces this order with a 64-bit
   atomicMax on (scan position, value)). Leaves are numbered in Morton order
   of their corner on the fine lattice. */
static void build(H *m, int nb0, int Lm, const int *tin) {
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
  memcpy(tgt, tin, nf * nf * sizeof(int));
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
  for (long q = 0; q < nf * nf; q++) {
    long i = cmp1((unsigned)(q >> 1)), j = cmp1((unsigned)q);
    int l = lv[j * nf + i];
    long s = 1L << (Lm - l), n = (long)nb0 << l;
    if (i % s || j % s)
      continue;
    if (m->nleaf >= cap) {
      cap = cap ? 2 * cap : 1024;
      m->lev = (int *)realloc(m->lev, cap * sizeof(int));
      m->bx = (int *)realloc(m->bx, cap * sizeof(int));
      m->by = (int *)realloc(m->by, cap * sizeof(int));
    }
    m->lev[m->nleaf] = l;
    m->bx[m->nleaf] = i / s;
    m->by[m->nleaf] = j / s;
    m->map[l][(j / s) * n + i / s] = m->nleaf++;
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
/* soft: case enumeration, an impossible neighbour case only sets cerr */
static int soft, cerr;
static void cfail(const char *msg) {
  if (soft) {
    cerr = 1;
    return;
  }
  fprintf(stderr, "vm2: error: %s\n", msg);
  exit(1);
}
static void svcp(SV *d, const SV *s) {
  d->n = s->n;
  if (s->n > 0) {
    memcpy(d->i, s->i, s->n * sizeof(int));
    memcpy(d->w, s->w, s->n * sizeof(double));
  }
}
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
        svcp(&bufs[bd][dst + e], &bufs[bs][src + e]);
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
            cfail("sparse row too long");
            return;
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
    cfail("restriction stencil leaves the neighbourhood");
    return b * NB;
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
  cfail("restriction stencil meets a coarser block");
  return b * NB;
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
    cfail("restriction row too long");
    return;
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
        cfail("restriction row longer than SVM");
        return;
      }
      s->n = rn_;
      for (int e = 0; e < rn_; e++) {
        s->i[e] = ri_[e];
        s->w[e] = rw_[e];
      }
    }
}
/* ghost rows of block b of mesh m into bufs[0] (the lab of sparse rows) */
static void compile1(H *m, Tab *T, int b, int smooth, SV **bufs) {
  int W = T->hd[4], G = T->hd[2], C = T->hd[6], WW = W * W;
  SV *lab = bufs[0], *buf = bufs[1], *tmp = bufs[2];
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
    for (int i = -G; i < B + G; i++)
      if ((i < 0 || i >= B || j < 0 || j >= B) && lab[(j + G) * W + i + G].n < 0) {
        if (!soft)
          fprintf(stderr, "vm2: error: ghost (%d %d) of block %d not filled\n", i, j, b);
        cfail("ghost not filled");
        return;
      }
}
static int G, W, NG;
static int *hgpos; /* lab position of ghost row g of a block, g < NG */
static void compile(H *m, Tab *T, CSR *g, int smooth) {
  int W = T->hd[4], C = T->hd[6], WW = W * W, NG = WW - NB;
  SV *bufs[3] = {(SV *)hmalloc(WW * sizeof(SV)), (SV *)hmalloc(C * C * sizeof(SV)),
                 (SV *)hmalloc(W * C * sizeof(SV))};
  long cn = 1 << 22;
  g->nr = (long)m->nleaf * NG;
  g->row = (long *)hmalloc(g->nr * sizeof(long));
  g->p = (int *)hmalloc((g->nr + 1) * sizeof(int));
  g->c = (int *)hmalloc(cn * sizeof(int));
  g->w = (double *)hmalloc(cn * sizeof(double));
  g->nnz = 0;
  for (int b = 0; b < m->nleaf; b++) {
    compile1(m, T, b, smooth, bufs);
    for (int q = 0; q < NG; q++) {
      SV *s = &bufs[0][hgpos[q]];
      long r = (long)b * NG + q;
      while (g->nnz + s->n >= cn) {
        cn *= 2;
        g->c = (int *)realloc(g->c, cn * sizeof(int));
        g->w = (double *)realloc(g->w, cn * sizeof(double));
      }
      g->row[r] = (long)b * WW + hgpos[q];
      g->p[r] = g->nnz;
      for (int e = 0; e < s->n; e++) {
        g->c[g->nnz] = s->i[e];
        g->w[g->nnz++] = s->w[e];
      }
    }
  }
  g->p[g->nr] = g->nnz;
  for (int q = 0; q < 3; q++)
    free(bufs[q]);
}
/* Neighbour cases. A block's ghost rows depend only on the relation of its 8
   neighbours (same, coarse, fine) and, with a coarse neighbour, its parity.
   Rows are compiled once per case with symbolic columns slot * NB + node:
   slot k * 4 + q = quadrant q of neighbour k (q = 0 for same level), 16 = own,
   36 + g = coarse block at offset g = (ox + 1) * 3 + oy + 1 from the parent
   (two neighbour slots can share one coarse block). On the GPU a block reads
   its case rows and maps slots to blocks (slot table, NS per block). */
static int ncase, *hcidx, *hcrow, *hccol;
static double *hcw;
static long cnnz;
static int *dcidx, *dcrow, *dccol, *dgpos;
static double *dcw;
__host__ __device__ static int grp(int k, int par) {
  int cx = k / 3 - 1, cy = k % 3 - 1, a = (par & 1) + cx, b = (par >> 1) + cy;
  int ox = a >= 0 ? a / 2 : -1, oy = b >= 0 ? b / 2 : -1;
  return (ox + 1) * 3 + oy + 1;
}
static void cases(Tab *T, int smooth) {
  int WW = W * W, C = T->hd[6], NG = WW - NB;
  int rel[17 * 9], nbid[17 * 36], bxv[17], byv[17];
  H v;
  memset(&v, 0, sizeof v);
  memset(bxv, 0, sizeof bxv);
  memset(byv, 0, sizeof byv);
  v.rel = rel;
  v.nbid = nbid;
  v.bx = bxv;
  v.by = byv;
  SV *bufs[3] = {(SV *)hmalloc(WW * sizeof(SV)), (SV *)hmalloc(C * C * sizeof(SV)),
                 (SV *)hmalloc(W * C * sizeof(SV))};
  long cap = 1 << 22, rcap = 1024;
  hcidx = (int *)hmalloc(NCODE * sizeof(int));
  hcrow = (int *)hmalloc(rcap * (NG + 1) * sizeof(int));
  hccol = (int *)hmalloc(cap * sizeof(int));
  hcw = (double *)hmalloc(cap * sizeof(double));
  for (int k = 0; k < NCODE; k++)
    hcidx[k] = -1;
  ncase = 0;
  cnnz = 0;
  for (int par = 0; par < 4; par++)
    for (int code = 0; code < 6561; code++) {
      int r[9], anyc = 0, ok = 1;
      for (int k = 0, c = code; k < 9; k++) {
        if (k == 4) {
          r[k] = 0;
          continue;
        }
        r[k] = c % 3 + 1;
        c /= 3;
        anyc |= r[k] == 2;
      }
      if (!anyc && par)
        continue;
      for (int k = 0; k < 9; k++) {
        if (r[k] != 2)
          continue;
        if (grp(k, par) == 4)
          ok = 0;
        for (int k2 = 0; k2 < 9; k2++) {
          if (k2 != 4 && grp(k2, par) == grp(k, par) && r[k2] != 2)
            ok = 0;
          if (r[k2] == 3 && abs(k / 3 - k2 / 3) <= 1 && abs(k % 3 - k2 % 3) <= 1)
            ok = 0;
        }
      }
      if (!ok)
        continue;
      for (int k = 0; k < 9; k++) {
        rel[16 * 9 + k] = r[k];
        for (int q = 0; q < 4; q++)
          nbid[(16 * 9 + k) * 4 + q] = r[k] == 3 ? k * 4 + q : -1;
        if (r[k] == 1)
          nbid[(16 * 9 + k) * 4] = k * 4;
        if (r[k] == 2)
          nbid[(16 * 9 + k) * 4] = 36 + grp(k, par);
      }
      bxv[16] = par & 1;
      byv[16] = par >> 1;
      soft = 1;
      cerr = 0;
      compile1(&v, T, 16, smooth, bufs);
      soft = 0;
      if (cerr)
        continue;
      hcidx[par * 6561 + code] = ncase;
      if (ncase + 1 >= rcap) {
        rcap *= 2;
        hcrow = (int *)realloc(hcrow, rcap * (NG + 1) * sizeof(int));
      }
      for (int g = 0; g < NG; g++) {
        SV *s = &bufs[0][hgpos[g]];
        hcrow[(long)ncase * (NG + 1) + g] = (int)cnnz;
        while (cnnz + s->n >= cap) {
          cap *= 2;
          hccol = (int *)realloc(hccol, cap * sizeof(int));
          hcw = (double *)realloc(hcw, cap * sizeof(double));
        }
        for (int e = 0; e < s->n; e++) {
          hccol[cnnz] = s->i[e];
          hcw[cnnz++] = s->w[e];
        }
      }
      hcrow[(long)ncase * (NG + 1) + NG] = (int)cnnz;
      ncase++;
    }
  for (int q = 0; q < 3; q++)
    free(bufs[q]);
  dcidx = (int *)dmalloc(NCODE * sizeof(int));
  up(dcidx, hcidx, NCODE * sizeof(int));
  dcrow = (int *)dmalloc((long)ncase * (NG + 1) * sizeof(int));
  up(dcrow, hcrow, (long)ncase * (NG + 1) * sizeof(int));
  dccol = (int *)dmalloc((cnnz + 1) * sizeof(int));
  up(dccol, hccol, cnnz * sizeof(int));
  dcw = (double *)dmalloc((cnnz + 1) * sizeof(double));
  up(dcw, hcw, cnnz * sizeof(double));
  dgpos = (int *)dmalloc(NG * sizeof(int));
  up(dgpos, hgpos, NG * sizeof(int));
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
static double2 *gTm2m, *gTl2l, *gTm2l, *gTfar;
static int gPt;
static void tables(const double2 *F, int Pt) {
  size_t q = (size_t)Pt * Pt * sizeof(double2);
  gPt = Pt;
  gTm2m = (double2 *)dmalloc(4 * q);
  gTl2l = (double2 *)dmalloc(4 * q);
  gTm2l = (double2 *)dmalloc(49 * q);
  gTfar = (double2 *)dmalloc(q);
  up(gTm2m, F, 4 * q);
  up(gTl2l, F + 4L * Pt * Pt, 4 * q);
  up(gTm2l, F + 8L * Pt * Pt, 49 * q);
  up(gTfar, F + 57L * Pt * Pt, q);
}
/* particle-sized FMM arrays; tree and list arrays are set by setup (host
   mesh) or dmesh (GPU mesh) */
static void talloc(T *t, long N, int nleaf, int nbox, int p) {
  int P = p + 1;
  t->Pt = gPt;
  t->Tm2m = gTm2m;
  t->Tl2l = gTl2l;
  t->Tm2l = gTm2l;
  t->Tfar = gTfar;
  t->M = (double2 *)dmalloc((size_t)nbox * P * sizeof(double2));
  t->Lc = (double2 *)dmalloc((size_t)nbox * P * sizeof(double2));
  t->key = (int *)dmalloc(N * sizeof(int));
  t->start = (int *)dmalloc((nleaf + 1) * sizeof(int));
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
  t->of = (int *)dmalloc(N * sizeof(int));
  ck(cudaMemset(t->of, 0, N * sizeof(int)), "cudaMemset");
}
static int *hfine(H *m) {
  int Fn = 1 << m->Dm, *fine = (int *)hmalloc((size_t)Fn * Fn * sizeof(int));
  for (int c = 0; c < m->nleaf; c++) {
    int b = m->lbox[c], s = 1 << (m->Dm - m->bd[b]);
    for (int j = 0; j < s; j++)
      for (int i = 0; i < s; i++)
        fine[(long)(m->py[b] * s + j) * Fn + m->px[b] * s + i] = c;
  }
  return fine;
}
static void setup(H *m, T *t, long N, int p, double L) {
  t->N = N;
  t->p = p;
  t->L = L;
  t->nleaf = m->nleaf;
  t->nbox = m->nbox;
  t->Dm = m->Dm;
  t->Fn = 1 << m->Dm;
  int *fine = hfine(m);
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
  t->co = idup(m->co, m->nleaf + 1);
  t->cl = idup(m->cl, m->co[m->nleaf] + 1);
  t->cs = idup(m->cs, m->co[m->nleaf] + 1);
  talloc(t, N, m->nleaf, m->nbox, p);
}
static long nfpair;
static int graphing, ptime;
static double ptacc[8];
static cudaEvent_t pev[8];
static unsigned long long *dpair;
static void pmark(int k) {
  if (ptime && !graphing)
    ck(cudaEventRecord(pev[k], cudaStreamPerThread), "cudaEventRecord");
}
/* pair count (statistics) on the GPU: no copy of the leaf starts to the host */
__global__ static void kpair(T t, unsigned long long *acc) {
  int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= t.nleaf)
    return;
  long s2 = 0;
  for (int e = t.uo[c]; e < t.uo[c + 1]; e++)
    s2 += t.start[t.ul[e] + 1] - t.start[t.ul[e]];
  for (int e = t.xo[c]; e < t.xo[c + 1]; e++)
    s2 += t.start[t.xl[e] + 1] - t.start[t.xl[e]];
  for (int e = t.co[c]; e < t.co[c + 1]; e++)
    s2 += t.start[t.cl[e] + 1] - t.start[t.cl[e]];
  atomicAdd(acc, (unsigned long long)(s2 * (t.start[c + 1] - t.start[c])));
}
static int nozbar;
static void fmm(T *t, H *m, int mode, const double *x, const double *y,
                const double *g, double *ku, double *kv, double *kg) {
  long N = t->N;
  int P = t->p + 1;
  static int *k2 = NULL;
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
    }
    k2 = (int *)dmalloc(N * sizeof(int));
    perm0 = (long *)dmalloc(N * sizeof(long));
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
  kpair<<<grid(t->nleaf), nthread>>>(*t, dpair);
  launched("kpair");
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
static CSR gh;
static T t;
static long N;
static int n0, Lm, nosm, p = 24, mode; /* mode 0: GPU regrid, 1: host reference (-H), 2: both, compared (-X) */
static double L, h0, wband, cd, nu, ramp = 16, band = 16, *dsm2;
static int *dlev, *dbx, *dby, *dgp, *dgc, *hf, *drel, *dslot, *dcas, *dbmap, *dbleaf, *derr;
static long *dgr;
static double *dgw, *lab, *lab2, *lab3, *hh, *hx, *hy, *hg, *hu, *hv, *hgr, rgsum, *dmxd;
static St s1, s0, wk[4], S[3], u, dd;
struct GM {
  int nleaf, nbox, nu, nw, nv, nx, nc, lo[DM + 3];
};
static GM g;   /* sizes of the GPU-built mesh */
static Off bo; /* depth d boxes: entries [bo.o[d], bo.o[d] + 4^d) of the box flag/map arrays */
static PP pp;
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
               t.odg, t.w, t.ow, t.oh, t.ov, t.of, dlev, dbx, dby, dgp, dgc, dgr, dgw, lab, lab2,
               lab3, t.o, t.oo, t.co, t.cl, t.cs, drel, dslot, dcas, dbmap, dbleaf, dmxd};
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
  void *h[] = {hh, hx, hy, hg, hu, hv, hgr, hf};
  for (size_t k = 0; k < sizeof h / sizeof *h; k++)
    free(h[k]);
  memset(&t, 0, sizeof t);
  dlev = dbx = dby = dgp = dgc = drel = dslot = dcas = dbmap = dbleaf = hf = NULL;
  dgr = NULL;
  dgw = lab = lab2 = lab3 = hh = hx = hy = hg = hu = hv = hgr = dmxd = NULL;
}
static void hfree(H *a, CSR *g) {
  if (a->map)
    for (int l = 0; l <= a->Lm; l++)
      free(a->map[l]);
  if (a->bmap)
    for (int d = 0; d <= a->Dm; d++)
      free(a->bmap[d]);
  void *h[] = {a->map, a->bmap, a->lev, a->bx, a->by, a->bd, a->px, a->py, a->ch, a->par, a->leaf,
               a->lbox, a->uo, a->ul, a->us, a->vo, a->vl, a->vs, a->wo, a->wl, a->ws, a->xo,
               a->xl, a->xs, a->ro, a->rl, a->rs, a->co, a->cl, a->cs, a->rel, a->nbid, g->row, g->p, g->c, g->w};
  for (size_t k = 0; k < sizeof h / sizeof *h; k++)
    free(h[k]);
}
/* host reference: PoU weights, core widths, C lists (BORDER.md section 5) */
static void pou(void) {
  int nl = m.nleaf, *cand = (int *)hmalloc(nl * sizeof(int));
  double *dmx = (double *)hmalloc(nl * sizeof(double));
  for (int b = 0; b < nl; b++) {
    int l = m.lev[b], nc = 0;
    for (int c = 0; c < nl; c++)
      if (m.lev[c] < l && pcand(l, m.bx[b], m.by[b], m.lev[c], m.bx[c], m.by[c], pp))
        cand[nc++] = c;
    dmx[b] = 0;
    for (int q = 0; q < NB; q++) {
      long n = (long)b * NB + q, gi = (long)m.bx[b] * B + q % B, gj = (long)m.by[b] * B + q / B;
      pnode(l, hx[n], hy[n], gi, gj, cand, nc, m.lev, m.bx, m.by, pp, &hh[n], &hgr[n]);
      dmx[b] = fmax(dmx[b], hgr[n]);
    }
  }
  int nc = 0, cc = 0, *cl = NULL, *cs = NULL;
  m.co = (int *)hmalloc((nl + 1) * sizeof(int));
  for (int c = 0; c < nl; c++) {
    m.co[c] = nc;
    double R = XM(7.0, dmx[c]);
    for (int e = 0; e < nl; e++)
      for (int sc = 0; sc < 9; sc++) {
        if (!cnear(m.lev[c], m.bx[c], m.by[c], m.lev[e], m.bx[e], m.by[e], sc, R, h0, L))
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
  m.co[nl] = nc;
  m.cl = cl ? cl : (int *)hmalloc(sizeof(int));
  m.cs = cs ? cs : (int *)hmalloc(sizeof(int));
  free(cand);
  free(dmx);
}
/* host reference mesh: build, ghost rows, node positions, weights */
static void hmesh(const int *tgt) {
  build(&m, n0 / B, Lm, tgt);
  compile(&m, &tab, &gh, !nosm);
  N = (long)m.nleaf * NB;
  hh = (double *)hmalloc(N * sizeof(double));
  hx = (double *)hmalloc(N * sizeof(double));
  hy = (double *)hmalloc(N * sizeof(double));
  hgr = (double *)hmalloc(N * sizeof(double));
  for (long n = 0; n < N; n++) {
    int b = n / NB, i = n % NB % B, j = n % NB / B;
    double hb = h0 / (1 << m.lev[b]);
    hx[n] = (m.bx[b] * B + i) * hb;
    hy[n] = (m.by[b] * B + j) * hb;
  }
  pou();
}
/* arrays common to both paths, allocated after the mesh */
static void common(void) {
  size_t nl = N / NB;
  lab = (double *)dmalloc(nl * W * W * sizeof(double));
  lab2 = (double *)dmalloc(nl * W * W * sizeof(double));
  lab3 = (double *)dmalloc(nl * W * W * sizeof(double));
  s1 = salloc(N);
  s0 = salloc(N);
  u = salloc(N);
  dd = salloc(N);
  for (int q = 0; q < 4; q++)
    wk[q] = salloc(N);
  for (int q = 0; q < 3; q++)
    S[q] = salloc(N);
  hg = (double *)hmalloc(N * sizeof(double));
  hu = (double *)hmalloc(N * sizeof(double));
  hv = (double *)hmalloc(N * sizeof(double));
  hf = (int *)hmalloc(N * sizeof(int));
  nodes<<<grid(N), nthread>>>(N, h0, dlev, dbx, dby, s0.x, s0.y, t.h);
  launched("nodes");
}
/* host path: upload the host-built mesh */
static void hupload(void) {
  memset(&t, 0, sizeof t);
  setup(&m, &t, N, p, L);
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
  t.ov = dup(hh, N * sizeof(double));
  t.oh = dup(hgr, N * sizeof(double));
  common();
}
/* ---------------- GPU regrid ---------------- */
__global__ static void bmax(const double *d, double *o) {
  __shared__ double s[NB];
  int b = blockIdx.x, k = threadIdx.x;
  s[k] = d[(long)b * NB + k];
  __syncthreads();
  for (int h = NB / 2; h > 0; h /= 2) {
    if (k < h)
      s[k] = fmax(s[k], s[k + h]);
    __syncthreads();
  }
  if (k == 0)
    o[b] = fmax(0.0, s[0]);
}
__global__ static void amax(long N, const double *x, double *o) {
  __shared__ double a[1024];
  double s = 0;
  for (long i = threadIdx.x; i < N; i += blockDim.x)
    s = fmax(s, fabs(x[i]));
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
/* ghost rows from the case tables: row g of block b = case rows of cas[b]
   with symbolic slots mapped to blocks by the slot table */
__global__ static void labcase(long nr, int NG, const int *cas, const int *slot, const int *crow,
                               const int *ccol, const double *cw, const int *gpos, int WW, const double *v,
                               double *lb) {
  long r = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (r >= nr)
    return;
  int b = r / NG, q = r % NG;
  const int *rp = crow + (long)cas[b] * (NG + 1);
  double s = 0;
  for (int e = rp[q]; e < rp[q + 1]; e++) {
    int y = ccol[e];
    s = fma(cw[e], v[(long)slot[b * NS + y / NB] * NB + y % NB], s);
  }
  lb[(long)b * WW + gpos[q]] = s;
}
__global__ static void kflag(int nl, int Lm, double er, double kap, double tau, const int *lev, const double *db,
                             const double *sb, const double *wm, int *nlv, int *fl) {
  int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= nl)
    return;
  double d = db[b] / wm[0] * exp(kap * sb[b] * tau), ec = er / 256;
  int l = lev[b];
  fl[b] = d > er;
  nlv[b] = d > er && l < Lm ? l + 1 : (d < ec && l > 0 ? l - 1 : l);
}
/* coarsen only if all four siblings are leaves that want to coarsen */
__global__ static void kcoarse(int nl, int nb0, int Lm, const int *lev, const int *bx, const int *by,
                               const int *fine, const int *nlv, int *o) {
  int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= nl)
    return;
  int l = lev[b], v = nlv[b];
  if (v < l) {
    int px = bx[b] & ~1, py = by[b] & ~1, ok = 1;
    for (int q = 0; q < 4; q++) {
      int c = lmapd(l, px + (q & 1), py + (q >> 1), nb0, Lm, fine, lev);
      ok &= c >= 0 && nlv[c] < l;
    }
    if (!ok)
      v = l;
  }
  o[b] = v;
}
__global__ static void ktg0(long n, const int *fine, const int *nlv, int *tg) {
  long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (q < n)
    tg[q] = nlv[fine[q]];
}
/* a flagged block asks for its new level over nbuf new-level block widths */
__global__ static void ktg1(int nl, int nb0, int Lm, int nbuf, const int *lev, const int *bx, const int *by,
                            const int *fl, const int *nlv, int *tg) {
  int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= nl || !fl[b])
    return;
  long nf = (long)nb0 << Lm, s = 1L << (Lm - lev[b]), r = (1L << (Lm - nlv[b])) * nbuf;
  for (long j = (long)by[b] * s - r; j < (long)(by[b] + 1) * s + r; j++)
    for (long i = (long)bx[b] * s - r; i < (long)(bx[b] + 1) * s + r; i++)
      atomicMax(&tg[((j % nf + nf) % nf) * nf + (i % nf + nf) % nf], nlv[b]);
}
__global__ static void kpyr(long n, const int *src, int *dst) {
  long k = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (k >= n * n)
    return;
  long X = k % n, Y = k / n, n2 = 2 * n;
  int a = max(src[2 * Y * n2 + 2 * X], src[2 * Y * n2 + 2 * X + 1]),
      c = max(src[(2 * Y + 1) * n2 + 2 * X], src[(2 * Y + 1) * n2 + 2 * X + 1]);
  dst[k] = max(a, c);
}
/* leaf level of a fine cell: the coarsest l whose block has max target <= l */
__global__ static void klv(long nf, int nb0, int Lm, Off po, const int *pyr, int *lv) {
  long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (q >= nf * nf)
    return;
  long i = q % nf, j = q / nf;
  int l = 0;
  for (; l < Lm; l++)
    if (pyr[po.o[l] + (j >> (Lm - l)) * ((long)nb0 << l) + (i >> (Lm - l))] <= l)
      break;
  lv[q] = l;
}
/* buffer rule and 2:1 balance, one thread per leaf (its corner cell). The
   host pass scans leaves in row-major order of the corner, buffer writes
   before the 2:1 write, and the last write wins. Each write here is the key
   (2 q + 1 + own) << 32 | level, q = corner index: the atomicMax keeps the
   last write of the host order. */
__global__ static void kchk(long nf, int Lm, int nbuf, const int *lv, unsigned long long *tw, int *chg) {
  long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (q >= nf * nf)
    return;
  long i = q % nf, j = q / nf;
  int l = lv[q], mx = 0, ch = 0;
  long s = 1L << (Lm - l);
  if (i % s || j % s)
    return;
  long r = l > 1 ? s * nbuf : s;
  unsigned long long kb = (unsigned long long)(2 * q + 1) << 32, ko = (unsigned long long)(2 * q + 2) << 32;
  for (long jj = j - r; jj < j + s + r; jj++)
    for (long ii = i - r; ii < i + s + r; ii++) {
      long qq = ((jj + nf) % nf) * nf + (ii + nf) % nf;
      if (jj >= j - s && jj < j + 2 * s && ii >= i - s && ii < i + 2 * s)
        mx = max(mx, lv[qq]);
      if (lv[qq] < l - 1) {
        atomicMax(&tw[qq], kb | (unsigned)(l - 1));
        ch = 1;
      }
    }
  if (mx > l + 1) {
    for (long jj = j; jj < j + s; jj++)
      for (long ii = i; ii < i + s; ii++)
        atomicMax(&tw[jj * nf + ii], ko | (unsigned)(mx - 1));
    ch = 1;
  }
  if (ch)
    *chg = 1;
}
__global__ static void kapply(long n, const unsigned long long *tw, int *tg) {
  long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (q < n && tw[q])
    tg[q] = (int)(tw[q] & 0xffffffffu);
}
/* leaf flags in Morton order of the fine lattice */
__global__ static void kleaf(long nf, int Lm, const int *lv, int *fl) {
  long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (q >= nf * nf)
    return;
  long i = cmp1((unsigned)(q >> 1)), j = cmp1((unsigned)q), s = 1L << (Lm - lv[j * nf + i]);
  fl[q] = i % s == 0 && j % s == 0;
}
__global__ static void kleafw(long nf, int Lm, const int *lv, const int *fl, const int *ix, int *lev, int *bx,
                              int *by) {
  long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (q >= nf * nf || !fl[q])
    return;
  long i = cmp1((unsigned)(q >> 1)), j = cmp1((unsigned)q);
  int l = lv[j * nf + i], c = ix[q];
  long s = 1L << (Lm - l);
  lev[c] = l;
  bx[c] = i / s;
  by[c] = j / s;
}
__global__ static void kfine(long nf, int Lm, const int *lv, const int *ix, int *fine) {
  long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (q >= nf * nf)
    return;
  long i = q % nf, j = q / nf, s = 1L << (Lm - lv[q]);
  fine[q] = ix[mort(i / s * s, j / s * s)];
}
/* neighbour relations, slot table and neighbour case of every leaf */
__global__ static void krel(int nl, int nb0, int Lm, const int *lev, const int *bx, const int *by,
                            const int *fine, const int *cidx, int *rel, int *slot, int *cas, int *err) {
  int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= nl)
    return;
  int l = lev[c], r[9], sl[NS], par = (bx[c] & 1) + 2 * (by[c] & 1), anyc = 0;
  for (int k = 0; k < NS; k++)
    sl[k] = -1;
  sl[16] = c;
  for (int k = 0; k < 9; k++) {
    long X = bx[c] + k / 3 - 1, Y = by[c] + k % 3 - 1, n = (long)nb0 << l;
    int p;
    r[k] = 0;
    if (k == 4)
      continue;
    if ((p = lmapd(l, X, Y, nb0, Lm, fine, lev)) >= 0) {
      r[k] = 1;
      sl[k * 4] = p;
    } else if ((p = lmapd(l - 1, ((X % n) + n) % n >> 1, ((Y % n) + n) % n >> 1, nb0, Lm, fine, lev)) >= 0) {
      r[k] = 2;
      anyc = 1;
      sl[36 + grp(k, par)] = p;
    } else {
      r[k] = 3;
      for (int q = 0; q < 4; q++) {
        p = lmapd(l + 1, 2 * X + (q & 1), 2 * Y + (q >> 1), nb0, Lm, fine, lev);
        if (p < 0)
          atomicMax(err, 3);
        sl[k * 4 + q] = p;
      }
    }
  }
  int code = 0;
  for (int k = 8; k >= 0; k--)
    if (k != 4)
      code = code * 3 + r[k] - 1;
  code += (anyc ? par : 0) * 6561;
  for (int k = 0; k < 9; k++)
    rel[c * 9 + k] = r[k];
  for (int k = 0; k < NS; k++)
    slot[c * NS + k] = sl[k];
  cas[c] = cidx[code];
  if (cidx[code] < 0)
    atomicMax(err, 4);
}
/* boxes of the FMM tree: box (d, X, Y) exists iff the leaf covering its
   corner is at depth >= d; BFS numbering = rank in (depth, Morton) order */
__global__ static void kbox(long tot, Off bo, int Dm, int D0, long Fn, const int *fine, const int *lev, int *ex) {
  long k = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (k >= tot)
    return;
  int d = 0;
  while (k >= bo.o[d + 1])
    d++;
  long q = k - bo.o[d], X = cmp1((unsigned)(q >> 1)), Y = cmp1((unsigned)q);
  int c = fine[(Y << (Dm - d)) * Fn + (X << (Dm - d))];
  ex[k] = d <= D0 + lev[c];
}
__global__ static void klohi(int n, Off bo, const int *bix, int *lh) {
  int d = threadIdx.x;
  if (d < n)
    lh[d] = bix[bo.o[d]];
}
__global__ static void kboxw(long tot, Off bo, int Dm, int D0, long Fn, const int *fine, const int *lev,
                             const int *ex, const int *bix, int *bd, int *px, int *py, int *ch, int *par,
                             int *leaf, int *lbox, int *bmap) {
  long k = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (k >= tot)
    return;
  int d = 0;
  while (k >= bo.o[d + 1])
    d++;
  long q = k - bo.o[d], X = cmp1((unsigned)(q >> 1)), Y = cmp1((unsigned)q), n = 1L << d;
  long rm = bo.o[d] + Y * n + X;
  if (!ex[k]) {
    bmap[rm] = -1;
    return;
  }
  int b = bix[k], c = fine[(Y << (Dm - d)) * Fn + (X << (Dm - d))], dl = D0 + lev[c];
  bmap[rm] = b;
  bd[b] = d;
  px[b] = X;
  py[b] = Y;
  leaf[b] = d == dl ? c : -1;
  if (d == dl)
    lbox[c] = b;
  ch[b] = d < dl ? bix[bo.o[d + 1] + mort(2 * X, 2 * Y)] : -1;
  par[b] = d > 0 ? bix[bo.o[d - 1] + mort(X >> 1, Y >> 1)] : -1;
}
__global__ static void kcen(int nbox, double L, const int *bd, const int *px, const int *py, double *cx,
                            double *cy) {
  int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= nbox)
    return;
  double h = L / (double)(1L << bd[b]);
  cx[b] = (px[b] + 0.5) * h;
  cy[b] = (py[b] + 0.5) * h;
}
struct TV {
  int Dm, nleaf, nbox;
  const int *bd, *px, *py, *ch, *par, *leaf, *lbox, *bmap;
  Off bo;
};
__device__ static int findv(const TV &v, int d, int x, int y, int *sc) {
  int n = 1 << d, xx = ((x % n) + n) % n, yy = ((y % n) + n) % n;
  *sc = ((x - xx) / n + 1) * 3 + (y - yy) / n + 1;
  return v.bmap[v.bo.o[d] + (long)yy * n + xx];
}
/* U and W lists of a leaf (same traversal as the host build), count or fill */
__global__ static void kul(TV v, int fill, int *uc, int *wc, const int *uo, int *ul, int *us, const int *wo,
                           int *wl, int *ws, int *err) {
  int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= v.nleaf)
    return;
  int b = v.lbox[c], d = v.bd[b], nu = 0, nwn = 0, n1 = 0, ns = 0;
  int pl[UL1], ps[UL1], stk[SK], sts[SK];
  for (int dd = 0; dd < d; dd++) {
    int ax = v.px[b] >> (d - dd), ay = v.py[b] >> (d - dd);
    for (int dx = -1; dx <= 1; dx++)
      for (int dy = -1; dy <= 1; dy++) {
        int sc, e = findv(v, dd, ax + dx, ay + dy, &sc), dup = 0;
        if (e < 0 || v.ch[e] >= 0 || !adjq(v.Dm, v.bd, v.px, v.py, b, e, sc))
          continue;
        for (int k = 0; k < n1; k++)
          dup |= pl[k] == v.leaf[e] && ps[k] == sc;
        if (dup)
          continue;
        if (n1 >= UL1) {
          atomicMax(err, 5);
          return;
        }
        pl[n1] = v.leaf[e];
        ps[n1++] = sc;
        if (fill) {
          ul[uo[c] + nu] = v.leaf[e];
          us[uo[c] + nu] = sc;
        }
        nu++;
      }
  }
  for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++) {
      int sc, e = findv(v, d, v.px[b] + dx, v.py[b] + dy, &sc);
      if (e >= 0) {
        stk[ns] = e;
        sts[ns++] = sc;
      }
    }
  while (ns > 0) {
    int e = stk[--ns], sc = sts[ns];
    if (v.ch[e] < 0) {
      if (fill) {
        ul[uo[c] + nu] = v.leaf[e];
        us[uo[c] + nu] = sc;
      }
      nu++;
      continue;
    }
    for (int q = 0; q < 4; q++) {
      int k = v.ch[e] + q;
      if (adjq(v.Dm, v.bd, v.px, v.py, b, k, sc)) {
        if (ns >= SK) {
          atomicMax(err, 5);
          return;
        }
        stk[ns] = k;
        sts[ns++] = sc;
      } else {
        if (fill) {
          wl[wo[c] + nwn] = k;
          ws[wo[c] + nwn] = sc;
        }
        nwn++;
      }
    }
  }
  if (!fill) {
    uc[c] = nu;
    wc[c] = nwn;
  }
}
/* V list and X0 (X entries owned by this box) of every box */
__global__ static void kvx(TV v, int fill, int *vc, int *xc, const int *vo, int *vl, int *vs, const int *xo,
                           int *xl, int *xs, int *err) {
  int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= v.nbox)
    return;
  int d = v.bd[b], pb = v.par[b], nv = 0, nx = 0, pl[UL1], ps[UL1];
  if (d > 0) {
    for (int dx = -1; dx <= 1; dx++)
      for (int dy = -1; dy <= 1; dy++) {
        int sc, e = findv(v, d - 1, v.px[pb] + dx, v.py[pb] + dy, &sc);
        if (e < 0 || v.ch[e] < 0)
          continue;
        for (int q = 0; q < 4; q++)
          if (!adjq(v.Dm, v.bd, v.px, v.py, b, v.ch[e] + q, sc)) {
            if (fill) {
              vl[vo[b] + nv] = v.ch[e] + q;
              vs[vo[b] + nv] = sc;
            }
            nv++;
          }
      }
    for (int dd = 0; dd < d; dd++) {
      int ax = v.px[b] >> (d - dd), ay = v.py[b] >> (d - dd);
      for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++) {
          int sc, e = findv(v, dd, ax + dx, ay + dy, &sc), dup = 0;
          if (e < 0 || v.ch[e] >= 0)
            continue;
          if (!adjq(v.Dm, v.bd, v.px, v.py, pb, e, sc) || adjq(v.Dm, v.bd, v.px, v.py, b, e, sc))
            continue;
          for (int k = 0; k < nx; k++)
            dup |= pl[k] == v.leaf[e] && ps[k] == sc;
          if (dup)
            continue;
          if (nx >= UL1) {
            atomicMax(err, 5);
            return;
          }
          pl[nx] = v.leaf[e];
          ps[nx] = sc;
          if (fill) {
            xl[xo[b] + nx] = v.leaf[e];
            xs[xo[b] + nx] = sc;
          }
          nx++;
        }
    }
  }
  if (!fill) {
    vc[b] = nv;
    xc[b] = nx;
  }
}
/* X list of a leaf: X0 of its ancestors and itself, shallow first */
__global__ static void kxl(TV v, int fill, const int *x0o, const int *x0l, const int *x0s, int *xc, const int *xo,
                           int *xl, int *xs, int *own) {
  int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= v.nleaf)
    return;
  int anc[DM + 2], na = 0, n = 0;
  for (int a = v.lbox[c]; a >= 0; a = v.par[a])
    anc[na++] = a;
  for (int k = na - 1; k >= 0; k--)
    for (int e = x0o[anc[k]]; e < x0o[anc[k] + 1]; e++) {
      if (fill) {
        xl[xo[c] + n] = x0l[e];
        xs[xo[c] + n] = x0s[e];
        own[xo[c] + n] = c;
      }
      n++;
    }
  if (!fill)
    xc[c] = n;
}
__global__ static void iotai(long n, int *a) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i < n)
    a[i] = (int)i;
}
__global__ static void kr(long n, const int *sv, const int *own, const int *xs, int *rl, int *rs) {
  long k = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (k >= n)
    return;
  int e = sv[k];
  rl[k] = own[e];
  rs[k] = (2 - xs[e] / 3) * 3 + 2 - xs[e] % 3;
}
/* C list: leaves within 7 delta_max, found on the fine lattice, sorted (e, sc) */
__global__ static void kcl(T t, PP q, const int *lev, const int *bx, const int *by, const double *dmx, int fill,
                           int *cc, int *err) {
  int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= t.nleaf)
    return;
  int l = lev[c], key[CL], nk = 0, last = -1;
  long nf = (long)q.nb0 << q.Lm, s = 1L << (q.Lm - l);
  double R = XM(7.0, dmx[c]);
  long mg = (long)ceil(R / (B * q.h0 / (1 << q.Lm))) + 1;
  for (long jj = (long)by[c] * s - mg; jj < (long)(by[c] + 1) * s + mg; jj++)
    for (long ii = (long)bx[c] * s - mg; ii < (long)(bx[c] + 1) * s + mg; ii++) {
      long wx = fdivl(ii, nf), wy = fdivl(jj, nf);
      if (wx < -1 || wx > 1 || wy < -1 || wy > 1)
        continue;
      int e = t.fine[(jj - wy * nf) * nf + ii - wx * nf], kk = e * 9 + (int)((wx + 1) * 3 + wy + 1);
      if (kk == last)
        continue;
      last = kk;
      int lo = 0, hi = nk;
      while (lo < hi) {
        int md = (lo + hi) / 2;
        if (key[md] < kk)
          lo = md + 1;
        else
          hi = md;
      }
      if (lo < nk && key[lo] == kk)
        continue;
      if (nk >= CL) {
        atomicMax(err, 6);
        return;
      }
      for (int a = nk; a > lo; a--)
        key[a] = key[a - 1];
      key[lo] = kk;
      nk++;
    }
  int n = 0;
  for (int a = 0; a < nk; a++) {
    int e = key[a] / 9, sc = key[a] % 9;
    if (!cnear(l, bx[c], by[c], lev[e], bx[e], by[e], sc, R, q.h0, q.L))
      continue;
    int in = 0, pse = 1;
    for (int k = t.uo[c]; k < t.uo[c + 1]; k++)
      in |= t.ul[k] == e && t.us[k] == sc;
    for (int k = t.xo[c]; k < t.xo[c + 1]; k++)
      in |= t.xl[k] == e && t.xs[k] == sc;
    for (int k = t.ro[c]; k < t.ro[c + 1]; k++)
      pse &= !(t.rl[k] == e && t.rs[k] == sc);
    if (in)
      continue;
    if (fill) {
      t.cl[t.co[c] + n] = e;
      t.cs[t.co[c] + n] = sc + 9 * pse;
    }
    n++;
  }
  if (!fill)
    cc[c] = n;
}
/* PoU weights and core widths: one CUDA block per mesh block, thread 0
   collects the coarser candidate blocks (fine lattice around the block) */
__global__ static void kpou(PP q, const int *lev, const int *bx, const int *by, const int *fine, double *P,
                            double *D, double *dmx, int *err) {
  __shared__ int cand[PC];
  __shared__ int nc;
  __shared__ double red[NB];
  int b = blockIdx.x, k = threadIdx.x, l = lev[b];
  if (k == 0) {
    int n = 0, ovf = 0;
    if (l > 0) {
      long nf = (long)q.nb0 << q.Lm, s = 1L << (q.Lm - l);
      long mg = (long)ceil(q.Rm / (B * q.h0 / (1 << q.Lm))) + 1;
      int last = -1;
      for (long jj = (long)by[b] * s - mg; jj < (long)(by[b] + 1) * s + mg; jj++)
        for (long ii = (long)bx[b] * s - mg; ii < (long)(bx[b] + 1) * s + mg; ii++) {
          int c = fine[((jj % nf + nf) % nf) * nf + (ii % nf + nf) % nf];
          if (c == last || lev[c] >= l)
            continue;
          last = c;
          int lo = 0, hi = n;
          while (lo < hi) {
            int md = (lo + hi) / 2;
            if (cand[md] < c)
              lo = md + 1;
            else
              hi = md;
          }
          if (lo < n && cand[lo] == c)
            continue;
          if (n >= PC) {
            ovf = 1;
            continue;
          }
          for (int a = n; a > lo; a--)
            cand[a] = cand[a - 1];
          cand[lo] = c;
          n++;
        }
      int w = 0;
      for (int e = 0; e < n; e++)
        if (pcand(l, bx[b], by[b], lev[cand[e]], bx[cand[e]], by[cand[e]], q))
          cand[w++] = cand[e];
      n = w;
    }
    if (ovf)
      atomicMax(err, 6);
    nc = n;
  }
  __syncthreads();
  int i = k % B, j = k / B;
  long n = (long)b * NB + k;
  double hb = XD(q.h0, (double)(1 << l)), Pv, Dv;
  double x = XM((double)(bx[b] * B + i), hb), y = XM((double)(by[b] * B + j), hb);
  pnode(l, x, y, (long)bx[b] * B + i, (long)by[b] * B + j, cand, nc, lev, bx, by, q, &Pv, &Dv);
  P[n] = Pv;
  D[n] = Dv;
  red[k] = Dv;
  __syncthreads();
  for (int h = NB / 2; h > 0; h /= 2) {
    if (k < h)
      red[k] = fmax(red[k], red[k + h]);
    __syncthreads();
  }
  if (k == 0)
    dmx[b] = fmax(0.0, red[0]);
}
/* old -> new values: copy/injection from the old block that covers the new
   node, or DD8 prediction (k levels) from the old lab */
__global__ static void kxfer(long N, int Lm, int nb0, int W, int G, const int *lev, const int *bx, const int *by,
                             const int *ofine, const int *olev, const int *obx, const int *oby, const double *olab,
                             double *o, int *err) {
  long n = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (n >= N)
    return;
  int b = n / NB, i = n % NB % B, j = n % NB / B, l = lev[b], lo = -1, ob = -1;
  long gx = (long)bx[b] * B + i, gy = (long)by[b] * B + j;
  for (int q = 0; q <= Lm && lo < 0; q++) {
    long X = q >= l ? gx << (q - l) : gx >> (l - q), Y = q >= l ? gy << (q - l) : gy >> (l - q);
    int p = lmapd(q, X / B, Y / B, nb0, Lm, ofine, olev);
    if (p >= 0) {
      lo = q;
      ob = p;
    }
  }
  if (lo < 0) {
    atomicMax(err, 1);
    return;
  }
  double v;
  if (lo >= l) {
    long X = gx << (lo - l), Y = gy << (lo - l);
    v = olab[(long)ob * W * W + (Y % B + G) * W + X % B + G];
  } else {
    int k = l - lo, e = 0;
    long X = gx - ((long)obx[ob] * B << k), Y = gy - ((long)oby[ob] * B << k);
    v = predv(olab, W, G, ob, k, X, Y, &e);
    if (e)
      atomicMax(err, e);
  }
  o[n] = v;
}
__global__ static void kdiag(long N, const double *a, const double *o, const double *u, const double *v,
                             double *out) {
  __shared__ double s[5][1024];
  double c = 0, ca = 0, z = 0, e = 0, wm = 0;
  for (long n = threadIdx.x; n < N; n += blockDim.x) {
    double ao = a[n], oo = o[n];
    c += oo * ao;
    ca += fabs(oo * ao);
    z += 0.5 * oo * oo * ao;
    e += 0.5 * (u[n] * u[n] + v[n] * v[n]) * ao;
    wm = fmax(wm, fabs(oo));
  }
  int k = threadIdx.x;
  s[0][k] = c;
  s[1][k] = ca;
  s[2][k] = z;
  s[3][k] = e;
  s[4][k] = wm;
  __syncthreads();
  for (int h = blockDim.x / 2; h > 0; h /= 2) {
    if (k < h) {
      for (int q = 0; q < 4; q++)
        s[q][k] += s[q][k + h];
      s[4][k] = fmax(s[4][k], s[4][k + h]);
    }
    __syncthreads();
  }
  if (k == 0)
    for (int q = 0; q < 5; q++)
      out[q] = s[q][0];
}
static void *stmp;
static size_t stn;
static void stgrow(size_t need) {
  if (need > stn) {
    dfree(stmp);
    stmp = dmalloc(need);
    stn = need;
  }
}
static void xscan(const int *in, int *out, long n) {
  size_t need = 0;
  ck(cub::DeviceScan::ExclusiveSum(NULL, need, in, out, (int)n, cudaStreamPerThread), "ExclusiveSum size");
  stgrow(need);
  ck(cub::DeviceScan::ExclusiveSum(stmp, need, in, out, (int)n, cudaStreamPerThread), "ExclusiveSum");
}
static int dget(const int *p) {
  int v;
  ck(cudaMemcpy(&v, p, sizeof v, cudaMemcpyDeviceToHost), "cudaMemcpy");
  return v;
}
/* counts c[0, n) -> offsets o[0, n]; returns the total (one int to the host) */
static int offs(int *c, int *o, long n) {
  ck(cudaMemset(c + n, 0, sizeof(int)), "cudaMemset");
  xscan(c, o, n + 1);
  return dget(o + n);
}
static void chkerr(const char *where) {
  static const char *msg[] = {"", "regrid: no old block covers a new node", "transfer stencil outside the lab",
                              "mesh is not 2:1 balanced", "neighbour case not in the case table",
                              "list buffer overflow", "candidate list overflow"};
  int e = dget(derr);
  if (e) {
    fprintf(stderr, "vm2: error: %s: %s\n", where, e < 7 ? msg[e] : "?");
    exit(1);
  }
}
static void labf(const double *v, double *lb) {
  labown<<<grid(N), nthread>>>(N, W, G, v, lb);
  if (mode == 1)
    labfill<<<grid(gh.nr), nthread>>>(gh.nr, dgr, dgp, dgc, dgw, v, lb);
  else
    labcase<<<grid((long)(N / NB) * NG), nthread>>>((long)(N / NB) * NG, NG, dcas, dslot, dcrow, dccol, dcw, dgpos,
                                                    W * W, v, lb);
  launched("lab");
}
/* GPU: target levels tg (fine lattice) -> leaves (Morton order), fine map */
static void dbuild(int *tg) {
  int nb0 = n0 / B;
  long nf = (long)nb0 << Lm, nf2 = nf * nf;
  Off po;
  po.o[0] = 0;
  for (int l = 0; l < Lm; l++)
    po.o[l + 1] = po.o[l] + ((long)nb0 << l) * ((long)nb0 << l);
  int *lv = (int *)dmalloc(nf2 * sizeof(int)), *pyr = (int *)dmalloc((po.o[Lm] + 1) * sizeof(int)),
      *chg = (int *)dmalloc(sizeof(int));
  unsigned long long *tw = (unsigned long long *)dmalloc(nf2 * sizeof(unsigned long long));
  for (;;) {
    for (int l = Lm - 1; l >= 0; l--) {
      long n = (long)nb0 << l;
      kpyr<<<grid(n * n), nthread>>>(n, l + 1 == Lm ? tg : pyr + po.o[l + 1], pyr + po.o[l]);
    }
    klv<<<grid(nf2), nthread>>>(nf, nb0, Lm, po, pyr, lv);
    ck(cudaMemset(chg, 0, sizeof(int)), "cudaMemset");
    ck(cudaMemset(tw, 0, nf2 * sizeof(unsigned long long)), "cudaMemset");
    kchk<<<grid(nf2), nthread>>>(nf, Lm, nbuf, lv, tw, chg);
    kapply<<<grid(nf2), nthread>>>(nf2, tw, tg);
    launched("balance");
    if (!dget(chg))
      break;
  }
  int *fl = (int *)dmalloc((nf2 + 1) * sizeof(int)), *ix = (int *)dmalloc((nf2 + 1) * sizeof(int));
  kleaf<<<grid(nf2), nthread>>>(nf, Lm, lv, fl);
  g.nleaf = offs(fl, ix, nf2);
  dlev = (int *)dmalloc(g.nleaf * sizeof(int));
  dbx = (int *)dmalloc(g.nleaf * sizeof(int));
  dby = (int *)dmalloc(g.nleaf * sizeof(int));
  t.fine = (int *)dmalloc(nf2 * sizeof(int));
  kleafw<<<grid(nf2), nthread>>>(nf, Lm, lv, fl, ix, dlev, dbx, dby);
  kfine<<<grid(nf2), nthread>>>(nf, Lm, lv, ix, t.fine);
  launched("leaves");
  dfree(lv);
  dfree(pyr);
  dfree(chg);
  dfree(tw);
  dfree(fl);
  dfree(ix);
}
/* GPU: neighbour cases, FMM tree and lists, weights, C lists */
static void dmesh(void) {
  int nl = g.nleaf, nb0 = n0 / B, D0 = m.D0, Dm = m.Dm, nbox;
  long nf = (long)nb0 << Lm, tot = bo.o[Dm + 1];
  N = (long)nl * NB;
  drel = (int *)dmalloc(nl * 9 * sizeof(int));
  dslot = (int *)dmalloc((long)nl * NS * sizeof(int));
  dcas = (int *)dmalloc(nl * sizeof(int));
  krel<<<grid(nl), nthread>>>(nl, nb0, Lm, dlev, dbx, dby, t.fine, dcidx, drel, dslot, dcas, derr);
  launched("krel");
  int *ex = (int *)dmalloc((tot + 1) * sizeof(int)), *bix = (int *)dmalloc((tot + 1) * sizeof(int)),
      *lh = (int *)dmalloc((DM + 3) * sizeof(int));
  kbox<<<grid(tot), nthread>>>(tot, bo, Dm, D0, nf, t.fine, dlev, ex);
  g.nbox = nbox = offs(ex, bix, tot);
  klohi<<<1, 32>>>(Dm + 2, bo, bix, lh);
  ck(cudaMemcpy(g.lo, lh, (Dm + 2) * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
  t.bdep = (int *)dmalloc(nbox * sizeof(int));
  t.bpx = (int *)dmalloc(nbox * sizeof(int));
  t.bpy = (int *)dmalloc(nbox * sizeof(int));
  t.bch = (int *)dmalloc(nbox * sizeof(int));
  t.bpar = (int *)dmalloc(nbox * sizeof(int));
  dbleaf = (int *)dmalloc(nbox * sizeof(int));
  t.lbox = (int *)dmalloc(nl * sizeof(int));
  dbmap = (int *)dmalloc(tot * sizeof(int));
  kboxw<<<grid(tot), nthread>>>(tot, bo, Dm, D0, nf, t.fine, dlev, ex, bix, t.bdep, t.bpx, t.bpy, t.bch, t.bpar,
                                dbleaf, t.lbox, dbmap);
  launched("kboxw");
  dfree(ex);
  dfree(bix);
  dfree(lh);
  t.cx = (double *)dmalloc(nbox * sizeof(double));
  t.cy = (double *)dmalloc(nbox * sizeof(double));
  kcen<<<grid(nbox), nthread>>>(nbox, L, t.bdep, t.bpx, t.bpy, t.cx, t.cy);
  TV v;
  v.Dm = Dm;
  v.nleaf = nl;
  v.nbox = nbox;
  v.bd = t.bdep;
  v.px = t.bpx;
  v.py = t.bpy;
  v.ch = t.bch;
  v.par = t.bpar;
  v.leaf = dbleaf;
  v.lbox = t.lbox;
  v.bmap = dbmap;
  v.bo = bo;
  int *c1 = (int *)dmalloc((nbox + nl + 1) * sizeof(int)), *c2 = (int *)dmalloc((nbox + nl + 1) * sizeof(int));
  t.uo = (int *)dmalloc((nl + 1) * sizeof(int));
  t.wo = (int *)dmalloc((nl + 1) * sizeof(int));
  kul<<<grid2(nl, 128), 128>>>(v, 0, c1, c2, NULL, NULL, NULL, NULL, NULL, NULL, derr);
  g.nu = offs(c1, t.uo, nl);
  g.nw = offs(c2, t.wo, nl);
  t.ul = (int *)dmalloc((g.nu + 1) * sizeof(int));
  t.us = (int *)dmalloc((g.nu + 1) * sizeof(int));
  t.wl = (int *)dmalloc((g.nw + 1) * sizeof(int));
  t.ws = (int *)dmalloc((g.nw + 1) * sizeof(int));
  kul<<<grid2(nl, 128), 128>>>(v, 1, NULL, NULL, t.uo, t.ul, t.us, t.wo, t.wl, t.ws, derr);
  launched("kul");
  int *x0o = (int *)dmalloc((nbox + 1) * sizeof(int));
  t.vo = (int *)dmalloc((nbox + 1) * sizeof(int));
  kvx<<<grid2(nbox, 128), 128>>>(v, 0, c1, c2, NULL, NULL, NULL, NULL, NULL, NULL, derr);
  g.nv = offs(c1, t.vo, nbox);
  int nx0 = offs(c2, x0o, nbox);
  t.vl = (int *)dmalloc((g.nv + 1) * sizeof(int));
  t.vs = (int *)dmalloc((g.nv + 1) * sizeof(int));
  int *x0l = (int *)dmalloc((nx0 + 1) * sizeof(int)), *x0s = (int *)dmalloc((nx0 + 1) * sizeof(int));
  kvx<<<grid2(nbox, 128), 128>>>(v, 1, NULL, NULL, t.vo, t.vl, t.vs, x0o, x0l, x0s, derr);
  launched("kvx");
  t.xo = (int *)dmalloc((nl + 1) * sizeof(int));
  kxl<<<grid(nl), nthread>>>(v, 0, x0o, x0l, x0s, c1, NULL, NULL, NULL, NULL);
  g.nx = offs(c1, t.xo, nl);
  t.xl = (int *)dmalloc((g.nx + 1) * sizeof(int));
  t.xs = (int *)dmalloc((g.nx + 1) * sizeof(int));
  int *own = (int *)dmalloc((g.nx + 1) * sizeof(int)), *sk = (int *)dmalloc((g.nx + 1) * sizeof(int)),
      *sv = (int *)dmalloc((g.nx + 1) * sizeof(int)), *iv = (int *)dmalloc((g.nx + 1) * sizeof(int));
  kxl<<<grid(nl), nthread>>>(v, 1, x0o, x0l, x0s, NULL, t.xo, t.xl, t.xs, own);
  launched("kxl");
  t.ro = (int *)dmalloc((nl + 1) * sizeof(int));
  t.rl = (int *)dmalloc((g.nx + 1) * sizeof(int));
  t.rs = (int *)dmalloc((g.nx + 1) * sizeof(int));
  if (g.nx > 0) {
    int kb = 1;
    while ((1L << kb) <= nl)
      kb++;
    iotai<<<grid(g.nx), nthread>>>(g.nx, iv);
    size_t need = 0;
    ck(cub::DeviceRadixSort::SortPairs(NULL, need, t.xl, sk, iv, sv, g.nx, 0, kb, cudaStreamPerThread), "SortPairs size");
    stgrow(need);
    ck(cub::DeviceRadixSort::SortPairs(stmp, need, t.xl, sk, iv, sv, g.nx, 0, kb, cudaStreamPerThread), "SortPairs");
    kr<<<grid(g.nx), nthread>>>(g.nx, sv, own, t.xs, t.rl, t.rs);
  }
  starts<<<grid(nl + 1), nthread>>>(g.nx, nl, sk, t.ro);
  launched("R lists");
  void *tmp[] = {x0o, x0l, x0s, own, sk, sv, iv};
  for (size_t k = 0; k < sizeof tmp / sizeof *tmp; k++)
    dfree(tmp[k]);
  t.N = N;
  t.p = p;
  t.L = L;
  t.nleaf = nl;
  t.nbox = nbox;
  t.Dm = Dm;
  t.Fn = 1 << Dm;
  t.cd = cd;
  t.nu = nu;
  talloc(&t, N, nl, nbox, p);
  t.ov = (double *)dmalloc(N * sizeof(double));
  t.oh = (double *)dmalloc(N * sizeof(double));
  dmxd = (double *)dmalloc(nl * sizeof(double));
  kpou<<<nl, NB>>>(pp, dlev, dbx, dby, t.fine, t.ov, t.oh, dmxd, derr);
  launched("kpou");
  t.co = (int *)dmalloc((nl + 1) * sizeof(int));
  kcl<<<grid2(nl, 128), 128>>>(t, pp, dlev, dbx, dby, dmxd, 0, c1, derr);
  g.nc = offs(c1, t.co, nl);
  t.cl = (int *)dmalloc((g.nc + 1) * sizeof(int));
  t.cs = (int *)dmalloc((g.nc + 1) * sizeof(int));
  kcl<<<grid2(nl, 128), 128>>>(t, pp, dlev, dbx, dby, dmxd, 1, c1, derr);
  launched("kcl");
  dfree(c1);
  dfree(c2);
  common();
  if (mode == 0) {
    m.nleaf = nl;
    m.nbox = nbox;
    for (int d = 0; d <= Dm + 1; d++) {
      m.lo[d] = g.lo[d];
      m.hi[d] = d <= Dm ? g.lo[d + 1] : nbox;
    }
  }
  chkerr("mesh");
}
/* -X: compare the host reference with the GPU mesh, bitwise */
static long nbad, ncmp;
static void cmpi(const char *s, const int *h, const int *d, long n) {
  if (n <= 0)
    return;
  int *a = (int *)hmalloc(n * sizeof(int));
  down(a, d, n * sizeof(int));
  long k = 0;
  while (k < n && a[k] == h[k])
    k++;
  if (k < n) {
    fprintf(stderr, "vm2: check: %s differs at %ld of %ld (host %d, gpu %d)\n", s, k, n, h[k], a[k]);
    nbad++;
  }
  free(a);
}
static void cmpd(const char *s, const double *h, const double *d, long n) {
  double *a = (double *)hmalloc(n * sizeof(double));
  down(a, d, n * sizeof(double));
  long k = 0, nd = 0;
  double mx = 0;
  for (long q = 0; q < n; q++)
    if (memcmp(&a[q], &h[q], sizeof(double))) {
      if (!nd)
        k = q;
      nd++;
      mx = fmax(mx, fabs(a[q] - h[q]));
    }
  if (nd) {
    fprintf(stderr, "vm2: check: %s: %ld of %ld differ (first %ld: host %.17g gpu %.17g, max diff %.3e)\n", s, nd,
            n, k, h[k], a[k], mx);
    nbad++;
  }
  free(a);
}
static void meshcmp(void) {
  int nl = m.nleaf, WW = W * W;
  long nbad0 = nbad;
  ncmp++;
  if (g.nleaf != m.nleaf || g.nbox != m.nbox) {
    fprintf(stderr, "vm2: check: leaves %d / %d, boxes %d / %d (host / gpu)\n", m.nleaf, g.nleaf, m.nbox, g.nbox);
    nbad++;
    return;
  }
  cmpi("lev", m.lev, dlev, nl);
  cmpi("bx", m.bx, dbx, nl);
  cmpi("by", m.by, dby, nl);
  int *hf2 = hfine(&m);
  cmpi("fine", hf2, t.fine, (1L << m.Dm) * (1L << m.Dm));
  free(hf2);
  cmpi("rel", m.rel, drel, (long)nl * 9);
  for (int d = 0; d <= m.Dm; d++)
    if (m.lo[d] != g.lo[d] || m.hi[d] != g.lo[d + 1]) {
      fprintf(stderr, "vm2: check: depth %d boxes [%d %d) / [%d %d)\n", d, m.lo[d], m.hi[d], g.lo[d], g.lo[d + 1]);
      nbad++;
    }
  cmpi("bd", m.bd, t.bdep, m.nbox);
  cmpi("px", m.px, t.bpx, m.nbox);
  cmpi("py", m.py, t.bpy, m.nbox);
  cmpi("ch", m.ch, t.bch, m.nbox);
  cmpi("par", m.par, t.bpar, m.nbox);
  cmpi("leaf", m.leaf, dbleaf, m.nbox);
  cmpi("lbox", m.lbox, t.lbox, nl);
  const char *nm[] = {"U", "W", "X", "R", "C", "V"};
  int *ho[] = {m.uo, m.wo, m.xo, m.ro, m.co, m.vo}, *hl[] = {m.ul, m.wl, m.xl, m.rl, m.cl, m.vl},
      *hs[] = {m.us, m.ws, m.xs, m.rs, m.cs, m.vs}, *dof[] = {t.uo, t.wo, t.xo, t.ro, t.co, t.vo},
      *dl[] = {t.ul, t.wl, t.xl, t.rl, t.cl, t.vl}, *ds[] = {t.us, t.ws, t.xs, t.rs, t.cs, t.vs};
  for (int k = 0; k < 6; k++) {
    int n = k < 5 ? nl : m.nbox;
    char s[16];
    snprintf(s, sizeof s, "%s offsets", nm[k]);
    cmpi(s, ho[k], dof[k], n + 1);
    snprintf(s, sizeof s, "%s list", nm[k]);
    cmpi(s, hl[k], dl[k], ho[k][n]);
    snprintf(s, sizeof s, "%s shift", nm[k]);
    cmpi(s, hs[k], ds[k], ho[k][n]);
  }
  cmpd("PoU weight", hh, t.ov, N);
  cmpd("core width", hgr, t.oh, N);
  int *cs = (int *)hmalloc(nl * sizeof(int)), *sl = (int *)hmalloc((long)nl * NS * sizeof(int));
  down(cs, dcas, nl * sizeof(int));
  down(sl, dslot, (long)nl * NS * sizeof(int));
  long bad = 0;
  for (int b = 0; b < nl; b++)
    for (int q = 0; q < NG; q++) {
      long r = (long)b * NG + q, e0 = gh.p[r], e1 = gh.p[r + 1];
      long f0 = hcrow[(long)cs[b] * (NG + 1) + q], f1 = hcrow[(long)cs[b] * (NG + 1) + q + 1];
      int ok = gh.row[r] == (long)b * WW + hgpos[q] && e1 - e0 == f1 - f0;
      for (long k = 0; ok && k < e1 - e0; k++) {
        int y = hccol[f0 + k];
        ok = sl[b * NS + y / NB] * NB + y % NB == gh.c[e0 + k] && !memcmp(&hcw[f0 + k], &gh.w[e0 + k], sizeof(double));
      }
      if (!ok && !bad++)
        fprintf(stderr, "vm2: check: ghost row %d of block %d (case %d) differs\n", q, b, cs[b]);
    }
  if (bad)
    nbad++;
  free(cs);
  free(sl);
  if (nbad == nbad0)
    fprintf(stderr, "vm2: check %ld: %d blocks, %d boxes, %ld ghost rows, lists, weights: identical\n", ncmp, nl,
            m.nbox, gh.nr);
}
/* host reference: refinement flags and fine-lattice targets */
static int *htarget(const double *db, const double *sb, double wm, double er, double kap, double tau) {
  int nl0 = m.nleaf;
  long nf = (long)m.nb0 << Lm;
  int *tg = (int *)hmalloc(nf * nf * sizeof(int)), *nlv = (int *)hmalloc(nl0 * sizeof(int)),
      *fl = (int *)hmalloc(nl0 * sizeof(int)), *nlv0 = (int *)hmalloc(nl0 * sizeof(int));
  double ec = er / 256;
  for (int b = 0; b < nl0; b++) {
    double d = db[b] / wm * exp(kap * sb[b] * tau);
    int l = m.lev[b];
    fl[b] = d > er;
    nlv0[b] = d > er && l < Lm ? l + 1 : (d < ec && l > 0 ? l - 1 : l);
  }
  for (int b = 0; b < nl0; b++) {
    int l = m.lev[b], px = m.bx[b] & ~1, py = m.by[b] & ~1, ok = 1;
    nlv[b] = nlv0[b];
    if (nlv0[b] >= l)
      continue;
    for (int q = 0; q < 4; q++) {
      int *c = lmap(&m, l, px + (q & 1), py + (q >> 1));
      ok &= c != NULL && *c >= 0 && nlv0[*c] < l;
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
  free(nlv);
  free(nlv0);
  free(fl);
  return tg;
}
/* host reference: old -> new values (new mesh m, old mesh mo, old lab ol) */
static double *htransfer(H *mo, const double *ol) {
  double *o = (double *)hmalloc(N * sizeof(double));
  for (long n = 0; n < N; n++) {
    int b = n / NB, i = n % NB % B, j = n % NB / B, l = m.lev[b], *p, lo = -1, ob = -1, err = 0;
    long gx = (long)m.bx[b] * B + i, gy = (long)m.by[b] * B + j;
    for (int q = 0; q <= mo->Lm && lo < 0; q++) {
      long X = q >= l ? gx << (q - l) : gx >> (l - q), Y = q >= l ? gy << (q - l) : gy >> (l - q);
      if ((p = lmap(mo, q, X / B, Y / B)) != NULL && *p >= 0) {
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
      o[n] = ol[(long)ob * W * W + (Y % B + G) * W + X % B + G];
    } else {
      int k = l - lo;
      long X = gx - ((long)mo->bx[ob] * B << k), Y = gy - ((long)mo->by[ob] * B << k);
      o[n] = predv(ol, W, G, ob, k, X, Y, &err);
    }
    if (err) {
      fprintf(stderr, "vm2: error: transfer stencil outside the lab\n");
      exit(1);
    }
  }
  return o;
}
static void regrid(double er, double kap, double tau, int ic, double rho, double eps) {
  int nl0 = m.nleaf, nb0 = n0 / B;
  long nf = (long)nb0 << Lm, nf2 = nf * nf;
  double *ddb = (double *)dmalloc(nl0 * sizeof(double)), *dsb = (double *)dmalloc(nl0 * sizeof(double)),
         *dwm = (double *)dmalloc(sizeof(double));
  labf(s0.g, lab);
  detail<<<grid(N), nthread>>>(N, W, G, lab, wk[0].y);
  bmax<<<nl0, NB>>>(wk[0].y, ddb);
  amax<<<1, 1024>>>(N, s0.g, dwm);
  launched("detail");
  if (kap > 0) {
    rhs(&t, &m, 0, s0, u);
    labf(u.x, lab2);
    labf(u.y, lab3);
    strain<<<grid(N), nthread>>>(N, W, G, h0, dlev, lab2, lab3, wk[0].g);
    bmax<<<nl0, NB>>>(wk[0].g, dsb);
    launched("strain");
  } else
    ck(cudaMemset(dsb, 0, nl0 * sizeof(double)), "cudaMemset");
  int *dtg = NULL, *htg = NULL;
  if (mode != 1) {
    int *nlv = (int *)dmalloc(nl0 * sizeof(int)), *nlv2 = (int *)dmalloc(nl0 * sizeof(int)),
        *fl = (int *)dmalloc(nl0 * sizeof(int));
    kflag<<<grid(nl0), nthread>>>(nl0, Lm, er, kap, tau, dlev, ddb, dsb, dwm, nlv, fl);
    kcoarse<<<grid(nl0), nthread>>>(nl0, nb0, Lm, dlev, dbx, dby, t.fine, nlv, nlv2);
    dtg = (int *)dmalloc(nf2 * sizeof(int));
    ktg0<<<grid(nf2), nthread>>>(nf2, t.fine, nlv2, dtg);
    ktg1<<<grid(nl0), nthread>>>(nl0, nb0, Lm, nbuf, dlev, dbx, dby, fl, nlv2, dtg);
    launched("targets");
    dfree(nlv);
    dfree(nlv2);
    dfree(fl);
  }
  if (mode != 0) {
    double *db = (double *)hmalloc(nl0 * sizeof(double)), *sb = (double *)hmalloc(nl0 * sizeof(double)), wm;
    down(db, ddb, nl0 * sizeof(double));
    down(sb, dsb, nl0 * sizeof(double));
    down(&wm, dwm, sizeof(double));
    htg = htarget(db, sb, wm, er, kap, tau);
    if (mode == 2)
      cmpi("targets", htg, dtg, nf2);
    free(db);
    free(sb);
  }
  dfree(ddb);
  dfree(dsb);
  dfree(dwm);
  double *holab = NULL;
  if (mode != 0 && !ic) {
    holab = (double *)hmalloc((size_t)nl0 * W * W * sizeof(double));
    down(holab, lab, (size_t)nl0 * W * W * sizeof(double));
  }
  double *olab = lab;
  int *ofine = t.fine, *olev = dlev, *obx = dbx, *oby = dby;
  lab = NULL;
  t.fine = dlev = dbx = dby = NULL;
  H mo = m;
  CSR go = gh;
  meshfree();
  if (mode == 1) {
    hmesh(htg);
    hupload();
  } else {
    dbuild(dtg);
    dmesh();
    if (mode == 2) {
      hmesh(htg);
      meshcmp();
    }
  }
  if (ic) {
    icw<<<grid(N), nthread>>>(N, pair_, rho, eps, s0.x, s0.y, s0.g);
    launched("icw");
  } else {
    if (mode != 1) {
      kxfer<<<grid(N), nthread>>>(N, Lm, nb0, W, G, dlev, dbx, dby, ofine, olev, obx, oby, olab, s0.g, derr);
      launched("kxfer");
      chkerr("transfer");
    }
    if (mode != 0) {
      double *ht = htransfer(&mo, holab);
      if (mode == 1)
        up(s0.g, ht, N * sizeof(double));
      else
        cmpd("transfer", ht, s0.g, N);
      free(ht);
    }
    sums3<<<1, 1024>>>(N, t.ov, s0.g, s0.g, dsm2);
    gfix<<<grid(N), nthread>>>(N, rgsum, dsm2, s0.g);
    launched("gfix");
  }
  dfree(olab);
  dfree(ofine);
  dfree(olev);
  dfree(obx);
  dfree(oby);
  dfree(dtg);
  free(htg);
  free(holab);
  if (mode != 0)
    hfree(&mo, &go);
}
/* host copies of the mesh for diagnostics and the final error (not in the step) */
static int *vlev, *vbx, *vby, *vrel, vown;
static double *vh, *vx, *vy;
static void hostsync(void) {
  if (vown) {
    void *h[] = {vlev, vbx, vby, vrel, vh, vx, vy};
    for (size_t k = 0; k < sizeof h / sizeof *h; k++)
      free(h[k]);
    vown = 0;
  }
  if (mode == 1) {
    vlev = m.lev;
    vbx = m.bx;
    vby = m.by;
    vrel = m.rel;
    vh = hh;
    vx = hx;
    vy = hy;
    return;
  }
  int nl = N / NB;
  vlev = (int *)hmalloc(nl * sizeof(int));
  vbx = (int *)hmalloc(nl * sizeof(int));
  vby = (int *)hmalloc(nl * sizeof(int));
  vrel = (int *)hmalloc(nl * 9 * sizeof(int));
  vh = (double *)hmalloc(N * sizeof(double));
  vx = (double *)hmalloc(N * sizeof(double));
  vy = (double *)hmalloc(N * sizeof(double));
  down(vlev, dlev, nl * sizeof(int));
  down(vbx, dbx, nl * sizeof(int));
  down(vby, dby, nl * sizeof(int));
  down(vrel, drel, nl * 9 * sizeof(int));
  down(vh, t.ov, N * sizeof(double));
  down(vx, s0.x, N * sizeof(double));
  down(vy, s0.y, N * sizeof(double));
  vown = 1;
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
              "[-i <reference> -I <reference N>] [-c (velocity check)] [-R <regrid interval> -E <eps_r> -K <c3 factor>] "
              "[-H (host regrid, reference)] [-X (GPU regrid checked against the host)]\n");
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
    if (argv[0][1] == 'H') {
      mode = 1;
      continue;
    }
    if (argv[0][1] == 'X') {
      mode = 2;
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
  NG = W * W - NB;
  hgpos = (int *)hmalloc(NG * sizeof(int));
  for (int j = -G, q = 0; j < B + G; j++)
    for (int i = -G; i < B + G; i++)
      if (i < 0 || i >= B || j < 0 || j >= B)
        hgpos[q++] = (j + G) * W + i + G;
  ptime = getenv("VM2TIME") != NULL;
  for (int k = 0; k < 8; k++)
    ck(cudaEventCreate(&pev[k]), "cudaEventCreate");
  h0 = L / n0;
  if (p + 1 > tab.Pt) {
    fprintf(stderr, "vm2: error: -p %d exceeds the table (%d)\n", p, tab.Pt - 1);
    exit(1);
  }
  if (mode != 1 && (chk || tst || verb || lts)) {
    fprintf(stderr, "vm2: -c, -T, -v and -l use the host mesh path (-H)\n");
    mode = 1;
  }
  {
    int nb0 = n0 / B, D0 = 0;
    while ((1 << D0) < nb0)
      D0++;
    if ((1 << D0) != nb0 || D0 + Lm > DM) {
      fprintf(stderr, "vm2: error: -n / %d must be a power of two, depth <= %d\n", B, DM);
      exit(1);
    }
    if (mode != 1 && nb0 < 4) {
      fprintf(stderr, "vm2: error: the GPU regrid needs -n >= %d (use -H)\n", 4 * B);
      exit(1);
    }
    m.nb0 = nb0;
    m.Lm = Lm;
    m.D0 = D0;
    m.Dm = D0 + Lm;
    bo.o[0] = 0;
    for (int d = 0; d <= m.Dm; d++)
      bo.o[d + 1] = bo.o[d] + (1L << (2 * d));
  }
  pp.h0 = h0;
  pp.L = L;
  pp.ramp = ramp;
  pp.band = band;
  pp.cd = cd;
  pp.Rm = ramp * h0 + (8 + band) * h0 / 2;
  pp.nb0 = n0 / B;
  pp.Lm = Lm;
  tables(tab.F, tab.Pt);
  derr = (int *)dmalloc(sizeof(int));
  ck(cudaMemset(derr, 0, sizeof(int)), "cudaMemset");
  dpair = (unsigned long long *)dmalloc(sizeof(unsigned long long));
  ck(cudaMemset(dpair, 0, sizeof(unsigned long long)), "cudaMemset");
  dsm2 = (double *)dmalloc(3 * sizeof(double));
  double *ddg = (double *)dmalloc(5 * sizeof(double)), rgt = 0, rgt0 = 0;
  if (mode != 1) {
    double c0 = now();
    cases(&tab, !nosm);
    fprintf(stderr, "vm2: %d neighbour cases, %ld ghost entries, %.2f s\n", ncase, cnnz, now() - c0);
  }
  {
    long nf2 = ((long)(n0 / B) << Lm) * ((long)(n0 / B) << Lm);
    int *tg0;
    if (rg) {
      tg0 = (int *)hmalloc(nf2 * sizeof(int));
      for (long k = 0; k < nf2; k++)
        tg0[k] = 0;
    } else
      tg0 = bandtg(n0 / B, Lm, wband);
    if (mode == 1) {
      hmesh(tg0);
      hupload();
    } else {
      int *dtg = idup(tg0, nf2);
      dbuild(dtg);
      dfree(dtg);
      dmesh();
      if (mode == 2) {
        hmesh(tg0);
        meshcmp();
      }
    }
    free(tg0);
  }
  icw<<<grid(N), nthread>>>(N, pair_, rho, eps, s0.x, s0.y, s0.g);
  launched("icw");
  if (rg)
    for (int k = 0; k <= Lm; k++) {
      ck(cudaDeviceSynchronize(), "sync");
      double r0 = now();
      regrid(er, kap, rg * H0, 1, rho, eps);
      ck(cudaDeviceSynchronize(), "sync");
      rgt0 += now() - r0;
    }
  {
    double sm[3];
    sums3<<<1, 1024>>>(N, t.ov, s0.g, s0.g, dsm2);
    down(sm, dsm2, sizeof sm);
    rgsum = sm[1];
  }
  hostsync();
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
      a += vh[n];
      double r = vh[n] * (n0 << vlev[n / NB]) * (n0 << vlev[n / NB]) / (L * L);
      mn = fmin(mn, r);
      mx = fmax(mx, r);
      qa += vh[n] * cos(vx[n] + 0.3) * cos(vx[n] + 0.3) * sin(vy[n]) * sin(vy[n]);
      qb += vh[n] * exp(sin(vx[n] + 0.1)) * cos(vy[n] - 0.2) * cos(vy[n] - 0.2);
    }
    fprintf(stderr, "vm2: sum of volumes / area - 1 = %.3e, V/h^2 in [%.6f, %.6f], quadrature errors %.3e %.3e\n",
            a / (L * L) - 1, mn, mx, qa / (pi * pi) - 1, qb / (2 * pi * pi * 1.2660658777520083) - 1);
  }
  if (mode == 1)
    fprintf(stderr, "vm2: blocks %d (levels %d) particles %ld boxes %d ghost rows %ld nnz %ld U %d V %d W %d X %d\n",
            m.nleaf, Lm + 1, N, m.nbox, gh.nr, gh.nnz, m.uo[m.nleaf], m.vo[m.nbox], m.wo[m.nleaf], m.xo[m.nleaf]);
  else {
    int nl = N / NB, *cs = (int *)hmalloc(nl * sizeof(int));
    long nz = 0;
    down(cs, dcas, nl * sizeof(int));
    for (int b = 0; b < nl; b++)
      nz += hcrow[(long)cs[b] * (NG + 1) + NG] - hcrow[(long)cs[b] * (NG + 1)];
    free(cs);
    fprintf(stderr, "vm2: blocks %d (levels %d) particles %ld boxes %d ghost rows %ld nnz %ld U %d V %d W %d X %d\n",
            nl, Lm + 1, N, g.nbox, (long)nl * NG, nz, g.nu, g.nv, g.nw, g.nx);
  }
  if (chk) {
    down(hg, s0.g, N * sizeof(double));
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
  ck(cudaProfilerStart(), "cudaProfilerStart");
  for (long st = 0;; st++) {
    if (st % nd == 0 || lts) {
      rhs(&t, &m, 0, s0, u);
      nfull++;
    }
    if (st % nd == 0) {
      double dg5[5];
      unsigned long long np;
      kdiag<<<1, 1024>>>(N, t.ov, s0.g, u.x, u.y, ddg);
      launched("kdiag");
      down(dg5, ddg, sizeof dg5);
      down(&np, dpair, sizeof np);
      printf("% 10ld % .16e % .16e % .16e % .16e % .16e % .16e %d %ld %ld %ld %ld %ld % .16e %ld\n", st,
             st * H0, dg5[0], dg5[1], dg5[2], dg5[3], dg5[4], nfst, nfull, nslow, nfast, (long)np, nfpair,
             now() - wt, N);
      if (fflush(stdout) != 0 || ferror(stdout)) {
        fprintf(stderr, "vm2: error: fail to write stdout\n");
        exit(1);
      }
    }
    if (st >= nst)
      break;
    if (rg && st > 0 && st % rg == 0) {
      ck(cudaDeviceSynchronize(), "sync");
      double r0 = now();
      regrid(er, kap, rg * H0, 0, rho, eps);
      ck(cudaDeviceSynchronize(), "sync");
      rgt += now() - r0;
      nreg++;
      if (gexe != NULL) {
        ck(cudaGraphExecDestroy(gexe), "GraphExecDestroy");
        gexe = NULL;
      }
    }
    nodet += N;
    if (lts) {
      labf(u.x, lab);
      labf(u.y, lab2);
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
    long c0f = nfull, c0s = nslow, c0a = nfast;
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
      labf(s1.g, lab);
      labf(dd.x, lab2);
      labf(dd.y, lab3);
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
  ck(cudaProfilerStop(), "cudaProfilerStop");
  if (ptime)
    fprintf(stderr, "vm2: time [s] sort+gather %.3f fmm-up+m2l+down %.3f l2p+m2p %.3f p2p %.3f unsort+zbar %.3f\n",
            ptacc[0] / 1e3, ptacc[1] / 1e3, ptacc[2] / 1e3, ptacc[3] / 1e3, ptacc[4] / 1e3);
  fprintf(stderr, "vm2: regrid %s: %d initial calls %.4f s, %ld calls in the run %.4f s (%.3f ms per call), "
          "mean nodes %.0f, last %ld\n", mode == 0 ? "gpu" : (mode == 1 ? "host" : "gpu+host check"),
          rg ? Lm + 1 : 0, rgt0, nreg, rgt, nreg > 0 ? 1e3 * rgt / nreg : 0.0, nst > 0 ? nodet / nst : (double)N, N);
  if (mode == 2)
    fprintf(stderr, "vm2: check: %ld meshes compared, %ld mismatches\n", ncmp, nbad);
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
    hostsync();
    down(hg, s0.g, N * sizeof(double));
    double a2 = 0, b2 = 0, ai = 0, bi = 0, el[DM + 1] = {0}, eg = 0;
    for (long n = 0; n < N; n++) {
      int b = n / NB, i = n % NB % B, j = n % NB / B, s = rn / (n0 << vlev[b]);
      double f = ref[(long)(vby[b] * B + j) * s * rn + (long)(vbx[b] * B + i) * s],
             o = hg[n], a = vh[n];
      a2 += (o - f) * (o - f) * a;
      b2 += f * f * a;
      ai = fmax(ai, fabs(o - f));
      bi = fmax(bi, fabs(f));
      el[vlev[b]] = fmax(el[vlev[b]], fabs(o - f));
      int bd = 0;
      for (int k = 0; k < 9; k++)
        bd |= vrel[b * 9 + k] != 1 && k != 4;
      if (bd)
        eg = fmax(eg, fabs(o - f));
    }
    for (int l = 0; l <= Lm; l++) {
      int nbl = 0;
      for (int c = 0; c < N / NB; c++)
        nbl += vlev[c] == l;
      fprintf(stderr, "vm2: level %d blocks %d Linf %.3e\n", l, nbl, el[l] / bi);
    }
    fprintf(stderr, "vm2: Linf in blocks next to a level change %.3e\n", eg / bi);
    printf("t L2 Linf nodes_mean regrids\n% .16e % .16e % .16e % .16e %ld\n", nst * H0, sqrt(a2 / b2), ai / bi,
           nst > 0 ? nodet / nst : (double)N, nreg);
  }
  return 0;
}
