#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <thrust/execution_policy.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>

#ifndef M2LF
#define M2LF 3
#endif
#ifndef M2LT
#define M2LT 256
#endif
#ifndef NW
#define NW 64
#endif
enum { nthread = 256, PM = 40, DM = 20, STK = 4096, nwarp = NW, TB = 128 };
static const double pi = 3.141592653589793238;
__constant__ int clo[DM + 2], chi[DM + 2];
struct T {
  long N;
  int p, D, s, nc, nt, nr, per, n1, n2, n3, K;
  double delta, x0, y0, z0, S;
  int ord;
  double *x, *y, *z, *q, *cx, *cy, *cz, *phi, *gr;
  long *key, *kb;
  int *lev, *bx, *by, *bz, *par, *ch, *start, *count, *bid, *err;
  double2 *M, *L, *Gv, *Gm, *Gl;
  double *fa, *fb, *fc, *fn, *fs, *fr;
  int *uo, *ul, *vo, *vl, *wo, *wl, *xo, *xl;
};
static void ck(cudaError_t e, const char *what) {
  if (e != cudaSuccess) {
    fprintf(stderr, "fmm3: error: %s: %s\n", what, cudaGetErrorString(e));
    exit(1);
  }
}
static void launched(const char *what) { ck(cudaGetLastError(), what); }
static unsigned grid(long n) { return (n + nthread - 1) / nthread; }
__host__ __device__ static double2 mk(double a, double b) {
  double2 r;
  r.x = a;
  r.y = b;
  return r;
}
__host__ __device__ static double2 add(double2 a, double2 b) {
  return mk(a.x + b.x, a.y + b.y);
}
__host__ __device__ static double2 mul(double2 a, double2 b) {
  return mk(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
__host__ __device__ static double2 cmul(double2 a, double2 b) {
  return mk(a.x * b.x + a.y * b.y, a.x * b.y - a.y * b.x);
}
__host__ __device__ static double2 scl(double2 a, double s) {
  return mk(a.x * s, a.y * s);
}
__host__ __device__ static double2 cj(double2 a) { return mk(a.x, -a.y); }
__host__ __device__ static int ix(int n, int m) { return n * n + n + m; }
__host__ __device__ static long goff(int n) {
  return (long)n * (2 * n - 1) * (2 * n + 1) / 3;
}
__host__ __device__ static long morton(long ix_, long iy, long iz, int D) {
  long k = 0;
  for (int b = 0; b < D; b++) {
    k |= ((ix_ >> b) & 1L) << (3 * b + 2);
    k |= ((iy >> b) & 1L) << (3 * b + 1);
    k |= ((iz >> b) & 1L) << (3 * b);
  }
  return k;
}
__global__ static void keys(T t) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  long n = 1L << t.D;
  if (i >= t.N)
    return;
  long a = (long)((t.x[i] - t.x0) / t.S * (double)n);
  long b = (long)((t.y[i] - t.y0) / t.S * (double)n);
  long c = (long)((t.z[i] - t.z0) / t.S * (double)n);
  a = a < t.n1 * n - 1 ? a : t.n1 * n - 1;
  b = b < t.n2 * n - 1 ? b : t.n2 * n - 1;
  c = c < t.n3 * n - 1 ? c : t.n3 * n - 1;
  long r = ((a >> t.D) * t.n2 + (b >> t.D)) * t.n3 + (c >> t.D);
  t.key[i] = (r << (3 * t.D)) | morton(a & (n - 1), b & (n - 1), c & (n - 1), t.D);
}
__global__ static void gather(long N, const long *perm, const double *a,
                              double *b, int st, int off) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i < N)
    b[i] = a[perm[i] * st + off];
}
__device__ static long lower(const long *k, long n, long v) {
  long lo = 0, hi = n;
  while (lo < hi) {
    long m = (lo + hi) / 2;
    if (k[m] < v)
      lo = m + 1;
    else
      hi = m;
  }
  return lo;
}
__global__ static void counts(T t, int lo, int hi, int l, int *split) {
  int b = lo + blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= hi)
    return;
  int sh = 3 * (t.D - l);
  long msk = (1L << l) - 1;
  long kb = ((((long)t.bx[b] >> l) * t.n2 + ((long)t.by[b] >> l)) * t.n3 + ((long)t.bz[b] >> l)) << (3 * l) | morton(t.bx[b] & msk, t.by[b] & msk, t.bz[b] & msk, l);
  long k0 = kb << sh;
  long a0 = lower(t.key, t.N, k0);
  long a1 = lower(t.key, t.N, k0 + (1L << sh));
  double h = t.S / (double)(1L << l);
  t.start[b] = (int)a0;
  t.count[b] = (int)(a1 - a0);
  t.lev[b] = l;
  t.kb[b] = kb;
  t.cx[b] = t.x0 + h * (t.bx[b] + 0.5);
  t.cy[b] = t.y0 + h * (t.by[b] + 0.5);
  t.cz[b] = t.z0 + h * (t.bz[b] + 0.5);
  t.ch[b] = -1;
  split[b - lo] = (a1 - a0 > t.s && l < t.D) ? 1 : 0;
}
__global__ static void children(T t, int lo, int hi, const int *split,
                                const int *sc) {
  int b = lo + blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= hi || !split[b - lo])
    return;
  int c = hi + 8 * sc[b - lo];
  t.ch[b] = c;
  for (int q = 0; q < 8; q++) {
    t.bx[c + q] = 2 * t.bx[b] + ((q >> 2) & 1);
    t.by[c + q] = 2 * t.by[b] + ((q >> 1) & 1);
    t.bz[c + q] = 2 * t.bz[b] + (q & 1);
    t.par[c + q] = b;
  }
}
__device__ static int find(T t, int l, long a, long b, long c) {
  long n1 = (long)t.n1 << l, n2 = (long)t.n2 << l, n3 = (long)t.n3 << l, msk = (1L << l) - 1;
  int sx = 0, sy = 0, sz = 0;
  if (t.per) {
    sx = a >= 0 ? a / n1 : -((-a + n1 - 1) / n1);
    sy = b >= 0 ? b / n2 : -((-b + n2 - 1) / n2);
    sz = c >= 0 ? c / n3 : -((-c + n3 - 1) / n3);
    a -= sx * n1;
    b -= sy * n2;
    c -= sz * n3;
    if (sx < -2 || sx > 2 || sy < -2 || sy > 2 || sz < -2 || sz > 2)
      return -1;
  } else if (a < 0 || b < 0 || c < 0 || a >= n1 || b >= n2 || c >= n3)
    return -1;
  long k = ((((a >> l) * t.n2 + (b >> l)) * t.n3 + (c >> l)) << (3 * l)) | morton(a & msk, b & msk, c & msk, l);
  int lo = clo[l], hi = chi[l];
  while (lo < hi) {
    int m = (lo + hi) / 2;
    if (t.kb[m] < k)
      lo = m + 1;
    else
      hi = m;
  }
  if (lo >= chi[l] || t.kb[lo] != k)
    return -1;
  return lo * 125 + (sx + 2) * 25 + (sy + 2) * 5 + (sz + 2);
}
__device__ static int adj(T t, int ca, int b) {
  int a = ca / 125, sc = ca % 125;
  long sa = t.D - t.lev[a], sb = t.D - t.lev[b];
  long A[3] = {t.bx[a], t.by[a], t.bz[a]}, B[3] = {t.bx[b], t.by[b], t.bz[b]};
  long sh[3] = {(long)(sc / 25 - 2) * t.n1, (long)((sc / 5) % 5 - 2) * t.n2, (long)(sc % 5 - 2) * t.n3};
  for (int k = 0; k < 3; k++) {
    long a0 = (A[k] << sa) + (sh[k] << t.D), a1 = ((A[k] + 1) << sa) + (sh[k] << t.D);
    if (a0 > ((B[k] + 1) << sb) || (B[k] << sb) > a1)
      return 0;
  }
  return 1;
}
__device__ static void put(int *out, int *n, int v) {
  if (out != NULL)
    out[*n] = v;
  (*n)++;
}
__device__ static void descend(T t, int b, int *uo, int *nu, int *wo, int *nw) {
  int st[STK], ns = 0, l = t.lev[b];
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      for (int k = -1; k <= 1; k++) {
        int c = find(t, l, (long)t.bx[b] + i, (long)t.by[b] + j, (long)t.bz[b] + k);
        if (c >= 0)
          st[ns++] = c;
      }
  while (ns > 0) {
    int c = st[--ns], cb = c / 125, sc = c % 125;
    if (t.ch[cb] < 0) {
      put(uo, nu, c);
      continue;
    }
    for (int q = 0; q < 8; q++) {
      int d = (t.ch[cb] + q) * 125 + sc;
      if (adj(t, d, b)) {
        if (ns >= STK) {
          atomicExch(t.err, 1);
          return;
        }
        st[ns++] = d;
      } else
        put(wo, nw, d);
    }
  }
}
__device__ static void ucoarse(T t, int b, int *out, int *n) {
  int lb = t.lev[b];
  for (int l = 0; l < lb; l++) {
    long a = (long)t.bx[b] >> (lb - l), bb = (long)t.by[b] >> (lb - l), c = (long)t.bz[b] >> (lb - l);
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
        for (int k = -1; k <= 1; k++) {
          int d = find(t, l, a + i, bb + j, c + k);
          if (d >= 0 && t.ch[d / 125] < 0 && adj(t, d, b))
            put(out, n, d);
        }
  }
}
__device__ static void vlist(T t, int b, int *out, int *n) {
  int l = t.lev[b];
  if (l == 0 && t.per && t.K > 1) {
    for (int i = -t.K; i <= t.K; i++)
      for (int j = -t.K; j <= t.K; j++)
        for (int k = -t.K; k <= t.K; k++) {
          if (abs(i) <= 1 && abs(j) <= 1 && abs(k) <= 1)
            continue;
          int c = find(t, 0, (long)t.bx[b] + i, (long)t.by[b] + j, (long)t.bz[b] + k);
          if (c >= 0)
            put(out, n, c);
        }
    return;
  }
  if (l < (t.per ? 1 : 2))
    return;
  int P = t.par[b];
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      for (int k = -1; k <= 1; k++) {
        int c = find(t, l - 1, (long)t.bx[P] + i, (long)t.by[P] + j, (long)t.bz[P] + k);
        if (c < 0 || t.ch[c / 125] < 0)
          continue;
        for (int q = 0; q < 8; q++) {
          int d = (t.ch[c / 125] + q) * 125 + c % 125;
          if (!adj(t, d, b))
            put(out, n, d);
        }
      }
}
__device__ static void xlist(T t, int b, int *out, int *n) {
  int lb = t.lev[b];
  if (lb == 0)
    return;
  int P = t.par[b];
  for (int l = 0; l < lb; l++) {
    long a = (long)t.bx[b] >> (lb - l), bb = (long)t.by[b] >> (lb - l), c = (long)t.bz[b] >> (lb - l);
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
        for (int k = -1; k <= 1; k++) {
          int d = find(t, l, a + i, bb + j, c + k);
          if (d >= 0 && t.ch[d / 125] < 0 && adj(t, d, P) && !adj(t, d, b))
            put(out, n, d);
        }
  }
}
__device__ static void shv(T t, int code, double *sx, double *sy, double *sz) {
  int sc = code % 125;
  *sx = (sc / 25 - 2) * t.n1 * t.S;
  *sy = ((sc / 5) % 5 - 2) * t.n2 * t.S;
  *sz = (sc % 5 - 2) * t.n3 * t.S;
}
__device__ static void offs(T t, int code, int b, int *ox, int *oy, int *oz) {
  int s = code / 125, sc = code % 125, l = t.lev[b];
  *ox = t.bx[s] + (sc / 25 - 2) * (t.n1 << l) - t.bx[b];
  *oy = t.by[s] + ((sc / 5) % 5 - 2) * (t.n2 << l) - t.by[b];
  *oz = t.bz[s] + (sc % 5 - 2) * (t.n3 << l) - t.bz[b];
}
__global__ static void lists(T t, int nb, int fill) {
  int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= nb)
    return;
  int nu = 0, nv = 0, nw = 0, nx = 0;
  if (t.ch[b] < 0) {
    ucoarse(t, b, fill ? t.ul + t.uo[b] : NULL, &nu);
    descend(t, b, fill ? t.ul + t.uo[b] : NULL, &nu,
            fill ? t.wl + t.wo[b] : NULL, &nw);
  }
  vlist(t, b, fill ? t.vl + t.vo[b] : NULL, &nv);
  xlist(t, b, fill ? t.xl + t.xo[b] : NULL, &nx);
  if (!fill) {
    t.uo[b] = nu;
    t.vo[b] = nv;
    t.wo[b] = nw;
    t.xo[b] = nx;
  }
}
__global__ static void setbid(T t, int nb) {
  int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b >= nb || t.ch[b] >= 0)
    return;
  for (int i = t.start[b]; i < t.start[b] + t.count[b]; i++)
    t.bid[i] = b;
}
__global__ static void p2m(T t, int nb) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int P = t.p + 1;
  if (j >= (long)nb * t.nc * P)
    return;
  int m = j % P, c = (j / P) % t.nc, b = j / P / t.nc;
  if (t.ch[b] >= 0)
    return;
  double h = t.S / (double)(1L << t.lev[b]);
  double2 acc[PM + 1];
  for (int n = m; n <= t.p; n++)
    acc[n] = mk(0, 0);
  for (int i = t.start[b]; i < t.start[b] + t.count[b]; i++) {
    double X = (t.x[i] - t.cx[b]) / h, Y = (t.y[i] - t.cy[b]) / h, Z = (t.z[i] - t.cz[b]) / h,
           r2 = X * X + Y * Y + Z * Z, qq = t.q[(long)i * t.nc + c];
    double2 w = mk(X, Y), u = mk(1, 0), u1, u2;
    for (int k = 1; k <= m; k++)
      u = scl(mul(u, w), sqrt((2.0 * k - 1) / (2.0 * k)));
    acc[m] = add(acc[m], scl(u, qq));
    u2 = mk(0, 0);
    u1 = u;
    for (int n = m + 1; n <= t.p; n++) {
      double2 v = scl(add(scl(u1, (2 * n - 1) * Z), scl(u2, -sqrt((double)(n - 1) * (n - 1) - (double)m * m) * r2)),
                      1 / sqrt((double)n * n - (double)m * m));
      acc[n] = add(acc[n], scl(v, qq));
      u2 = u1;
      u1 = v;
    }
  }
  double2 *Mb = t.M + ((long)b * t.nc + c) * t.nt;
  for (int n = m; n <= t.p; n++) {
    Mb[ix(n, m)] = acc[n];
    if (m > 0)
      Mb[ix(n, -m)] = scl(cj(acc[n]), m % 2 ? -1 : 1);
  }
}
__global__ static void p2l(T t, int nb) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int P = t.p + 1;
  if (j >= (long)nb * t.nc * P)
    return;
  int m = j % P, c = (j / P) % t.nc, b = j / P / t.nc;
  if (t.xo[b] == t.xo[b + 1])
    return;
  double h = t.S / (double)(1L << t.lev[b]);
  double2 acc[PM + 1];
  for (int n = m; n <= t.p; n++)
    acc[n] = mk(0, 0);
  for (int e = t.xo[b]; e < t.xo[b + 1]; e++) {
    int cl = t.xl[e] / 125;
    double ssx, ssy, ssz;
    shv(t, t.xl[e], &ssx, &ssy, &ssz);
    for (int i = t.start[cl]; i < t.start[cl] + t.count[cl]; i++) {
      double X = (t.cx[b] - (t.x[i] + ssx)) / h, Y = (t.cy[b] - (t.y[i] + ssy)) / h, Z = (t.cz[b] - (t.z[i] + ssz)) / h,
             r2 = X * X + Y * Y + Z * Z, qq = t.q[(long)i * t.nc + c];
      double2 w = scl(mk(X, Y), 1 / r2), u = mk(1 / sqrt(r2), 0), u1, u2;
      for (int k = 1; k <= m; k++)
        u = scl(mul(u, w), sqrt((2.0 * k - 1) / (2.0 * k)));
      acc[m] = add(acc[m], scl(u, qq));
      u2 = mk(0, 0);
      u1 = u;
      for (int n = m + 1; n <= t.p; n++) {
        double2 v = scl(add(scl(u1, (2 * n - 1) * Z / r2), scl(u2, -sqrt((double)(n - 1) * (n - 1) - (double)m * m) / r2)),
                        1 / sqrt((double)n * n - (double)m * m));
        acc[n] = add(acc[n], scl(v, qq));
        u2 = u1;
        u1 = v;
      }
    }
  }
  double2 *Lb = t.L + ((long)b * t.nc + c) * t.nt;
  for (int n = m; n <= t.p; n++) {
    Lb[ix(n, m)] = add(Lb[ix(n, m)], acc[n]);
    if (m > 0)
      Lb[ix(n, -m)] = add(Lb[ix(n, -m)], scl(cj(acc[n]), m % 2 ? -1 : 1));
  }
}
__device__ static void hnm(int j, int *n, int *m) {
  int a = (int)((sqrt(8.0 * j + 1) - 1) / 2);
  while (a * (a + 1) / 2 > j)
    a--;
  while ((a + 1) * (a + 2) / 2 <= j)
    a++;
  *n = a;
  *m = j - a * (a + 1) / 2;
}
__device__ static void rfw(const double2 *G, const double2 *A, double2 *B, int nh) {
  for (int j = threadIdx.x; j < nh; j += blockDim.x) {
    int n, m;
    hnm(j, &n, &m);
    int w = 2 * n + 1;
    const double2 *g = G + goff(n) + (m + n) * w + n;
    double2 v = mul(g[0], A[n * n + n]);
    for (int l = 1; l <= n; l++) {
      double2 a = A[n * n + n + l], na = scl(cj(a), l % 2 ? -1 : 1);
      v = add(v, add(mul(g[l], a), mul(g[-l], na)));
    }
    B[n * n + n + m] = v;
  }
}
__device__ static void rbk(const double2 *G, const double2 *A, double2 *C, int nh, int acc) {
  for (int j = threadIdx.x; j < nh; j += blockDim.x) {
    int n, l;
    hnm(j, &n, &l);
    int w = 2 * n + 1;
    const double2 *g = G + goff(n) + l + n;
    double2 v = cmul(g[n * w], A[n * n + n]);
    for (int m = 1; m <= n; m++) {
      double2 a = A[n * n + n + m], na = scl(cj(a), m % 2 ? -1 : 1);
      v = add(v, add(cmul(g[(m + n) * w], a), cmul(g[(n - m) * w], na)));
    }
    C[j] = acc ? add(C[j], v) : v;
  }
}
__device__ static void flush(T t, double2 *dst, const double2 *C, int nh, int acc) {
  for (int j = threadIdx.x; j < nh; j += blockDim.x) {
    int n, m;
    hnm(j, &n, &m);
    double2 v = C[j];
    dst[n * n + n + m] = acc ? add(dst[n * n + n + m], v) : v;
    if (m > 0)
      dst[n * n + n - m] = acc ? add(dst[n * n + n - m], scl(cj(v), m % 2 ? -1 : 1)) : scl(cj(v), m % 2 ? -1 : 1);
  }
}
__global__ static void m2m(T t, int lo, int hi, int l) {
  extern __shared__ double2 sh[];
  int b = lo + blockIdx.x, nt = t.nt, nh = (t.p + 1) * (t.p + 2) / 2, P1 = t.p + 1;
  if (b >= hi || t.ch[b] < 0)
    return;
  double2 *A = sh, *B = sh + nt, *C = sh + 2 * nt;
  double r = 0.25 * sqrt(3.0);
  for (int c = 0; c < t.nc; c++) {
    for (int j = threadIdx.x; j < nh; j += blockDim.x)
      C[j] = mk(0, 0);
    for (int q = 0; q < 8; q++) {
      const double2 *Mc = t.M + ((long)(t.ch[b] + q) * t.nc + c) * nt;
      const double2 *G = t.Gm + q * goff(t.p + 1);
      __syncthreads();
      for (int j = threadIdx.x; j < nt; j += blockDim.x)
        A[j] = Mc[j];
      __syncthreads();
      rfw(G, A, B, nh);
      __syncthreads();
      for (int j = threadIdx.x; j < nh; j += blockDim.x) {
        int n, m;
        hnm(j, &n, &m);
        double2 v = mk(0, 0);
        const double *f = t.fa + ((long)m * P1 + n) * P1;
        double cr = 1, ch = 1;
        for (int k = m; k < n; k++)
          cr *= r;
        for (int k = 0; k < m; k++)
          ch *= 0.5;
        for (int k = m; k <= n; k++) {
          v = add(v, scl(B[k * k + k + m], f[k] * cr * ch));
          cr /= r;
          ch *= 0.5;
        }
        A[n * n + n + m] = v;
      }
      __syncthreads();
      rbk(G, A, C, nh, 1);
    }
    __syncthreads();
    flush(t, t.M + ((long)b * t.nc + c) * nt, C, nh, 0);
    __syncthreads();
  }
}
__global__ static void m2lh(T t, int nb) {
  extern __shared__ double2 sh[];
  int b = blockIdx.x, nt = t.nt, nh = (t.p + 1) * (t.p + 2) / 2, P1 = t.p + 1;
  if (b >= nb || t.vo[b] == t.vo[b + 1])
    return;
  double2 *A = sh, *B = sh + nt, *C = sh + 2 * nt;
  for (int c = 0; c < t.nc; c++) {
    for (int j = threadIdx.x; j < nh; j += blockDim.x)
      C[j] = mk(0, 0);
    for (int e = t.vo[b]; e < t.vo[b + 1]; e++) {
      int s = t.vl[e] / 125, ox, oy, oz;
      offs(t, t.vl[e], b, &ox, &oy, &oz);
      int o = (ox + 3) * 49 + (oy + 3) * 7 + (oz + 3);
      const double2 *G = t.Gv + o * goff(t.p + 1);
      const double2 *Ms = t.M + ((long)s * t.nc + c) * nt;
      double ir = 1 / sqrt((double)(ox * ox + oy * oy + oz * oz));
      __syncthreads();
      for (int j = threadIdx.x; j < nt; j += blockDim.x)
        A[j] = Ms[j];
      __syncthreads();
      rfw(G, A, B, nh);
      __syncthreads();
      for (int j = threadIdx.x; j < nh; j += blockDim.x) {
        int n, m;
        hnm(j, &n, &m);
        double cc = pow(ir, n + m + 1);
        double2 v = mk(0, 0);
#if M2LF == 1
        for (int k = m; k <= t.p; k++) {
          v = add(v, scl(B[k * k + k + m], t.fn[((long)n * P1 + k) * P1 + m] * cc));
          cc *= ir;
        }
#elif M2LF == 3
        double f = t.fs[n * P1 + m] * cc;
        for (int k = m; k <= t.p; k++) {
          v = add(v, scl(B[k * k + k + m], f));
          f *= (n + k + 1) * ir * t.fr[k * P1 + m];
        }
#elif M2LF == 2
        double f = t.fs[n * P1 + m] * cc;
        for (int k = m; k <= t.p; k++) {
          v = add(v, scl(B[k * k + k + m], f));
          f *= (n + k + 1) * ir / sqrt((double)(k - m + 1) * (k + m + 1));
        }
#else
        const double *f = t.fc + ((long)m * P1 + n) * P1;
        for (int k = m; k <= t.p; k++) {
          v = add(v, scl(B[k * k + k + m], f[k] * cc));
          cc *= ir;
        }
#endif
        A[n * n + n + m] = scl(v, m % 2 ? -1 : 1);
      }
      __syncthreads();
      rbk(G, A, C, nh, 1);
    }
    __syncthreads();
    flush(t, t.L + ((long)b * t.nc + c) * nt, C, nh, 1);
    __syncthreads();
  }
}
__global__ static void l2l(T t, int lo, int hi, int l) {
  extern __shared__ double2 sh[];
  int b = lo + blockIdx.x, nt = t.nt, nh = (t.p + 1) * (t.p + 2) / 2, P1 = t.p + 1;
  if (b >= hi)
    return;
  double2 *A = sh, *B = sh + nt, *C = sh + 2 * nt;
  int q = ((t.bx[b] & 1) << 2) | ((t.by[b] & 1) << 1) | (t.bz[b] & 1);
  const double2 *G = t.Gl + q * goff(t.p + 1);
  double r = 0.25 * sqrt(3.0);
  for (int c = 0; c < t.nc; c++) {
    const double2 *Fp = t.L + ((long)t.par[b] * t.nc + c) * nt;
    __syncthreads();
    for (int j = threadIdx.x; j < nt; j += blockDim.x)
      A[j] = Fp[j];
    __syncthreads();
    rfw(G, A, B, nh);
    __syncthreads();
    for (int j = threadIdx.x; j < nh; j += blockDim.x) {
      int n, m;
      hnm(j, &n, &m);
      double2 v = mk(0, 0);
      const double *f = t.fb + ((long)m * P1 + n) * P1;
      double cr = 1, ch = 0.5;
      for (int k = 0; k < n; k++)
        ch *= 0.5;
      for (int k = n; k <= t.p; k++) {
        v = add(v, scl(B[k * k + k + m], f[k] * cr * ch));
        cr *= r;
      }
      A[n * n + n + m] = v;
    }
    __syncthreads();
    rbk(G, A, C, nh, 0);
    __syncthreads();
    flush(t, t.L + ((long)b * t.nc + c) * nt, C, nh, 1);
    __syncthreads();
  }
}
__global__ static void dense(int p, const double2 *Ih, const double *fe, double2 *Tm) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int nt = (p + 1) * (p + 1), nh = (p + 1) * (p + 2) / 2;
  if (j >= (long)nh * nt)
    return;
  int r = j % nh, cidx = j / nh, n, m;
  hnm(r, &n, &m);
  int k = (int)sqrt((double)cidx);
  while (k * k > cidx)
    k--;
  while ((k + 1) * (k + 1) <= cidx)
    k++;
  int l = cidx - k * k - k, q = m - l, J = n + k;
  Tm[j] = scl(Ih[J * J + J + q], fe[j] * ((l < 0 ? -l : l) % 2 ? -1 : 1));
}
__global__ static void gath(T t, const int *src, int K, int c, double2 *X) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (j >= (long)K * t.nt)
    return;
  X[j] = t.M[((long)src[j / t.nt] * t.nc + c) * t.nt + j % t.nt];
}
__global__ static void scat(T t, const int *tg, int K, int c, const double2 *Y) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int nh = (t.p + 1) * (t.p + 2) / 2, n, m;
  if (j >= (long)K * nh)
    return;
  int col = j / nh, r = j % nh;
  hnm(r, &n, &m);
  double2 v = Y[j], *Lb = t.L + ((long)tg[col] * t.nc + c) * t.nt;
  Lb[n * n + n + m] = add(Lb[n * n + n + m], v);
  if (m > 0)
    Lb[n * n + n - m] = add(Lb[n * n + n - m], scl(cj(v), m % 2 ? -1 : 1));
}
__global__ static void farop(T t, int nr, const double2 *Tf) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int nh = (t.p + 1) * (t.p + 2) / 2, n, m;
  if (j >= (long)nr * t.nc * nh)
    return;
  int i = j % nh, c = (j / nh) % t.nc, b = j / nh / t.nc;
  hnm(i, &n, &m);
  double2 v = mk(0, 0);
  for (int s = 0; s < nr; s++) {
    int dx = t.bx[s] - t.bx[b] + t.n1 - 1, dy = t.by[s] - t.by[b] + t.n2 - 1, dz = t.bz[s] - t.bz[b] + t.n3 - 1;
    const double2 *Tm = Tf + (((long)dx * (2 * t.n2 - 1) + dy) * (2 * t.n3 - 1) + dz) * nh * t.nt + (long)i * t.nt;
    const double2 *Ms = t.M + ((long)s * t.nc + c) * t.nt;
    for (int k = 0; k < t.nt; k++)
      v = add(v, mul(Tm[k], Ms[k]));
  }
  double2 *Lb = t.L + ((long)b * t.nc + c) * t.nt;
  Lb[n * n + n + m] = add(Lb[n * n + n + m], v);
  if (m > 0)
    Lb[n * n + n - m] = add(Lb[n * n + n - m], scl(cj(v), m % 2 ? -1 : 1));
}
__device__ static double2 Fat(const double2 *F, int p, int n, int m) {
  if (n > p || m > n || m < -n)
    return mk(0, 0);
  return F[ix(n, m)];
}
__global__ static void l2p(T t) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i >= t.N)
    return;
  int b = t.bid[i], p = t.p;
  double h = t.S / (double)(1L << t.lev[b]);
  double X = (t.cx[b] - t.x[i]) / h, Y = (t.cy[b] - t.y[i]) / h, Z = (t.cz[b] - t.z[i]) / h,
         r2 = X * X + Y * Y + Z * Z;
  double2 w = mk(X, Y);
  for (int c = 0; c < t.nc; c++) {
    const double2 *F = t.L + ((long)b * t.nc + c) * t.nt;
    double ph = 0, gx = 0, gy = 0, gz = 0;
    double2 um = mk(1, 0);
    for (int m = 0; m <= p; m++) {
      if (m > 0)
        um = scl(mul(um, w), sqrt((2.0 * m - 1) / (2.0 * m)));
      double2 u1 = um, u2 = mk(0, 0);
      for (int n = m; n <= p; n++) {
        double2 u;
        if (n == m)
          u = um;
        else {
          u = scl(add(scl(u1, (2 * n - 1) * Z), scl(u2, -sqrt((double)(n - 1) * (n - 1) - (double)m * m) * r2)),
                  1 / sqrt((double)n * n - (double)m * m));
          u2 = u1;
          u1 = u;
        }
        for (int sg = 0; sg < (m > 0 ? 2 : 1); sg++) {
          int mm = sg ? -m : m;
          double2 cu = cj(sg ? scl(cj(u), m % 2 ? -1 : 1) : u);
          ph += mul(cu, F[ix(n, mm)]).x;
          if (n < p) {
            double ap = sqrt((double)(n + mm + 1) * (n + mm + 2)), am = sqrt((double)(n - mm + 1) * (n - mm + 2));
            double2 fp = Fat(F, p, n + 1, mm + 1), fm = Fat(F, p, n + 1, mm - 1), fz = Fat(F, p, n + 1, mm);
            gx -= 0.5 * mul(cu, mk(ap * fp.x - am * fm.x, ap * fp.y - am * fm.y)).x;
            gy -= 0.5 * mul(mul(mk(0, -1), cu), mk(ap * fp.x + am * fm.x, ap * fp.y + am * fm.y)).x;
            gz -= sqrt((double)(n + 1) * (n + 1) - (double)mm * mm) * mul(cu, fz).x;
          }
        }
      }
    }
    t.phi[i * t.nc + c] = ph / h;
    t.gr[(i * t.nc + c) * 3] = gx / (h * h);
    t.gr[(i * t.nc + c) * 3 + 1] = gy / (h * h);
    t.gr[(i * t.nc + c) * 3 + 2] = gz / (h * h);
  }
}
__device__ static double2 Mat(const double2 *M, int p, int n, int m) {
  if (n < 0 || n > p || m > n || m < -n)
    return mk(0, 0);
  return M[ix(n, m)];
}
__global__ static void m2p(T t) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i >= t.N)
    return;
  int b = t.bid[i], p = t.p;
  if (t.wo[b] == t.wo[b + 1])
    return;
  for (int c = 0; c < t.nc; c++) {
    double ph = 0, gx = 0, gy = 0, gz = 0;
    for (int e = t.wo[b]; e < t.wo[b + 1]; e++) {
      int d = t.wl[e] / 125;
      double ssx, ssy, ssz, h = t.S / (double)(1L << t.lev[d]), sp = 0, sx = 0, sy = 0, sz = 0;
      shv(t, t.wl[e], &ssx, &ssy, &ssz);
      const double2 *M = t.M + ((long)d * t.nc + c) * t.nt;
      double X = (t.x[i] - (t.cx[d] + ssx)) / h, Y = (t.y[i] - (t.cy[d] + ssy)) / h, Z = (t.z[i] - (t.cz[d] + ssz)) / h,
             r2 = X * X + Y * Y + Z * Z;
      double2 w = scl(mk(X, Y), 1 / r2), um = mk(1 / sqrt(r2), 0);
      for (int m = 0; m <= p + 1; m++) {
        if (m > 0)
          um = scl(mul(um, w), sqrt((2.0 * m - 1) / (2.0 * m)));
        double2 u1 = um, u2 = mk(0, 0);
        for (int n = m; n <= p + 1; n++) {
          double2 u;
          if (n == m)
            u = um;
          else {
            u = scl(add(scl(u1, (2 * n - 1) * Z / r2), scl(u2, -sqrt((double)(n - 1) * (n - 1) - (double)m * m) / r2)),
                    1 / sqrt((double)n * n - (double)m * m));
            u2 = u1;
            u1 = u;
          }
          for (int sg = 0; sg < (m > 0 ? 2 : 1); sg++) {
            int mm = sg ? -m : m;
            double2 cu = cj(sg ? scl(cj(u), m % 2 ? -1 : 1) : u);
            if (n <= p)
              sp += mul(cu, M[ix(n, mm)]).x;
            if (n >= 1) {
              double2 mu = Mat(M, p, n - 1, mm + 1), md = Mat(M, p, n - 1, mm - 1), mz = Mat(M, p, n - 1, mm);
              double au = n - mm - 1 >= 0 ? sqrt((double)(n - mm - 1) * (n - mm)) : 0, ad = n + mm - 1 >= 0 ? sqrt((double)(n + mm - 1) * (n + mm)) : 0;
              sx += 0.5 * mul(cu, mk(au * mu.x - ad * md.x, au * mu.y - ad * md.y)).x;
              sy += 0.5 * mul(mul(mk(0, -1), cu), mk(au * mu.x + ad * md.x, au * mu.y + ad * md.y)).x;
              sz -= sqrt((double)n * n - (double)mm * mm) * mul(cu, mz).x;
            }
          }
        }
      }
      ph += sp / h;
      gx += sx / (h * h);
      gy += sy / (h * h);
      gz += sz / (h * h);
    }
    t.phi[i * t.nc + c] += ph;
    t.gr[(i * t.nc + c) * 3] += gx;
    t.gr[(i * t.nc + c) * 3 + 1] += gy;
    t.gr[(i * t.nc + c) * 3 + 2] += gz;
  }
}
__global__ static void quad(T t, const double *qd) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i >= t.N)
    return;
  double V = (double)t.n1 * t.n2 * t.n3 * t.S * t.S * t.S, x = t.x[i], y = t.y[i], z = t.z[i];
  for (int c = 0; c < t.nc; c++) {
    const double *a = qd + 5 * c;
    t.phi[i * t.nc + c] += 2 * pi / (3 * V) * (a[0] * (x * x + y * y + z * z) - 2 * (x * a[1] + y * a[2] + z * a[3]) + a[4]);
    t.gr[(i * t.nc + c) * 3] += 4 * pi / (3 * V) * (a[0] * x - a[1]);
    t.gr[(i * t.nc + c) * 3 + 1] += 4 * pi / (3 * V) * (a[0] * y - a[2]);
    t.gr[(i * t.nc + c) * 3 + 2] += 4 * pi / (3 * V) * (a[0] * z - a[3]);
  }
}
__device__ static void nsum(double *s, double *c, double v) {
  double t = *s + v;
  if (fabs(*s) >= fabs(v))
    *c += (*s - t) + v;
  else
    *c += (v - t) + *s;
  *s = t;
}
__device__ static double core8(double r, double d, int ord) {
  if (d <= 0)
    return 1;
  double a = r / d, x = a * a;
  if (a < 0.3) {
    const double c8[12] = {35.0 / 4, -231.0 / 20, 429.0 / 56, -715.0 / 216, 1105.0 / 1056, -323.0 / 1248,
                           2261.0 / 43200, -437.0 / 48960, 115.0 / 87552, -115.0 / 677376, 29.0 / 1483776,
                           -899.0 / 443520000};
    const double c2[12] = {4.0 / 3, -4.0 / 5, 2.0 / 7, -2.0 / 27, 1.0 / 66, -1.0 / 390,
                           1.0 / 2700, -1.0 / 21420, 1.0 / 191520, -1.0 / 1905120, 1.0 / 20865600,
                           -1.0 / 249480000};
    const double *c = ord == 8 ? c8 : c2;
    double s = c[11];
    for (int k = 10; k >= 0; k--)
      s = s * x + c[k];
    return s * a * x / sqrt(pi);
  }
  if (ord == 8)
    return erf(a) - 2 / sqrt(pi) * a * exp(-x) * (1 - 89.0 / 24 * x + 5.0 / 3 * x * x - x * x * x / 6);
  return erf(a) - 2 / sqrt(pi) * a * exp(-x);
}
__device__ static int lmort(int x, int y, int z, int B) {
  int i = 0;
  for (int b = 0; (1 << b) < B; b++)
    i |= (((x >> b) & 1) << (3 * b + 2)) | (((y >> b) & 1) << (3 * b + 1)) | (((z >> b) & 1) << (3 * b));
  return i;
}
__global__ static void fpatch(int nbt, int k0, int nl, int B, const int *srca, const int *start, const double *q, double *P) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int M = 4 * B;
  long M3 = (long)M * M * M;
  if (j >= (long)nbt * 3 * M3)
    return;
  int bb = j / (3 * M3), c = (j / M3) % 3;
  long r = j % M3;
  int px = r / (M * M), py = (r / M) % M, pz = r % M;
  double v = 0;
  if (px < 3 * B && py < 3 * B && pz < 3 * B) {
    int o = (px / B) * 9 + (py / B) * 3 + pz / B, s = srca[(long)o * nl + k0 + bb];
    if (s >= 0)
      v = q[((long)start[s] + lmort(px % B, py % B, pz % B, B)) * 3 + c];
  }
  P[j] = v;
}
__global__ static void fmulk(int nbt, long Mc3, const double2 *Gh, const double2 *A, double2 *U) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (j >= (long)nbt * Mc3)
    return;
  int bb = j / Mc3;
  long r = j % Mc3;
  const double2 *ax = A + ((long)bb * 3) * Mc3, *ay = ax + Mc3, *az = ay + Mc3;
  double2 gx = Gh[r], gy = Gh[Mc3 + r], gz = Gh[2 * Mc3 + r];
  double2 *ux = U + ((long)bb * 3) * Mc3, *uy = ux + Mc3, *uz = uy + Mc3;
  double2 X = ax[r], Y = ay[r], Z = az[r];
  ux[r] = mk(Y.x * gz.x - Y.y * gz.y - (Z.x * gy.x - Z.y * gy.y), Y.x * gz.y + Y.y * gz.x - (Z.x * gy.y + Z.y * gy.x));
  uy[r] = mk(Z.x * gx.x - Z.y * gx.y - (X.x * gz.x - X.y * gz.y), Z.x * gx.y + Z.y * gx.x - (X.x * gz.y + X.y * gz.x));
  uz[r] = mk(X.x * gy.x - X.y * gy.y - (Y.x * gx.x - Y.y * gx.y), X.x * gy.y + X.y * gy.x - (Y.x * gx.y + Y.y * gx.x));
}
__global__ static void fext(int nbt, int k0, int B, const int *lst, const double *R, double *ul) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int B3 = B * B * B, M = 4 * B;
  long M3 = (long)M * M * M;
  if (j >= (long)nbt * B3 * 3)
    return;
  int bb = j / (3 * B3), c = (j / B3) % 3, i = j % B3;
  int x = 0, y = 0, z = 0;
  for (int b = 0; (1 << b) < B; b++) {
    x |= ((i >> (3 * b + 2)) & 1) << b;
    y |= ((i >> (3 * b + 1)) & 1) << b;
    z |= ((i >> (3 * b)) & 1) << b;
  }
  double v = R[((long)bb * 3 + c) * M3 + ((long)(x + B) * M + y + B) * M + z + B] / (double)M3;
  ul[((long)lst[k0 + bb] + i) * 3 + c] += v;
}
__global__ static void lgath(int nl, int B3, const int *src, const int *start, const double *q, double *X) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (j >= (long)nl * B3 * 3)
    return;
  int i = j % B3, col = j / B3, b = col / 3, c = col % 3, s = src[b];
  X[j] = s < 0 ? 0 : q[((long)start[s] + i) * 3 + c];
}
__device__ static long bidx(T t, int l, long x, long y, long z) {
  long msk = (1L << l) - 1, r = ((x >> l) * t.n2 + (y >> l)) * t.n3 + (z >> l);
  return (r << (3 * l)) | morton(x & msk, y & msk, z & msk, l);
}
__device__ static double2 iget(const double2 *I, int J, int q) {
  if (q >= 0)
    return I[J * (J + 1) / 2 + q];
  return scl(cj(I[J * (J + 1) / 2 - q]), q % 2 ? -1 : 1);
}
__global__ static void mkA(int p, const double2 *I, const double *fe, double *A) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  int nh = (p + 1) * (p + 2) / 2, H = 16 * nh, nt = (p + 1) * (p + 1), nj = (2 * p + 1) * (2 * p + 2) / 2;
  if (j >= 26L * H * H)
    return;
  int R = j % H, C = (j / H) % H, d = j / H / H;
  d += d >= 13;
  int q = R / (2 * nh), rr = R % (2 * nh), qs = C / (2 * nh), cc = C % (2 * nh);
  int ox = 2 * (d / 9 - 1) + ((qs >> 2) & 1) - ((q >> 2) & 1), oy = 2 * ((d / 3) % 3 - 1) + ((qs >> 1) & 1) - ((q >> 1) & 1), oz = 2 * (d % 3 - 1) + (qs & 1) - (q & 1);
  double v = 0;
  if (abs(ox) > 1 || abs(oy) > 1 || abs(oz) > 1) {
    int n, m, k, l;
    hnm(rr / 2, &n, &m);
    hnm(cc / 2, &k, &l);
    const double2 *Io = I + (long)((ox + 3) * 49 + (oy + 3) * 7 + oz + 3) * nj;
    double2 T1 = scl(iget(Io, n + k, m - l), fe[(long)ix(n, m) * nt + ix(k, l)] * (l % 2 ? -1 : 1)), T2 = mk(0, 0);
    if (l > 0)
      T2 = scl(iget(Io, n + k, m + l), fe[(long)ix(n, m) * nt + ix(k, -l)]);
    if (rr % 2 == 0)
      v = cc % 2 == 0 ? T1.x + T2.x : T2.y - T1.y;
    else
      v = cc % 2 == 0 ? T1.y + T2.y : T1.x - T2.x;
  }
  A[j] = v;
}
__global__ static void mpack(T t, int l, int lo, int3 ax, int3 D, double *S) {
  int nh = (t.p + 1) * (t.p + 2) / 2;
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x, nc = (long)(D.x + 2) * (D.y + 2) * (D.z + 2);
  if (j >= nc * t.nc * 8 * nh)
    return;
  int h = j % nh, q = (j / nh) % 8, c = (j / nh / 8) % t.nc, a[3] = {ax.x, ax.y, ax.z};
  long cell = j / nh / 8 / t.nc, u[3], P[3];
  u[2] = cell % (D.z + 2);
  u[1] = (cell / (D.z + 2)) % (D.y + 2);
  u[0] = cell / ((long)(D.z + 2) * (D.y + 2));
  for (int k = 0; k < 3; k++)
    P[a[k]] = u[k] - 1;
  long NX = (long)t.n1 << l, NY = (long)t.n2 << l, NZ = (long)t.n3 << l;
  long x = 2 * P[0] + ((q >> 2) & 1), y = 2 * P[1] + ((q >> 1) & 1), z = 2 * P[2] + (q & 1);
  double2 v = mk(0, 0);
  if (t.per) {
    x = (x + NX) % NX;
    y = (y + NY) % NY;
    z = (z + NZ) % NZ;
  }
  if (x >= 0 && y >= 0 && z >= 0 && x < NX && y < NY && z < NZ) {
    int n, m;
    hnm(h, &n, &m);
    v = t.M[((lo + bidx(t, l, x, y, z)) * t.nc + c) * t.nt + n * n + n + m];
  }
  S[j * 2] = v.x;
  S[j * 2 + 1] = v.y;
}
__global__ static void munpack(T t, int l, int lo, int3 ax, int3 D, const double *S) {
  int nh = (t.p + 1) * (t.p + 2) / 2;
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (j >= (long)D.x * D.y * D.z * t.nc * 8 * nh)
    return;
  int h = j % nh, q = (j / nh) % 8, c = (j / nh / 8) % t.nc, n, m, a[3] = {ax.x, ax.y, ax.z};
  long cell = j / nh / 8 / t.nc, u[3], P[3];
  u[2] = cell % D.z;
  u[1] = (cell / D.z) % D.y;
  u[0] = cell / ((long)D.z * D.y);
  for (int k = 0; k < 3; k++)
    P[a[k]] = u[k];
  long cp = ((u[0] + 1) * (D.y + 2) + u[1] + 1) * (D.z + 2) + u[2] + 1;
  long b = lo + bidx(t, l, 2 * P[0] + ((q >> 2) & 1), 2 * P[1] + ((q >> 1) & 1), 2 * P[2] + (q & 1));
  const double *s = S + (((cp * t.nc + c) * 8 + q) * nh + h) * 2;
  double2 v = mk(s[0], s[1]), *Lb = t.L + (b * t.nc + c) * t.nt;
  hnm(h, &n, &m);
  Lb[n * n + n + m] = add(Lb[n * n + n + m], v);
  if (m > 0)
    Lb[n * n + n - m] = add(Lb[n * n + n - m], scl(cj(v), m % 2 ? -1 : 1));
}
__global__ static void lpack(T t, int l, int lo, int B3, int K, int3 ax, int3 D, double *X) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (j >= (long)(D.x + 2 * K) * (D.y + 2 * K) * (D.z + 2 * K) * 3 * B3)
    return;
  int i = j % B3, c = (j / B3) % 3, a[3] = {ax.x, ax.y, ax.z};
  long cell = j / B3 / 3, u[3], P[3];
  u[2] = cell % (D.z + 2 * K);
  u[1] = (cell / (D.z + 2 * K)) % (D.y + 2 * K);
  u[0] = cell / ((long)(D.z + 2 * K) * (D.y + 2 * K));
  for (int k = 0; k < 3; k++)
    P[a[k]] = u[k] - K;
  long NX = (long)t.n1 << l, NY = (long)t.n2 << l, NZ = (long)t.n3 << l, x = P[0], y = P[1], z = P[2];
  if (t.per) {
    x = (x % NX + NX) % NX;
    y = (y % NY + NY) % NY;
    z = (z % NZ + NZ) % NZ;
  }
  double v = 0;
  if (x >= 0 && y >= 0 && z >= 0 && x < NX && y < NY && z < NZ)
    v = t.q[((long)t.start[lo + bidx(t, l, x, y, z)] + i) * 3 + c];
  X[j] = v;
}
__global__ static void lcombp(T t, int l, int lo, int B3, int K, int3 ax, int3 D, const double *Y, double *ul) {
  long j = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (j >= (long)D.x * D.y * D.z * B3)
    return;
  int i = j % B3, a[3] = {ax.x, ax.y, ax.z};
  long cell = j / B3, u[3], P[3];
  u[2] = cell % D.z;
  u[1] = (cell / D.z) % D.y;
  u[0] = cell / ((long)D.z * D.y);
  for (int k = 0; k < 3; k++)
    P[a[k]] = u[k];
  long cp = ((u[0] + K) * (D.y + 2 * K) + u[1] + K) * (D.z + 2 * K) + u[2] + K;
  const double *y = Y + cp * 6 * B3 + i;
  long o = ((long)t.start[lo + bidx(t, l, P[0], P[1], P[2])] + i) * 3;
  ul[o] += y[5 * B3] - y[2 * B3];
  ul[o + 1] += y[B3] - y[4 * B3];
  ul[o + 2] += y[3 * B3] - y[0];
}
__global__ static void p2p(T t, int nb) {
  __shared__ double sx[nwarp], sy[nwarp], sz[nwarp], sq[3 * nwarp];
  int b = blockIdx.x;
  if (b >= nb || t.ch[b] >= 0)
    return;
  for (int t0 = 0; t0 < t.count[b]; t0 += blockDim.x) {
    int i = t.start[b] + t0 + threadIdx.x, ok = t0 + threadIdx.x < t.count[b];
    double xi = ok ? t.x[i] : 0, yi = ok ? t.y[i] : 0, zi = ok ? t.z[i] : 0;
    double ph[3] = {0, 0, 0}, g[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, cp[3] = {0, 0, 0}, cg[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int e = t.uo[b]; e < t.uo[b + 1]; e++) {
      int c = t.ul[e] / 125;
      double ssx, ssy, ssz;
      shv(t, t.ul[e], &ssx, &ssy, &ssz);
      for (int s0 = 0; s0 < t.count[c]; s0 += blockDim.x) {
        int js = t.start[c] + s0 + threadIdx.x;
        __syncthreads();
        if (s0 + threadIdx.x < t.count[c]) {
          sx[threadIdx.x] = t.x[js] + ssx;
          sy[threadIdx.x] = t.y[js] + ssy;
          sz[threadIdx.x] = t.z[js] + ssz;
          for (int k = 0; k < t.nc; k++)
            sq[threadIdx.x * 3 + k] = t.q[(long)js * t.nc + k];
        }
        __syncthreads();
        int m = t.count[c] - s0 < (int)blockDim.x ? t.count[c] - s0
                                                   : (int)blockDim.x;
        if (ok)
          for (int j = 0; j < m; j++) {
            double dx = xi - sx[j], dy = yi - sy[j], dz = zi - sz[j],
                   r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > 0) {
              double r = sqrt(r2), ir = 1 / r, g3 = core8(r, t.delta, t.ord) * ir * ir * ir;
              for (int k = 0; k < t.nc; k++) {
                double qq = sq[j * 3 + k];
#ifdef KP
                nsum(&ph[k], &cp[k], qq * ir);
                nsum(&g[3 * k], &cg[3 * k], -qq * g3 * dx);
                nsum(&g[3 * k + 1], &cg[3 * k + 1], -qq * g3 * dy);
                nsum(&g[3 * k + 2], &cg[3 * k + 2], -qq * g3 * dz);
#else
                ph[k] += qq * ir;
                g[3 * k] -= qq * g3 * dx;
                g[3 * k + 1] -= qq * g3 * dy;
                g[3 * k + 2] -= qq * g3 * dz;
#endif
              }
            }
          }
      }
    }
    if (ok)
      for (int k = 0; k < t.nc; k++) {
        t.phi[(long)i * t.nc + k] += ph[k] + cp[k];
        for (int d = 0; d < 3; d++)
          t.gr[((long)i * t.nc + k) * 3 + d] += g[3 * k + d] + cg[3 * k + d];
      }
  }
}
__global__ static void check(T t, const long *pos, double *ref) {
  __shared__ double red[TB][13];
  long i = pos[blockIdx.x];
  double xi = t.x[i], yi = t.y[i], zi = t.z[i];
  double ph[3] = {0, 0, 0}, g[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, cp[3] = {0, 0, 0}, cg[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (long j = threadIdx.x; j < t.N; j += blockDim.x) {
    double dx = xi - t.x[j], dy = yi - t.y[j], dz = zi - t.z[j],
           r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > 0) {
      double r = sqrt(r2), ir = 1 / r, g3 = core8(r, t.delta, t.ord) * ir * ir * ir;
      for (int k = 0; k < t.nc; k++) {
        double qq = t.q[j * t.nc + k];
#ifdef KC
        nsum(&ph[k], &cp[k], qq * ir);
        nsum(&g[3 * k], &cg[3 * k], -qq * g3 * dx);
        nsum(&g[3 * k + 1], &cg[3 * k + 1], -qq * g3 * dy);
        nsum(&g[3 * k + 2], &cg[3 * k + 2], -qq * g3 * dz);
#else
        ph[k] += qq * ir;
        g[3 * k] -= qq * g3 * dx;
        g[3 * k + 1] -= qq * g3 * dy;
        g[3 * k + 2] -= qq * g3 * dz;
#endif
      }
    }
  }
  for (int k = 0; k < 3; k++)
    red[threadIdx.x][k] = ph[k] + cp[k];
  for (int k = 0; k < 9; k++)
    red[threadIdx.x][3 + k] = g[k] + cg[k];
  __syncthreads();
  for (int h = blockDim.x / 2; h > 0; h /= 2) {
    if ((int)threadIdx.x < h)
      for (int k = 0; k < 12; k++)
        red[threadIdx.x][k] += red[threadIdx.x + h][k];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    for (int k = 0; k < 12; k++)
      ref[blockIdx.x * 12 + k] = red[0][k];
}
__global__ static void inverse(long N, const long *perm, long *pos) {
  long i = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (i < N)
    pos[perm[i]] = i;
}
static double *Dh[PM + 2];
static int ncache;
static int cn[1 << 14];
static double cb[1 << 14], *cr[1 << 14];
static void dbet(int n, double beta, double *d) {
  int w = 2 * n + 1;
  for (int i = 0; i < ncache; i++)
    if (cn[i] == n && cb[i] == beta) {
      memcpy(d, cr[i], w * w * sizeof(double));
      return;
    }
  const double *D = Dh[n];
  double cs[2 * PM + 1], sn[2 * PM + 1];
  for (int k = 0; k < w; k++) {
    cs[k] = cos((k - n) * beta);
    sn[k] = sin((k - n) * beta);
  }
  for (int i = 0; i < w; i++)
    for (int j = 0; j < w; j++) {
      double re = 0, im = 0;
      for (int k = 0; k < w; k++) {
        double v = D[k * w + i] * D[k * w + j];
        re += v * cs[k];
        im += v * sn[k];
      }
      int q = (((j - i) % 4) + 4) % 4;
      d[i * w + j] = q == 0 ? re : q == 1 ? -im : q == 2 ? -re : im;
    }
  if (ncache < (1 << 14)) {
    cn[ncache] = n;
    cb[ncache] = beta;
    cr[ncache] = (double *)malloc(w * w * sizeof(double));
    memcpy(cr[ncache], d, w * w * sizeof(double));
    ncache++;
  }
}
static void rotn(double vx, double vy, double vz, int p, double2 *G) {
  double th = atan2(hypot(vx, vy), vz), ph = atan2(vy, vx);
  double *d = (double *)malloc((2 * p + 1) * (2 * p + 1) * sizeof(double));
  for (int n = 0; n <= p; n++) {
    int w = 2 * n + 1;
    dbet(n, -th, d);
    for (int i = 0; i < w; i++)
      for (int j = 0; j < w; j++)
        G[goff(n) + i * w + j] = mk(d[j * w + i] * cos(-(j - n) * ph), d[j * w + i] * sin(-(j - n) * ph));
  }
  free(d);
}
static void irrh(double x, double y, double z, int P, double2 *I) {
  double r2 = x * x + y * y + z * z;
  double2 w = mk(x / r2, y / r2), um = mk(1 / sqrt(r2), 0);
  for (int m = 0; m <= P; m++) {
    if (m > 0)
      um = scl(mul(um, w), sqrt((2.0 * m - 1) / (2.0 * m)));
    double2 u1 = um, u2 = mk(0, 0);
    for (int n = m; n <= P; n++) {
      double2 u;
      if (n == m)
        u = um;
      else {
        u = scl(add(scl(u1, (2 * n - 1) * z / r2), scl(u2, -sqrt((double)(n - 1) * (n - 1) - (double)m * m) / r2)), 1 / sqrt((double)n * n - (double)m * m));
        u2 = u1;
        u1 = u;
      }
      I[n * n + n + m] = u;
      if (m > 0)
        I[n * n + n - m] = scl(cj(u), m % 2 ? -1 : 1);
    }
  }
}
static void grow(void **a, size_t old, size_t nw) {
  void *b;
  ck(cudaMalloc(&b, nw), "cudaMalloc");
  if (old > 0)
    ck(cudaMemcpy(b, *a, old, cudaMemcpyDeviceToDevice), "cudaMemcpy");
  if (*a != NULL)
    ck(cudaFree(*a), "cudaFree");
  *a = b;
}
static int scan(int *d, int n) {
  int last, total;
  if (n == 0)
    return 0;
  long sum = thrust::reduce(thrust::device, d, d + n, 0L);
  if (sum > 2147483647L) {
    fprintf(stderr, "fmm3: error: list size %ld overflows int\n", sum);
    exit(1);
  }
  ck(cudaMemcpy(&last, d + n - 1, sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
  thrust::exclusive_scan(thrust::device, d, d + n, d);
  ck(cudaMemcpy(&total, d + n - 1, sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
  total += last;
  ck(cudaMemcpy(d + n, &total, sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
  return total;
}
static void *pmalloc(size_t n) {
  void *p;
  return cudaMallocHost(&p, n) == cudaSuccess ? p : NULL;
}
int main(int argc, char **argv) {
  T t;
  FILE *file;
  char *ipath, *opath, *end, *ppath, *tpath;
  double A0;
  double x, *host, *dx, *out;
  float ms;
  long N, *perm, *k2, nchk, *pos;
  int lo[DM + 2], hi[DM + 2], nl, nb, cap, timing, *split, *sc, err, biot, nin;
  cudaEvent_t ev[12];
  (void)argc;
  ipath = opath = ppath = tpath = NULL;
  A0 = 1;
  t.per = 0;
  t.n1 = t.n2 = t.n3 = 1;
  t.K = 1;
  t.p = 10;
  t.s = 40;
  t.D = DM;
  t.delta = 0;
  t.ord = 2;
  char *lpath = NULL, *fpath = NULL, *gpath = NULL;
  timing = biot = 0;
  int dbg = 0;
  nchk = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    if (argv[0][1] == 'h') {
      fprintf(stderr, "Usage: fmm3 -i <in.raw> -o <out.raw> [-p <order>] [-s <max per leaf>] [-D <max depth>] [-d <gauss core delta>] [-b (Biot-Savart)] [-t (timings)] [-c <m> (direct check of the first m particles)] [-P <far table from periodic.py> -A <cube edge> (periodic)] -T <tables from tables.py> [-L <near table from neartab.py> (lattice near field) [-F <leaf table from leaftab.py> (P2M/L2P as matrices)]] [-G <table from m2ltab.py> (M2L as GEMMs on uniform levels)]\n"
                      "  Laplace: in N x (x, y, z, q), out N x (phi, gx, gy, gz); -b: in N x (x, y, z, ax, ay, az), out N x (u, v, w)\n");
      exit(1);
    }
    if (argv[0][1] == 't') {
      timing = 1;
      continue;
    }
    if (argv[0][1] == 'b') {
      biot = 1;
      continue;
    }
    if (argv[0][1] == 'g') {
      dbg = 1;
      continue;
    }
    if (argv[1] == NULL) {
      fprintf(stderr, "fmm3: error: %s needs an argument\n", argv[0]);
      exit(1);
    }
    if (argv[0][1] == 'i') {
      ipath = *++argv;
      continue;
    }
    if (argv[0][1] == 'o') {
      opath = *++argv;
      continue;
    }
    if (argv[0][1] == 'P') {
      ppath = *++argv;
      continue;
    }
    if (argv[0][1] == 'T') {
      tpath = *++argv;
      continue;
    }
    if (argv[0][1] == 'L') {
      lpath = *++argv;
      continue;
    }
    if (argv[0][1] == 'F') {
      fpath = *++argv;
      continue;
    }
    if (argv[0][1] == 'G') {
      gpath = *++argv;
      continue;
    }
    x = strtod(argv[1], &end);
    if (*end != '\0') {
      fprintf(stderr, "fmm3: error: '%s' is not a number\n", argv[1]);
      exit(1);
    }
    switch (argv[0][1]) {
    case 'p':
      t.p = x;
      break;
    case 's':
      t.s = x;
      break;
    case 'D':
      t.D = x;
      break;
    case 'd':
      t.delta = x;
      break;
    case 'c':
      nchk = x;
      break;
    case 'A':
      A0 = x;
      break;
    case 'K':
      t.K = x;
      break;
    case 'O':
      t.ord = x;
      break;
    default:
      fprintf(stderr, "fmm3: error: unknown option '%s'\n", *argv);
      exit(1);
    }
    argv++;
  }
  if (ipath == NULL || opath == NULL || t.p < 1 || t.p > PM || t.s < 1 || t.D < 2 || t.D > DM) {
    fprintf(stderr, "fmm3: error: need -i, -o, 1 <= -p <= %d, -s >= 1, 2 <= -D <= %d\n", PM, DM);
    exit(1);
  }
  double2 *Tfh = NULL;
  long ntab = 0;
  if (ppath != NULL) {
    long hd[5];
    if ((file = fopen(ppath, "r")) == NULL || fread(hd, sizeof(long), 5, file) != 5) {
      fprintf(stderr, "fmm3: error: fail to read '%s'\n", ppath);
      exit(1);
    }
    if (hd[3] != t.p || hd[4] != t.K || nchk > 0 || A0 <= 0) {
      fprintf(stderr, "fmm3: error: table '%s' is for p = %ld, K = %ld (need matching -p, -K), no -c in periodic mode, -A > 0\n", ppath, hd[3], hd[4]);
      exit(1);
    }
    t.per = 1;
    t.n1 = hd[0];
    t.n2 = hd[1];
    t.n3 = hd[2];
    if (t.D > 19)
      t.D = 19;
    ntab = (long)(2 * t.n1 - 1) * (2 * t.n2 - 1) * (2 * t.n3 - 1) * ((t.p + 1) * (t.p + 2) / 2) * (t.p + 1) * (t.p + 1);
    if ((Tfh = (double2 *)malloc(ntab * sizeof(double2))) == NULL || fread(Tfh, sizeof(double2), ntab, file) != (size_t)ntab) {
      fprintf(stderr, "fmm3: error: fail to read '%s'\n", ppath);
      exit(1);
    }
    fclose(file);
  }
  t.nc = biot ? 3 : 1;
  nin = 3 + t.nc;
  t.nt = (t.p + 1) * (t.p + 1);
  if ((file = fopen(ipath, "r")) == NULL) {
    fprintf(stderr, "fmm3: error: fail to open '%s'\n", ipath);
    exit(1);
  }
  fseek(file, 0, SEEK_END);
  N = ftell(file) / (nin * sizeof(double));
  rewind(file);
  if (N < 1) {
    fprintf(stderr, "fmm3: error: '%s' has no particles\n", ipath);
    exit(1);
  }
  if ((host = (double *)malloc(nin * N * sizeof(double))) == NULL ||
      (out = (double *)malloc(4 * N * sizeof(double))) == NULL) {
    fprintf(stderr, "fmm3: error: malloc failed\n");
    exit(1);
  }
  if (fread(host, sizeof(double), nin * N, file) != (size_t)(nin * N)) {
    fprintf(stderr, "fmm3: error: fail to read '%s'\n", ipath);
    exit(1);
  }
  fclose(file);
  if (N > 2147483647L) {
    fprintf(stderr, "fmm3: error: N = %ld overflows int indices\n", N);
    exit(1);
  }
  t.N = N;
  double mn[3], mx[3];
  for (int k = 0; k < 3; k++)
    mn[k] = mx[k] = host[k];
  for (long i = 0; i < N; i++)
    for (int k = 0; k < 3; k++) {
      mn[k] = fmin(mn[k], host[nin * i + k]);
      mx[k] = fmax(mx[k], host[nin * i + k]);
    }
  t.x0 = mn[0];
  t.y0 = mn[1];
  t.z0 = mn[2];
  t.S = fmax(mx[0] - mn[0], fmax(mx[1] - mn[1], mx[2] - mn[2])) * (1 + 1e-12) + 1e-300;
  if (t.per) {
    double Lc[3] = {t.n1 * A0, t.n2 * A0, t.n3 * A0};
    for (long i = 0; i < N; i++)
      for (int k = 0; k < 3; k++) {
        double v = fmod(host[nin * i + k], Lc[k]);
        host[nin * i + k] = v < 0 ? v + Lc[k] : v;
      }
    t.x0 = t.y0 = t.z0 = 0;
    t.S = A0;
  }
  if (tpath == NULL) {
    fprintf(stderr, "fmm3: error: need -T <tables from tables.py>\n");
    exit(1);
  }
  long tp;
  double *fe = NULL;
  {
    FILE *tf = fopen(tpath, "r");
    int P1 = t.p + 1;
    if (tf == NULL || fread(&tp, sizeof(long), 1, tf) != 1 || tp != t.p) {
      fprintf(stderr, "fmm3: error: '%s' is not a table for p = %d\n", tpath, t.p);
      exit(1);
    }
    for (int n = 0; n <= t.p; n++) {
      int w = 2 * n + 1;
      Dh[n] = (double *)malloc(w * w * sizeof(double));
      if (fread(Dh[n], sizeof(double), w * w, tf) != (size_t)(w * w)) {
        fprintf(stderr, "fmm3: error: fail to read '%s'\n", tpath);
        exit(1);
      }
    }
    double *fh = (double *)malloc(3 * P1 * P1 * P1 * sizeof(double));
    fe = (double *)malloc((size_t)t.nt * t.nt * sizeof(double));
    if (fread(fh, sizeof(double), 3 * P1 * P1 * P1, tf) != (size_t)(3 * P1 * P1 * P1) || fread(fe, sizeof(double), (size_t)t.nt * t.nt, tf) != (size_t)t.nt * t.nt) {
      fprintf(stderr, "fmm3: error: fail to read '%s'\n", tpath);
      exit(1);
    }
    fclose(tf);
    double *ft = (double *)malloc(3 * P1 * P1 * P1 * sizeof(double));
    for (int a3 = 0; a3 < 3; a3++)
      for (int n = 0; n < P1; n++)
        for (int k = 0; k < P1; k++)
          for (int m = 0; m < P1; m++)
            ft[a3 * P1 * P1 * P1 + (m * P1 + n) * P1 + k] = fh[a3 * P1 * P1 * P1 + (n * P1 + k) * P1 + m];
    ck(cudaMalloc(&t.fn, P1 * P1 * P1 * sizeof(double)), "cudaMalloc");
    ck(cudaMemcpy(t.fn, fh + 2 * P1 * P1 * P1, P1 * P1 * P1 * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    double *fsh = (double *)malloc(P1 * P1 * sizeof(double));
    for (int n = 0; n < P1; n++)
      for (int m = 0; m < P1; m++)
        fsh[n * P1 + m] = m <= n ? fh[2 * P1 * P1 * P1 + (n * P1 + m) * P1 + m] : 0;
    double *frh = (double *)malloc(P1 * P1 * sizeof(double));
    for (int k = 0; k < P1; k++)
      for (int m = 0; m < P1; m++)
        frh[k * P1 + m] = m <= k ? 1 / sqrt((double)(k - m + 1) * (k + m + 1)) : 0;
    ck(cudaMalloc(&t.fr, P1 * P1 * sizeof(double)), "cudaMalloc");
    ck(cudaMemcpy(t.fr, frh, P1 * P1 * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    free(frh);
    ck(cudaMalloc(&t.fs, P1 * P1 * sizeof(double)), "cudaMalloc");
    ck(cudaMemcpy(t.fs, fsh, P1 * P1 * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    free(fsh);
    ck(cudaMalloc(&t.fa, 3 * P1 * P1 * P1 * sizeof(double)), "cudaMalloc");
    ck(cudaMemcpy(t.fa, ft, 3 * P1 * P1 * P1 * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    free(ft);
    t.fb = t.fa + P1 * P1 * P1;
    t.fc = t.fa + 2 * P1 * P1 * P1;
  }
  for (int i = 0; i < 12; i++)
    ck(cudaEventCreate(&ev[i]), "cudaEventCreate");
  ck(cudaMalloc(&dx, nin * N * sizeof(double)), "cudaMalloc");
  ck(cudaMemcpy(dx, host, nin * N * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
  ck(cudaMalloc(&t.x, N * sizeof(double)), "cudaMalloc");
  ck(cudaMalloc(&t.y, N * sizeof(double)), "cudaMalloc");
  ck(cudaMalloc(&t.z, N * sizeof(double)), "cudaMalloc");
  ck(cudaMalloc(&t.q, t.nc * N * sizeof(double)), "cudaMalloc");
  ck(cudaMalloc(&t.key, N * sizeof(long)), "cudaMalloc");
  ck(cudaMalloc(&k2, N * sizeof(long)), "cudaMalloc");
  ck(cudaMalloc(&perm, N * sizeof(long)), "cudaMalloc");
  ck(cudaMalloc(&pos, N * sizeof(long)), "cudaMalloc");
  ck(cudaMalloc(&t.bid, N * sizeof(int)), "cudaMalloc");
  ck(cudaMalloc(&t.phi, t.nc * N * sizeof(double)), "cudaMalloc");
  ck(cudaMalloc(&t.gr, 3 * t.nc * N * sizeof(double)), "cudaMalloc");
  ck(cudaMalloc(&t.err, sizeof(int)), "cudaMalloc");
  ck(cudaMemset(t.err, 0, sizeof(int)), "cudaMemset");
  ck(cudaEventRecord(ev[0]), "cudaEventRecord");
  thrust::sequence(thrust::device, perm, perm + N);
  gather<<<grid(N), nthread>>>(N, perm, dx, t.x, nin, 0);
  gather<<<grid(N), nthread>>>(N, perm, dx, t.y, nin, 1);
  gather<<<grid(N), nthread>>>(N, perm, dx, t.z, nin, 2);
  keys<<<grid(N), nthread>>>(t);
  launched("keys");
  ck(cudaMemcpy(k2, t.key, N * sizeof(long), cudaMemcpyDeviceToDevice), "cudaMemcpy");
  thrust::stable_sort_by_key(thrust::device, k2, k2 + N, perm);
  ck(cudaMemcpy(t.key, k2, N * sizeof(long), cudaMemcpyDeviceToDevice), "cudaMemcpy");
  gather<<<grid(N), nthread>>>(N, perm, dx, t.x, nin, 0);
  gather<<<grid(N), nthread>>>(N, perm, dx, t.y, nin, 1);
  gather<<<grid(N), nthread>>>(N, perm, dx, t.z, nin, 2);
  for (int k = 0; k < t.nc; k++) {
    double *tmp;
    ck(cudaMalloc(&tmp, N * sizeof(double)), "cudaMalloc");
    gather<<<grid(N), nthread>>>(N, perm, dx, tmp, nin, 3 + k);
    ck(cudaMemcpy2D(t.q + k, t.nc * sizeof(double), tmp, sizeof(double), sizeof(double), N, cudaMemcpyDeviceToDevice), "cudaMemcpy2D");
    ck(cudaFree(tmp), "cudaFree");
  }
  inverse<<<grid(N), nthread>>>(N, perm, pos);
  launched("gather");
  ck(cudaEventRecord(ev[1]), "cudaEventRecord");
  cap = 1 << 16;
  t.lev = t.bx = t.by = t.bz = t.par = t.ch = t.start = t.count = NULL;
  t.kb = NULL;
  t.cx = t.cy = t.cz = NULL;
  grow((void **)&t.lev, 0, cap * sizeof(int));
  grow((void **)&t.bx, 0, cap * sizeof(int));
  grow((void **)&t.by, 0, cap * sizeof(int));
  grow((void **)&t.bz, 0, cap * sizeof(int));
  grow((void **)&t.par, 0, cap * sizeof(int));
  grow((void **)&t.ch, 0, cap * sizeof(int));
  grow((void **)&t.start, 0, cap * sizeof(int));
  grow((void **)&t.count, 0, cap * sizeof(int));
  grow((void **)&t.kb, 0, cap * sizeof(long));
  grow((void **)&t.cx, 0, cap * sizeof(double));
  grow((void **)&t.cy, 0, cap * sizeof(double));
  grow((void **)&t.cz, 0, cap * sizeof(double));
  int nroot = t.n1 * t.n2 * t.n3, rb[3 * 512], rp[512];
  if (nroot > 512) {
    fprintf(stderr, "fmm3: error: too many root cubes\n");
    exit(1);
  }
  for (int r = 0; r < nroot; r++) {
    rb[r] = r / (t.n2 * t.n3);
    rb[nroot + r] = (r / t.n3) % t.n2;
    rb[2 * nroot + r] = r % t.n3;
    rp[r] = -1;
  }
  ck(cudaMemcpy(t.bx, rb, nroot * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
  ck(cudaMemcpy(t.by, rb + nroot, nroot * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
  ck(cudaMemcpy(t.bz, rb + 2 * nroot, nroot * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
  ck(cudaMemcpy(t.par, rp, nroot * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
  lo[0] = 0;
  hi[0] = nroot;
  nl = 0;
  for (int l = 0;; l++) {
    int n = hi[l] - lo[l], ns;
    ck(cudaMalloc(&split, n * sizeof(int)), "cudaMalloc");
    ck(cudaMalloc(&sc, (n + 1) * sizeof(int)), "cudaMalloc");
    counts<<<grid(n), nthread>>>(t, lo[l], hi[l], l, split);
    launched("counts");
    ck(cudaMemcpy(sc, split, n * sizeof(int), cudaMemcpyDeviceToDevice), "cudaMemcpy");
    ns = scan(sc, n);
    nl = l;
    if (ns == 0) {
      ck(cudaFree(split), "cudaFree");
      ck(cudaFree(sc), "cudaFree");
      break;
    }
    if ((long)hi[l] + 8L * ns > 2000000000L) {
      fprintf(stderr, "fmm3: error: too many boxes\n");
      exit(1);
    }
    if (hi[l] + 8 * ns > cap) {
      int nc = cap;
      while (nc < hi[l] + 8 * ns)
        nc *= 2;
      grow((void **)&t.lev, cap * sizeof(int), nc * sizeof(int));
      grow((void **)&t.bx, cap * sizeof(int), nc * sizeof(int));
      grow((void **)&t.by, cap * sizeof(int), nc * sizeof(int));
      grow((void **)&t.bz, cap * sizeof(int), nc * sizeof(int));
      grow((void **)&t.par, cap * sizeof(int), nc * sizeof(int));
      grow((void **)&t.ch, cap * sizeof(int), nc * sizeof(int));
      grow((void **)&t.start, cap * sizeof(int), nc * sizeof(int));
      grow((void **)&t.count, cap * sizeof(int), nc * sizeof(int));
      grow((void **)&t.kb, cap * sizeof(long), nc * sizeof(long));
      grow((void **)&t.cx, cap * sizeof(double), nc * sizeof(double));
      grow((void **)&t.cy, cap * sizeof(double), nc * sizeof(double));
      grow((void **)&t.cz, cap * sizeof(double), nc * sizeof(double));
      cap = nc;
    }
    children<<<grid(n), nthread>>>(t, lo[l], hi[l], split, sc);
    launched("children");
    lo[l + 1] = hi[l];
    hi[l + 1] = hi[l] + 8 * ns;
    ck(cudaFree(split), "cudaFree");
    ck(cudaFree(sc), "cudaFree");
  }
  nb = hi[nl];
  if ((long)nb * 125 > 2147483647L) {
    fprintf(stderr, "fmm3: error: %d boxes overflow the list encoding\n", nb);
    exit(1);
  }
  if ((long)t.n1 * t.n2 * t.n3 >= (1L << (63 - 3 * t.D))) {
    fprintf(stderr, "fmm3: error: Morton key overflows: %d roots at depth %d\n", t.n1 * t.n2 * t.n3, t.D);
    exit(1);
  }
  ck(cudaMemcpyToSymbol(clo, lo, sizeof(lo)), "cudaMemcpyToSymbol");
  ck(cudaMemcpyToSymbol(chi, hi, sizeof(hi)), "cudaMemcpyToSymbol");
  ck(cudaMalloc(&t.uo, (nb + 1) * sizeof(int)), "cudaMalloc");
  ck(cudaMalloc(&t.vo, (nb + 1) * sizeof(int)), "cudaMalloc");
  ck(cudaMalloc(&t.wo, (nb + 1) * sizeof(int)), "cudaMalloc");
  ck(cudaMalloc(&t.xo, (nb + 1) * sizeof(int)), "cudaMalloc");
  t.ul = t.vl = t.wl = t.xl = NULL;
  lists<<<grid(nb), nthread>>>(t, nb, 0);
  launched("lists count");
  int nU = scan(t.uo, nb), nV = scan(t.vo, nb), nW = scan(t.wo, nb), nX = scan(t.xo, nb);
  ck(cudaMalloc(&t.ul, (nU + 1) * sizeof(int)), "cudaMalloc");
  ck(cudaMalloc(&t.vl, (nV + 1) * sizeof(int)), "cudaMalloc");
  ck(cudaMalloc(&t.wl, (nW + 1) * sizeof(int)), "cudaMalloc");
  ck(cudaMalloc(&t.xl, (nX + 1) * sizeof(int)), "cudaMalloc");
  lists<<<grid(nb), nthread>>>(t, nb, 1);
  launched("lists fill");
  setbid<<<grid(nb), nthread>>>(t, nb);
  launched("setbid");
  ck(cudaMemcpy(&err, t.err, sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
  if (err) {
    fprintf(stderr, "fmm3: error: list stack overflow (raise STK)\n");
    exit(1);
  }
  ck(cudaEventRecord(ev[2]), "cudaEventRecord");
  cublasHandle_t hb = NULL;
  if ((lpath != NULL || gpath != NULL) && cublasCreate(&hb) != CUBLAS_STATUS_SUCCESS) {
    fprintf(stderr, "fmm3: error: cublasCreate failed\n");
    exit(1);
  }
  int lB = 0, lB3 = 0, lKc = 0, nlv = 0, ll = -1, nof = 0, ffm = 0, nbt = 0, *lst = NULL, *dst = NULL, *dsrca = NULL, *dlid = NULL;
  long lM3 = 0, lMc3 = 0;
  int3 lax = make_int3(0, 1, 2), lD = make_int3(0, 0, 0);
  double *ul = NULL, hl = 0, *dG = NULL, *dX = NULL, *dY = NULL, *dQ = NULL, *dP2 = NULL, *dV2 = NULL;
  double *fP = NULL, *fR = NULL;
  double2 *fA = NULL, *fU = NULL, *fG = NULL;
  cufftHandle pf, pb, pg;
  if (lpath != NULL) {
    long hd[3];
    double eh;
    FILE *lf = fopen(lpath, "r");
    if (lf == NULL || fread(hd, sizeof(long), 3, lf) != 3 || fread(&eh, sizeof(double), 1, lf) != 1) {
      fprintf(stderr, "fmm3: error: fail to read '%s'\n", lpath);
      exit(1);
    }
    lB = hd[0];
    lKc = hd[1];
    lB3 = lB * lB * lB;
    int R = hd[2], W = 2 * R + 1;
    long nw = (long)W * W * W * 3;
    double *tf = (double *)malloc(2 * nw * sizeof(double));
    if (tf == NULL || fread(tf, sizeof(double), 2 * nw, lf) != (size_t)(2 * nw) || !biot) {
      fprintf(stderr, "fmm3: error: fail to read '%s' (or no -b)\n", lpath);
      exit(1);
    }
    fclose(lf);
    int *chh = (int *)malloc(nb * sizeof(int)), *cnt2 = (int *)malloc(nb * sizeof(int)), *lvh = (int *)malloc(nb * sizeof(int)), *st2 = (int *)malloc(nb * sizeof(int)), *bh = (int *)malloc(3 * nb * sizeof(int));
    double *xs = (double *)pmalloc(3 * N * sizeof(double)), *cxh = (double *)malloc(3 * nb * sizeof(double));
    ck(cudaMemcpy(chh, t.ch, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(cnt2, t.count, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(lvh, t.lev, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(st2, t.start, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(bh, t.bx, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(bh + nb, t.by, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(bh + 2 * nb, t.bz, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(xs, t.x, N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(xs + N, t.y, N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(xs + 2 * N, t.z, N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(cxh, t.cx, nb * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(cxh + nb, t.cy, nb * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(cxh + 2 * nb, t.cz, nb * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    for (int b = 0; b < nb; b++)
      if (chh[b] < 0) {
        if (ll < 0)
          ll = lvh[b];
        if (lvh[b] != ll || cnt2[b] != lB3) {
          fprintf(stderr, "fmm3: error: -L needs uniform leaves of %d particles (leaf %d level %d count %d)\n", lB3, b, lvh[b], cnt2[b]);
          exit(1);
        }
        nlv++;
      }
    hl = t.S / (double)(1L << ll);
    double hh = hl / lB;
    int *lc = (int *)malloc(3 * lB3 * sizeof(int));
    for (int i = 0; i < lB3; i++) {
      int a3[3] = {0, 0, 0};
      for (int bb = 0; (1 << bb) < lB; bb++) {
        a3[0] |= ((i >> (3 * bb + 2)) & 1) << bb;
        a3[1] |= ((i >> (3 * bb + 1)) & 1) << bb;
        a3[2] |= ((i >> (3 * bb)) & 1) << bb;
      }
      for (int d = 0; d < 3; d++)
        lc[3 * i + d] = a3[d];
    }
    int *lid = (int *)malloc(nlv * sizeof(int));
    long NX = (long)t.n1 << ll, NY = (long)t.n2 << ll, NZ = (long)t.n3 << ll;
    int *map = (int *)malloc(NX * NY * NZ * sizeof(int));
    for (long i = 0; i < NX * NY * NZ; i++)
      map[i] = -1;
    for (int b = 0, k = 0; b < nb; b++)
      if (chh[b] < 0) {
        lid[k++] = b;
        map[((long)bh[b] * NY + bh[nb + b]) * NZ + bh[2 * nb + b]] = b;
        for (int i = 0; i < lB3; i++)
          for (int d = 0; d < 3; d++) {
            double v = (xs[d * N + st2[b] + i] - (cxh[d * nb + b] - hl / 2)) / hh - 0.5 - lc[3 * i + d];
            if (fabs(v) > 1e-6) {
              fprintf(stderr, "fmm3: error: -L needs cell-centred lattice particles (leaf %d node %d off %g)\n", b, i, v);
              exit(1);
            }
          }
      }
    int no2 = 2 * lKc + 1;
    nof = no2 * no2 * no2;
    double *Gh = (double *)malloc((size_t)3 * lB3 * lB3 * sizeof(double));
    ck(cudaMalloc(&dG, (size_t)nof * 3 * lB3 * lB3 * sizeof(double)), "cudaMalloc");
    for (int o = 0; o < nof; o++) {
      int ox = o / (no2 * no2) - lKc, oy = (o / no2) % no2 - lKc, oz = o % no2 - lKc;
      int dfl = abs(ox) > 1 || abs(oy) > 1 || abs(oz) > 1;
      for (int c = 0; c < 3; c++)
        for (int j = 0; j < lB3; j++)
          for (int i = 0; i < lB3; i++) {
            int dx = lc[3 * i] - ox * lB - lc[3 * j], dy = lc[3 * i + 1] - oy * lB - lc[3 * j + 1], dz = lc[3 * i + 2] - oz * lB - lc[3 * j + 2];
            Gh[(size_t)c * lB3 * lB3 + (size_t)j * lB3 + i] = tf[(size_t)dfl * nw + ((((size_t)(dx + R) * W + dy + R) * W + dz + R) * 3 + c)] / (hh * hh);
          }
      ck(cudaMemcpy(dG + (size_t)o * 3 * lB3 * lB3, Gh, (size_t)3 * lB3 * lB3 * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    }
    free(Gh);
    ck(cudaMalloc(&dst, nb * sizeof(int)), "cudaMalloc");
    ck(cudaMemcpy(dst, st2, nb * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
    ck(cudaMalloc(&ul, 3 * N * sizeof(double)), "cudaMalloc");
    ck(cudaMemset(ul, 0, 3 * N * sizeof(double)), "cudaMemset");
    ck(cudaMalloc(&lst, nlv * sizeof(int)), "cudaMalloc");
    int *lsth = (int *)malloc(nlv * sizeof(int));
    for (int k = 0; k < nlv; k++)
      lsth[k] = st2[lid[k]];
    ck(cudaMemcpy(lst, lsth, nlv * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
    ck(cudaMalloc(&dlid, nlv * sizeof(int)), "cudaMalloc");
    ck(cudaMemcpy(dlid, lid, nlv * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
    int *srca = (int *)malloc((size_t)nof * nlv * sizeof(int));
    for (int o = 0; o < nof; o++) {
      int ox = o / (no2 * no2) - lKc, oy = (o / no2) % no2 - lKc, oz = o % no2 - lKc;
      for (int k = 0; k < nlv; k++) {
        int b = lid[k];
        long X = bh[b] + ox, Y = bh[nb + b] + oy, Z = bh[2 * nb + b] + oz;
        if (t.per) {
          X = ((X % NX) + NX) % NX;
          Y = ((Y % NY) + NY) % NY;
          Z = ((Z % NZ) + NZ) % NZ;
          srca[(size_t)o * nlv + k] = map[(X * NY + Y) * NZ + Z];
        } else
          srca[(size_t)o * nlv + k] = (X < 0 || Y < 0 || Z < 0 || X >= NX || Y >= NY || Z >= NZ) ? -1 : map[(X * NY + Y) * NZ + Z];
      }
    }
    ck(cudaMalloc(&dsrca, (size_t)nof * nlv * sizeof(int)), "cudaMalloc");
    ck(cudaMemcpy(dsrca, srca, (size_t)nof * nlv * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
    ffm = getenv("FMM_NEAR") != NULL && !strcmp(getenv("FMM_NEAR"), "fft");
    if (ffm && lKc != 1) {
      fprintf(stderr, "fmm3: error: FMM_NEAR=fft needs a Kc = 1 table\n");
      exit(1);
    }
    int M = 4 * lB;
    lM3 = (long)M * M * M;
    lMc3 = (long)M * M * (M / 2 + 1);
    if (ffm) {
      nbt = nlv < 512 ? nlv : 512;
      double *kr = (double *)calloc(3 * lM3, sizeof(double)), *dkr;
      for (int dx = -(2 * lB - 1); dx <= 2 * lB - 1; dx++)
        for (int dy = -(2 * lB - 1); dy <= 2 * lB - 1; dy++)
          for (int dz = -(2 * lB - 1); dz <= 2 * lB - 1; dz++)
            for (int c = 0; c < 3; c++)
              kr[c * lM3 + ((long)((dx + M) % M) * M + (dy + M) % M) * M + (dz + M) % M] = tf[((((size_t)(dx + R) * W + dy + R) * W + dz + R) * 3 + c)] / (hh * hh);
      ck(cudaMalloc(&dkr, 3 * lM3 * sizeof(double)), "cudaMalloc");
      ck(cudaMemcpy(dkr, kr, 3 * lM3 * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
      ck(cudaMalloc(&fG, 3 * lMc3 * sizeof(double2)), "cudaMalloc");
      int nn[3] = {M, M, M};
      if (cufftPlanMany(&pg, 3, nn, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z, 3) != CUFFT_SUCCESS || cufftExecD2Z(pg, dkr, (cufftDoubleComplex *)fG) != CUFFT_SUCCESS ||
          cufftPlanMany(&pf, 3, nn, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z, 3 * nbt) != CUFFT_SUCCESS || cufftPlanMany(&pb, 3, nn, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D, 3 * nbt) != CUFFT_SUCCESS) {
        fprintf(stderr, "fmm3: error: cufft plan failed\n");
        exit(1);
      }
      ck(cudaMalloc(&fP, (size_t)nbt * 3 * lM3 * sizeof(double)), "cudaMalloc");
      ck(cudaMalloc(&fR, (size_t)nbt * 3 * lM3 * sizeof(double)), "cudaMalloc");
      ck(cudaMalloc(&fA, (size_t)nbt * 3 * lMc3 * sizeof(double2)), "cudaMalloc");
      ck(cudaMalloc(&fU, (size_t)nbt * 3 * lMc3 * sizeof(double2)), "cudaMalloc");
    } else {
      int P[3] = {(int)NX, (int)NY, (int)NZ}, a[3] = {0, 1, 2};
      for (int i = 0; i < 3; i++)
        for (int k = 0; k < 2; k++)
          if (P[a[k]] > P[a[k + 1]]) {
            int w = a[k];
            a[k] = a[k + 1];
            a[k + 1] = w;
          }
      lax = make_int3(a[0], a[1], a[2]);
      lD = make_int3(P[a[0]], P[a[1]], P[a[2]]);
      long ncp = (long)(lD.x + 2 * lKc) * (lD.y + 2 * lKc) * (lD.z + 2 * lKc);
      ck(cudaMalloc(&dX, (size_t)lB3 * 3 * ncp * sizeof(double)), "cudaMalloc");
      ck(cudaMalloc(&dY, (size_t)lB3 * 6 * ncp * sizeof(double)), "cudaMalloc");
    }
    if (fpath != NULL) {
      long fh[2];
      FILE *ff = fopen(fpath, "r");
      size_t np2 = (size_t)2 * t.nt * lB3, nv2 = (size_t)3 * lB3 * 6 * t.nt;
      double *P2 = (double *)malloc(np2 * sizeof(double)), *V2 = (double *)malloc(nv2 * sizeof(double));
      if (ff == NULL || fread(fh, sizeof(long), 2, ff) != 2 || fh[0] != lB || fh[1] != t.p || fread(P2, sizeof(double), np2, ff) != np2 || fread(V2, sizeof(double), nv2, ff) != nv2) {
        fprintf(stderr, "fmm3: error: '%s' is not a leaf table for B = %d, p = %d\n", fpath, lB, t.p);
        exit(1);
      }
      fclose(ff);
      if (ll != nl || nlv != hi[ll] - lo[ll]) {
        fprintf(stderr, "fmm3: error: -F needs every box of the last level to be a leaf\n");
        exit(1);
      }
      for (int k = 0; k < nlv; k++)
        if (lid[k] != lo[ll] + k || st2[lid[k]] != k * lB3) {
          fprintf(stderr, "fmm3: error: -F needs leaves in particle order (leaf %d)\n", k);
          exit(1);
        }
      ck(cudaMalloc(&dP2, np2 * sizeof(double)), "cudaMalloc");
      ck(cudaMemcpy(dP2, P2, np2 * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
      ck(cudaMalloc(&dV2, nv2 * sizeof(double)), "cudaMalloc");
      ck(cudaMemcpy(dV2, V2, nv2 * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
      ck(cudaMalloc(&dQ, (size_t)lB3 * 3 * nlv * sizeof(double)), "cudaMalloc");
      free(P2);
      free(V2);
    }
    free(tf);
    free(chh);
    free(cnt2);
    free(lvh);
    free(st2);
    free(bh);
    cudaFreeHost(xs);
    free(cxh);
    free(lc);
    free(lid);
    free(map);
    free(lsth);
    free(srca);
  }
  int m2g = 0, H = 16 * ((t.p + 1) * (t.p + 2) / 2), lg0 = t.per ? 1 : 2;
  double *dA = NULL, *dS = NULL, *dLs = NULL;
  int mode = getenv("FMM_M2L") != NULL ? atoi(getenv("FMM_M2L")) : gpath != NULL ? 3 : 1;
  if (mode == 3 && gpath == NULL) {
    fprintf(stderr, "fmm3: error: FMM_M2L=3 needs -G\n");
    exit(1);
  }
  if (mode == 3) {
    int full = 1;
    for (int l = 0; l <= nl; l++)
      full = full && (long)(hi[l] - lo[l]) == ((long)nroot << (3 * l));
    if (!full || nl < lg0) {
      fprintf(stderr, "fmm3: error: -G needs a uniform tree (every level full) of depth >= %d\n", lg0);
      exit(1);
    }
    long hp, nj = (long)(2 * t.p + 1) * (2 * t.p + 2) / 2;
    double2 *Ih = (double2 *)malloc(343 * nj * sizeof(double2)), *dI;
    FILE *gf = fopen(gpath, "r");
    if (gf == NULL || fread(&hp, sizeof(long), 1, gf) != 1 || hp != t.p || fread(Ih, sizeof(double2), 343 * nj, gf) != (size_t)(343 * nj)) {
      fprintf(stderr, "fmm3: error: '%s' is not an M2L table for p = %d\n", gpath, t.p);
      exit(1);
    }
    fclose(gf);
    double *dfe;
    ck(cudaMalloc(&dI, 343 * nj * sizeof(double2)), "cudaMalloc");
    ck(cudaMemcpy(dI, Ih, 343 * nj * sizeof(double2), cudaMemcpyHostToDevice), "cudaMemcpy");
    ck(cudaMalloc(&dfe, (size_t)t.nt * t.nt * sizeof(double)), "cudaMalloc");
    ck(cudaMemcpy(dfe, fe, (size_t)t.nt * t.nt * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    ck(cudaMalloc(&dA, 26 * (size_t)H * H * sizeof(double)), "cudaMalloc");
    mkA<<<grid(26L * H * H), nthread>>>(t.p, dI, dfe, dA);
    launched("mkA");
    long nc = (long)((t.n1 << nl) / 2 + 2) * ((t.n2 << nl) / 2 + 2) * ((t.n3 << nl) / 2 + 2) * t.nc * H;
    ck(cudaMalloc(&dS, nc * sizeof(double)), "cudaMalloc");
    ck(cudaMalloc(&dLs, nc * sizeof(double)), "cudaMalloc");
    ck(cudaFree(dI), "cudaFree");
    ck(cudaFree(dfe), "cudaFree");
    free(Ih);
    m2g = 1;
  }
  long gs = goff(t.p + 1);
  double2 *A = (double2 *)malloc(343 * gs * sizeof(double2)), *B = (double2 *)malloc(8 * gs * sizeof(double2));
  if (A == NULL || B == NULL) {
    fprintf(stderr, "fmm3: error: malloc failed\n");
    exit(1);
  }
  memset(A, 0, 343 * gs * sizeof(double2));
  for (int o = 0; o < 343; o++) {
    int ox = o / 49 - 3, oy = (o / 7) % 7 - 3, oz = o % 7 - 3;
    if (abs(ox) < 2 && abs(oy) < 2 && abs(oz) < 2)
      continue;
    rotn(-ox, -oy, -oz, t.p, A + o * gs);
  }
  ck(cudaMalloc(&t.Gv, 343 * gs * sizeof(double2)), "cudaMalloc");
  ck(cudaMemcpy(t.Gv, A, 343 * gs * sizeof(double2), cudaMemcpyHostToDevice), "cudaMemcpy");
  ck(cudaMalloc(&t.Gm, 8 * gs * sizeof(double2)), "cudaMalloc");
  ck(cudaMalloc(&t.Gl, 8 * gs * sizeof(double2)), "cudaMalloc");
  for (int q = 0; q < 8; q++)
    rotn(((q >> 2) & 1) - 0.5, ((q >> 1) & 1) - 0.5, (q & 1) - 0.5, t.p, B + q * gs);
  ck(cudaMemcpy(t.Gm, B, 8 * gs * sizeof(double2), cudaMemcpyHostToDevice), "cudaMemcpy");
  for (int q = 0; q < 8; q++)
    rotn(0.5 - ((q >> 2) & 1), 0.5 - ((q >> 1) & 1), 0.5 - (q & 1), t.p, B + q * gs);
  ck(cudaMemcpy(t.Gl, B, 8 * gs * sizeof(double2), cudaMemcpyHostToDevice), "cudaMemcpy");
  size_t shm = 3 * t.nt * sizeof(double2);
  ck(cudaFuncSetAttribute(m2m, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shm), "cudaFuncSetAttribute");
  ck(cudaFuncSetAttribute(m2lh, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shm), "cudaFuncSetAttribute");
  ck(cudaFuncSetAttribute(l2l, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shm), "cudaFuncSetAttribute");
  ck(cudaMalloc(&t.M, (size_t)nb * t.nc * t.nt * sizeof(double2)), "cudaMalloc");
  ck(cudaMalloc(&t.L, (size_t)nb * t.nc * t.nt * sizeof(double2)), "cudaMalloc");
  ck(cudaMemset(t.M, 0, (size_t)nb * t.nc * t.nt * sizeof(double2)), "cudaMemset");
  ck(cudaMemset(t.L, 0, (size_t)nb * t.nc * t.nt * sizeof(double2)), "cudaMemset");
  ck(cudaEventRecord(ev[3]), "cudaEventRecord");
  if (dP2 != NULL) {
    double one = 1, zero = 0;
    lgath<<<grid((long)nlv * lB3 * 3), nthread>>>(nlv, lB3, dlid, dst, t.q, dQ);
    if (cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_N, 2 * t.nt, 3 * nlv, lB3, &one, dP2, 2 * t.nt, dQ, lB3, &zero, (double *)(t.M + (size_t)lo[ll] * 3 * t.nt), 2 * t.nt) != CUBLAS_STATUS_SUCCESS) {
      fprintf(stderr, "fmm3: error: cublasDgemm failed\n");
      exit(1);
    }
  } else
    p2m<<<grid((long)nb * t.nc * (t.p + 1)), nthread>>>(t, nb);
  launched("p2m");
  for (int l = nl - 1; l >= 0; l--) {
    m2m<<<hi[l] - lo[l], TB, shm>>>(t, lo[l], hi[l], l);
    launched("m2m");
  }
  if (t.per) {
    double2 *Tf;
    ck(cudaMalloc(&Tf, ntab * sizeof(double2)), "cudaMalloc");
    ck(cudaMemcpy(Tf, Tfh, ntab * sizeof(double2), cudaMemcpyHostToDevice), "cudaMemcpy");
    farop<<<grid((long)nroot * t.nc * ((t.p + 1) * (t.p + 2) / 2)), nthread>>>(t, nroot, Tf);
    launched("farop");
    ck(cudaFree(Tf), "cudaFree");
  }
  ck(cudaEventRecord(ev[4]), "cudaEventRecord");
  if (m2g) {
    if (t.per && t.K > 1)
      m2lh<<<hi[0], M2LT, shm>>>(t, hi[0]);
    for (int l = lg0; l <= nl; l++) {
      int P[3] = {(t.n1 << l) / 2, (t.n2 << l) / 2, (t.n3 << l) / 2}, a[3] = {0, 1, 2};
      for (int i = 0; i < 3; i++)
        for (int k = 0; k < 2; k++)
          if (P[a[k]] > P[a[k + 1]]) {
            int w = a[k];
            a[k] = a[k + 1];
            a[k + 1] = w;
          }
      int3 ax = make_int3(a[0], a[1], a[2]), D = make_int3(P[a[0]], P[a[1]], P[a[2]]);
      long c0 = (((long)D.y + 2) * (D.z + 2) + D.z + 2 + 1) * t.nc, nh = H / 16;
      long ncol = (((long)D.y - 1) * (D.z + 2) + D.z) * t.nc, sb = ((long)D.y + 2) * (D.z + 2) * t.nc * H;
      mpack<<<grid((long)(D.x + 2) * (D.y + 2) * (D.z + 2) * t.nc * 8 * nh), nthread>>>(t, l, lo[l], ax, D, dS);
      for (int d = 0, k = 0; d < 27; d++) {
        if (d == 13)
          continue;
        int dd[3] = {d / 9 - 1, (d / 3) % 3 - 1, d % 3 - 1};
        long sh = ((long)dd[a[0]] * (D.y + 2) + dd[a[1]]) * (D.z + 2) + dd[a[2]];
        double one = 1, beta = k > 0;
        if (cublasDgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N, H, ncol, H, &one, dA + (size_t)k * H * H, H, 0, dS + (c0 + sh * t.nc) * H, H, sb, &beta, dLs + c0 * H, H, sb, D.x) != CUBLAS_STATUS_SUCCESS) {
          fprintf(stderr, "fmm3: error: cublasDgemm failed\n");
          exit(1);
        }
        k++;
      }
      munpack<<<grid((long)D.x * D.y * D.z * t.nc * 8 * nh), nthread>>>(t, l, lo[l], ax, D, dLs);
    }
  } else if (mode != 2)
    m2lh<<<nb, M2LT, shm>>>(t, nb);
  else {
    int nh = (t.p + 1) * (t.p + 2) / 2, P2 = 2 * t.p + 1;
    int *vo = (int *)malloc((nb + 1) * sizeof(int)), *vl = (int *)malloc((nV + 1) * sizeof(int)),
        *bh = (int *)malloc(3 * nb * sizeof(int)), *lvh = (int *)malloc(nb * sizeof(int));
    ck(cudaMemcpy(vo, t.vo, (nb + 1) * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(vl, t.vl, nV * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(bh, t.bx, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(bh + nb, t.by, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(bh + 2 * nb, t.bz, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(lvh, t.lev, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    int cnt[343], *pos = (int *)malloc(344 * sizeof(int)), *pt = (int *)malloc((nV + 1) * 3 * sizeof(int)), maxk = 0;
    memset(cnt, 0, sizeof(cnt));
    for (int b = 0; b < nb; b++)
      for (int e = vo[b]; e < vo[b + 1]; e++) {
        int sidx = vl[e] / 125, sc = vl[e] % 125;
        int ox = bh[sidx] + (sc / 25 - 2) * (t.n1 << lvh[b]) - bh[b], oy = bh[nb + sidx] + ((sc / 5) % 5 - 2) * (t.n2 << lvh[b]) - bh[nb + b], oz = bh[2 * nb + sidx] + (sc % 5 - 2) * (t.n3 << lvh[b]) - bh[2 * nb + b];
        cnt[(ox + 3) * 49 + (oy + 3) * 7 + (oz + 3)]++;
      }
    pos[0] = 0;
    for (int o = 0; o < 343; o++) {
      pos[o + 1] = pos[o] + cnt[o];
      maxk = cnt[o] > maxk ? cnt[o] : maxk;
      cnt[o] = 0;
    }
    for (int b = 0; b < nb; b++)
      for (int e = vo[b]; e < vo[b + 1]; e++) {
        int sidx = vl[e] / 125, sc = vl[e] % 125;
        int o = (bh[sidx] + (sc / 25 - 2) * (t.n1 << lvh[b]) - bh[b] + 3) * 49 + (bh[nb + sidx] + ((sc / 5) % 5 - 2) * (t.n2 << lvh[b]) - bh[nb + b] + 3) * 7 + (bh[2 * nb + sidx] + (sc % 5 - 2) * (t.n3 << lvh[b]) - bh[2 * nb + b] + 3);
        int q = pos[o] + cnt[o]++;
        pt[q] = b;
        pt[nV + q] = sidx;
        pt[2 * nV + q] = lvh[b];
      }
    int *dpt;
    double2 *Th = (double2 *)malloc((size_t)P2 * P2 * sizeof(double2)), *dTh, *Tm, *X, *Y;
    ck(cudaMalloc(&dpt, (size_t)3 * (nV + 1) * sizeof(int)), "cudaMalloc");
    ck(cudaMemcpy(dpt, pt, (size_t)3 * (nV + 1) * sizeof(int), cudaMemcpyHostToDevice), "cudaMemcpy");
    ck(cudaMalloc(&dTh, (size_t)P2 * P2 * sizeof(double2)), "cudaMalloc");
    double *dfe, *feh = (double *)malloc((size_t)t.nt * nh * sizeof(double));
    for (int ci = 0; ci < t.nt; ci++)
      for (int n = 0, r = 0; n <= t.p; n++)
        for (int m = 0; m <= n; m++, r++)
          feh[(long)ci * nh + r] = fe[(long)(n * n + n + m) * t.nt + ci];
    ck(cudaMalloc(&dfe, (size_t)t.nt * nh * sizeof(double)), "cudaMalloc");
    ck(cudaMemcpy(dfe, feh, (size_t)t.nt * nh * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    free(feh);
    ck(cudaMalloc(&Tm, (size_t)nh * t.nt * sizeof(double2)), "cudaMalloc");
    ck(cudaMalloc(&X, (size_t)t.nt * (maxk + 1) * sizeof(double2)), "cudaMalloc");
    ck(cudaMalloc(&Y, (size_t)nh * (maxk + 1) * sizeof(double2)), "cudaMalloc");
    cublasHandle_t hb;
    if (cublasCreate(&hb) != CUBLAS_STATUS_SUCCESS) {
      fprintf(stderr, "fmm3: error: cublasCreate failed\n");
      exit(1);
    }
    cuDoubleComplex one = {1, 0}, zero = {0, 0};
    for (int o = 0; o < 343; o++) {
      int K = pos[o + 1] - pos[o];
      if (K == 0)
        continue;
      int ox = o / 49 - 3, oy = (o / 7) % 7 - 3, oz = o % 7 - 3;
      irrh(-ox, -oy, -oz, 2 * t.p, Th);
      ck(cudaMemcpy(dTh, Th, (size_t)P2 * P2 * sizeof(double2), cudaMemcpyHostToDevice), "cudaMemcpy");
      dense<<<grid((long)nh * t.nt), nthread>>>(t.p, dTh, dfe, Tm);
      for (int c = 0; c < t.nc; c++) {
        gath<<<grid((long)K * t.nt), nthread>>>(t, dpt + nV + pos[o], K, c, X);
        if (cublasZgemm(hb, CUBLAS_OP_N, CUBLAS_OP_N, nh, K, t.nt, &one, (cuDoubleComplex *)Tm, nh, (cuDoubleComplex *)X, t.nt, &zero, (cuDoubleComplex *)Y, nh) != CUBLAS_STATUS_SUCCESS) {
          fprintf(stderr, "fmm3: error: cublasZgemm failed\n");
          exit(1);
        }
        scat<<<grid((long)K * nh), nthread>>>(t, dpt + pos[o], K, c, Y);
      }
    }
    launched("m2l dense");
    cublasDestroy(hb);
  }
  launched("m2l");
  ck(cudaEventRecord(ev[5]), "cudaEventRecord");
  p2l<<<grid((long)nb * t.nc * (t.p + 1)), nthread>>>(t, nb);
  launched("p2l");
  for (int l = 1; l <= nl; l++) {
    l2l<<<hi[l] - lo[l], TB, shm>>>(t, lo[l], hi[l], l);
    launched("l2l");
  }
  ck(cudaEventRecord(ev[6]), "cudaEventRecord");
  if (dbg) {
    size_t ne = (size_t)nb * t.nc * t.nt;
    double2 *hm = (double2 *)malloc(ne * sizeof(double2)), *hl = (double2 *)malloc(ne * sizeof(double2));
    int *lv = (int *)malloc(nb * sizeof(int));
    ck(cudaMemcpy(hm, t.M, ne * sizeof(double2), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(hl, t.L, ne * sizeof(double2), cudaMemcpyDeviceToHost), "cudaMemcpy");
    ck(cudaMemcpy(lv, t.lev, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
    long bad = 0, sub = 0;
    double gmx[2] = {0, 0}, gmn[2] = {1e308, 1e308};
    printf("DBG kind level n max|c| min|c|(nonzero)\n");
    for (int k2 = 0; k2 < 2; k2++) {
      double2 *a2 = k2 ? hl : hm;
      for (int l = 0; l <= nl; l++)
        for (int n = 0; n <= t.p; n++) {
          double mxv = 0, mnv = 1e308;
          for (int b = lo[l]; b < hi[l]; b++)
            for (int c = 0; c < t.nc; c++)
              for (int m = -n; m <= n; m++) {
                double2 v = a2[((size_t)b * t.nc + c) * t.nt + n * n + n + m];
                double av = hypot(v.x, v.y);
                if (!isfinite(v.x) || !isfinite(v.y))
                  bad++;
                if ((v.x != 0 && fabs(v.x) < 2.2250738585072014e-308) || (v.y != 0 && fabs(v.y) < 2.2250738585072014e-308))
                  sub++;
                if (av > mxv)
                  mxv = av;
                if (av > 0 && av < mnv)
                  mnv = av;
              }
          if (mxv > gmx[k2])
            gmx[k2] = mxv;
          if (mnv < gmn[k2])
            gmn[k2] = mnv;
          printf("DBG %c %d %d %.6e %.6e\n", k2 ? 'L' : 'M', l, n, mxv, mnv < 1e308 ? mnv : 0.0);
        }
    }
    printf("DBG summary max|M| %.6e min|M| %.6e max|L| %.6e min|L| %.6e nonfinite %ld subnormal %ld\n", gmx[0], gmn[0], gmx[1], gmn[1], bad, sub);
    free(hm);
    free(hl);
    free(lv);
  }
  if (dV2 != NULL) {
    double one = 1, a2 = 1 / (hl * hl);
    ck(cudaMemset(t.phi, 0, t.nc * N * sizeof(double)), "cudaMemset");
    ck(cudaMemset(t.gr, 0, 3 * t.nc * N * sizeof(double)), "cudaMemset");
    if (cublasDgemm(hb, CUBLAS_OP_N, CUBLAS_OP_N, 3 * lB3, nlv, 6 * t.nt, &a2, dV2, 3 * lB3, (double *)(t.L + (size_t)lo[ll] * 3 * t.nt), 6 * t.nt, &one, ul, 3 * lB3) != CUBLAS_STATUS_SUCCESS) {
      fprintf(stderr, "fmm3: error: cublasDgemm failed\n");
      exit(1);
    }
  } else
    l2p<<<grid(N), nthread>>>(t);
  launched("l2p");
  m2p<<<grid(N), nthread>>>(t);
  launched("m2p");
  ck(cudaEventRecord(ev[7]), "cudaEventRecord");
  if (lpath == NULL) {
    p2p<<<nb, nwarp>>>(t, nb);
    launched("p2p");
  } else {
    cudaEvent_t e0, e1;
    ck(cudaEventCreate(&e0), "cudaEventCreate");
    ck(cudaEventCreate(&e1), "cudaEventCreate");
    ck(cudaEventRecord(e0), "cudaEventRecord");
    if (ffm)
      for (int k0 = 0; k0 < nlv; k0 += nbt) {
        int nc2 = nlv - k0 < nbt ? nlv - k0 : nbt;
        fpatch<<<grid((long)nc2 * 3 * lM3), nthread>>>(nc2, k0, nlv, lB, dsrca, dst, t.q, fP);
        if (cufftExecD2Z(pf, fP, (cufftDoubleComplex *)fA) != CUFFT_SUCCESS) {
          fprintf(stderr, "fmm3: error: cufft forward failed\n");
          exit(1);
        }
        fmulk<<<grid((long)nc2 * lMc3), nthread>>>(nc2, lMc3, fG, fA, fU);
        if (cufftExecZ2D(pb, (cufftDoubleComplex *)fU, fR) != CUFFT_SUCCESS) {
          fprintf(stderr, "fmm3: error: cufft inverse failed\n");
          exit(1);
        }
        fext<<<grid((long)nc2 * lB3 * 3), nthread>>>(nc2, k0, lB, lst, fR, ul);
      }
    if (!ffm) {
      int K = lKc, no2 = 2 * K + 1;
      long Dy = lD.y + 2 * K, Dz = lD.z + 2 * K, c0 = ((long)K * Dy + K) * Dz + K;
      long ncol = ((long)lD.y - 1) * Dz + lD.z, sp = Dy * Dz;
      lpack<<<grid((lD.x + 2L * K) * Dy * Dz * 3 * lB3), nthread>>>(t, ll, lo[ll], lB3, K, lax, lD, dX);
      for (int o = 0; o < nof; o++) {
        int oo[3] = {o / (no2 * no2) - K, (o / no2) % no2 - K, o % no2 - K}, a[3] = {lax.x, lax.y, lax.z};
        long sh = ((long)oo[a[0]] * Dy + oo[a[1]]) * Dz + oo[a[2]];
        double one = 1, beta = o > 0;
        for (int c = 0; c < 3; c++)
          for (int k = 0; k < 2; k++)
            if (cublasDgemmStridedBatched(hb, CUBLAS_OP_N, CUBLAS_OP_N, lB3, ncol, lB3, &one, dG + ((size_t)o * 3 + c) * lB3 * lB3, lB3, 0, dX + ((c0 + sh) * 3 + (c + 1 + k) % 3) * lB3, 3 * lB3, sp * 3 * lB3, &beta, dY + (c0 * 6 + 2 * c + k) * lB3, 6 * lB3, sp * 6 * lB3, lD.x) != CUBLAS_STATUS_SUCCESS) {
              fprintf(stderr, "fmm3: error: cublasDgemm failed\n");
              exit(1);
            }
      }
      lcombp<<<grid((long)nlv * lB3), nthread>>>(t, ll, lo[ll], lB3, K, lax, lD, dY, ul);
    }
    ck(cudaEventRecord(e1), "cudaEventRecord");
    ck(cudaEventSynchronize(e1), "cudaEventSynchronize");
    float lm;
    ck(cudaEventElapsedTime(&lm, e0, e1), "cudaEventElapsedTime");
    fprintf(stderr, "fmm3: lattice near field B %d Kc %d leaves %d offsets %d: %.6e s\n", lB, lKc, nlv, nof, lm / 1000);
    launched("lattice");
  }
  if (t.per) {
    double qh[15], *qd;
    memset(qh, 0, sizeof(qh));
    for (long i = 0; i < N; i++)
      for (int c = 0; c < t.nc; c++) {
        double qq = host[nin * i + 3 + c], *xi = host + nin * i;
        qh[5 * c] += qq;
        qh[5 * c + 1] += qq * xi[0];
        qh[5 * c + 2] += qq * xi[1];
        qh[5 * c + 3] += qq * xi[2];
        qh[5 * c + 4] += qq * (xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
      }
    ck(cudaMalloc(&qd, sizeof(qh)), "cudaMalloc");
    ck(cudaMemcpy(qd, qh, sizeof(qh), cudaMemcpyHostToDevice), "cudaMemcpy");
    quad<<<grid(N), nthread>>>(t, qd);
    launched("quad");
  }
  ck(cudaEventRecord(ev[8]), "cudaEventRecord");
  double *ph = (double *)pmalloc(t.nc * N * sizeof(double)), *gr = (double *)pmalloc(3 * t.nc * N * sizeof(double));
  long *ph_ = (long *)pmalloc(N * sizeof(long));
  if (ph == NULL || gr == NULL || ph_ == NULL) {
    fprintf(stderr, "fmm3: error: malloc failed\n");
    exit(1);
  }
  ck(cudaMemcpy(ph, t.phi, t.nc * N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
  ck(cudaMemcpy(gr, t.gr, 3 * t.nc * N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
  ck(cudaMemcpy(ph_, perm, N * sizeof(long), cudaMemcpyDeviceToHost), "cudaMemcpy");
  double *ulv = NULL;
  if (ul != NULL) {
    ulv = (double *)pmalloc(3 * N * sizeof(double));
    ck(cudaMemcpy(ulv, ul, 3 * N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
  }
  int no = biot ? 3 : 4;
  for (long i = 0; i < N; i++) {
    long j = ph_[i];
    if (biot) {
      const double *g = gr + 9 * i;
      out[3 * j] = (g[7] - g[5]) / (4 * pi);
      out[3 * j + 1] = (g[2] - g[6]) / (4 * pi);
      out[3 * j + 2] = (g[3] - g[1]) / (4 * pi);
      if (ulv != NULL)
        for (int d = 0; d < 3; d++)
          out[3 * j + d] += ulv[3 * i + d];
    } else {
      out[4 * j] = ph[i];
      out[4 * j + 1] = gr[3 * i];
      out[4 * j + 2] = gr[3 * i + 1];
      out[4 * j + 3] = gr[3 * i + 2];
    }
  }
  ck(cudaEventRecord(ev[9]), "cudaEventRecord");
  ck(cudaEventSynchronize(ev[9]), "cudaEventSynchronize");
  if ((file = fopen(opath, "w")) == NULL) {
    fprintf(stderr, "fmm3: error: fail to open '%s'\n", opath);
    exit(1);
  }
  if (fwrite(out, sizeof(double), no * N, file) != (size_t)(no * N)) {
    fprintf(stderr, "fmm3: error: fail to write '%s'\n", opath);
    exit(1);
  }
  if (ferror(file) || fclose(file) != 0) {
    fprintf(stderr, "fmm3: error: fail to write '%s'\n", opath);
    exit(1);
  }
  int *cnt = (int *)malloc(nb * sizeof(int)), *chh = (int *)malloc(nb * sizeof(int)), *uoh = (int *)malloc((nb + 1) * sizeof(int)), *ulh = (int *)malloc((nU + 1) * sizeof(int)), ne = 0, nlf = 0;
  long npair = 0;
  ck(cudaMemcpy(cnt, t.count, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
  ck(cudaMemcpy(chh, t.ch, nb * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
  ck(cudaMemcpy(uoh, t.uo, (nb + 1) * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
  ck(cudaMemcpy(ulh, t.ul, nU * sizeof(int), cudaMemcpyDeviceToHost), "cudaMemcpy");
  for (int b = 0; b < nb; b++) {
    ne += cnt[b] > 0;
    nlf += chh[b] < 0 && cnt[b] > 0;
    if (chh[b] < 0)
      for (int e = uoh[b]; e < uoh[b + 1]; e++)
        npair += (long)cnt[b] * cnt[ulh[e] / 125];
  }
  if (timing) {
    const char *nm[] = {"sort", "tree", "tables", "up", "m2l", "down", "l2p_m2p", "p2p", "out"};
    printf("N p s boxes nonempty leaves depth U V W X p2p_pairs");
    for (int i = 0; i < 9; i++)
      printf(" %s", nm[i]);
    printf(" total\n%ld %d %d %d %d %d %d %d %d %d %d %ld", N, t.p, t.s, nb, ne, nlf, nl, nU, nV, nW, nX, npair);
    for (int i = 0; i < 9; i++) {
      ck(cudaEventElapsedTime(&ms, ev[i], ev[i + 1]), "cudaEventElapsedTime");
      printf(" %.6e", ms / 1000);
    }
    ck(cudaEventElapsedTime(&ms, ev[0], ev[9]), "cudaEventElapsedTime");
    printf(" %.6e\n", ms / 1000);
  }
  if (nchk > 0) {
    double *dref, *ref, e2 = 0, n2 = 0, ep = 0, pm = 0, eg = 0, gm = 0;
    if (nchk > N)
      nchk = N;
    ck(cudaMalloc(&dref, 12 * nchk * sizeof(double)), "cudaMalloc");
    if ((ref = (double *)malloc(12 * nchk * sizeof(double))) == NULL) {
      fprintf(stderr, "fmm3: error: malloc failed\n");
      exit(1);
    }
    check<<<nchk, TB>>>(t, pos, dref);
    launched("check");
    ck(cudaMemcpy(ref, dref, 12 * nchk * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    for (long j = 0; j < nchk; j++) {
      const double *r = ref + 12 * j;
      if (biot) {
        double u[3] = {(r[3 + 7] - r[3 + 5]) / (4 * pi), (r[3 + 2] - r[3 + 6]) / (4 * pi), (r[3 + 3] - r[3 + 1]) / (4 * pi)};
        for (int k = 0; k < 3; k++) {
          eg = fmax(eg, fabs(u[k] - out[3 * j + k]));
          gm = fmax(gm, fabs(u[k]));
          e2 += (u[k] - out[3 * j + k]) * (u[k] - out[3 * j + k]);
          n2 += u[k] * u[k];
        }
      } else {
        e2 += (r[0] - out[4 * j]) * (r[0] - out[4 * j]);
        n2 += r[0] * r[0];
        ep = fmax(ep, fabs(r[0] - out[4 * j]));
        pm = fmax(pm, fabs(r[0]));
        for (int k = 0; k < 3; k++) {
          eg = fmax(eg, fabs(r[3 + k] - out[4 * j + 1 + k]));
          gm = fmax(gm, fabs(r[3 + k]));
        }
      }
    }
    if (biot)
      printf("check first %ld: velocity relL2 %.16e maxrel %.16e\n", nchk, sqrt(e2 / n2), eg / gm);
    else
      printf("check first %ld: potential relL2 %.16e maxrel %.16e gradient maxrel %.16e\n", nchk, sqrt(e2 / n2), ep / pm, eg / gm);
  }
  return 0;
}
