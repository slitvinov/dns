#include <cufft.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

enum { nthread = 256, nblock = 1024 };
static const double pi = 3.141592653589793238;
static const double KT[3][4][8] = {
    {{1, 0, -5. / 2, 3. / 2}, {0, -1. / 2, 1, -1. / 2}},
    {{1, 0, -5. / 4, -35. / 12, 21. / 4, -25. / 12},
     {0, -2. / 3, 2. / 3, 13. / 8, -8. / 3, 25. / 24},
     {0, 1. / 12, -1. / 24, -3. / 8, 13. / 24, -5. / 24}},
    {{1, 0, -49. / 36, 0, -959. / 144, 2569. / 144, -727. / 48, 623. / 144},
     {0, -3. / 4, 3. / 4, 13. / 48, 95. / 24, -429. / 40, 1091. / 120, -623. / 240},
     {0, 3. / 20, -3. / 40, -1. / 6, -191. / 144, 2581. / 720, -2183. / 720, 623. / 720},
     {0, -1. / 60, 1. / 180, 1. / 48, 7. / 36, -37. / 72, 13. / 30, -89. / 720}}};
static const int KS[3] = {2, 3, 4};
#define K8 KT[2]
__constant__ double cK[4][8];
__constant__ int cS;
struct G {
  long nx, ny, nz, nzf, nf, N;
  double ax, ay, az, hx, hy, hz, s, nu, fp;
  int f;
};
static void cuda(cudaError_t e, const char *what) {
  if (e != cudaSuccess) {
    fprintf(stderr, "vm3: error: %s: %s\n", what, cudaGetErrorString(e));
    exit(1);
  }
}
static void cufft(cufftResult r, const char *what) {
  if (r != CUFFT_SUCCESS) {
    fprintf(stderr, "vm3: error: %s: cufft error %d\n", what, (int)r);
    exit(1);
  }
}
static void launched(const char *what) { cuda(cudaGetLastError(), what); }
static unsigned grid(long n) { return (n + nthread - 1) / nthread; }
static double *dalloc(long n) {
  double *a;
  cuda(cudaMalloc(&a, n * sizeof(double)), "cudaMalloc");
  return a;
}
static double now(void) {
  struct timeval v;
  gettimeofday(&v, NULL);
  return v.tv_sec + 1e-6 * v.tv_usec;
}
__host__ __device__ static double kw(const double (*c)[8], int S, double x) {
  double a = fabs(x);
  if (a >= S)
    return 0;
  int j = (int)a;
  double u = a - j, w = c[j][7];
  for (int e = 6; e >= 0; e--)
    w = w * u + c[j][e];
  return w;
}
__device__ static double W(double x) { return kw(cK, cS, x); }
__host__ __device__ static long wrap(long i, long n) {
  i %= n;
  return i < 0 ? i + n : i;
}
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
__device__ static int nyq(G g, long i, long j, long k) {
  return 2 * i == g.nx || 2 * j == g.ny || 2 * k == g.nz;
}
__device__ static double p36(double x) {
  double x2 = x * x, x4 = x2 * x2, x8 = x4 * x4, x16 = x8 * x8;
  return x16 * x16 * x4;
}
__device__ static double mask(G g, long i, long j, long k) {
  if (g.f == 2) {
    double x = 2.0 * labs(wave(i, g.nx)) / g.nx, y = 2.0 * labs(wave(j, g.ny)) / g.ny,
           z = 2.0 * k / g.nz;
    return exp(-36 * (p36(x) + p36(y) + p36(z)));
  }
  if (g.f == 1)
    return 3 * labs(wave(i, g.nx)) < g.nx && 3 * labs(wave(j, g.ny)) < g.ny &&
           3 * k < g.nz;
  return 1;
}
__global__ static void filt(G g, double2 *o) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.nf;
       l += (long)gridDim.x * blockDim.x) {
    long i, j, k;
    double kx, ky, kz;
    mode(g, l, &i, &j, &k, &kx, &ky, &kz);
    double m = mask(g, i, j, k) * g.s;
    o[l].x *= m;
    o[l].y *= m;
  }
}
__global__ static void spec(G g, int m, int a, int c, const double2 *w0,
                            const double2 *w1, const double2 *w2, double2 *o) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.nf;
       l += (long)gridDim.x * blockDim.x) {
    long i, j, k;
    double kx, ky, kz;
    mode(g, l, &i, &j, &k, &kx, &ky, &kz);
    double kk = kx * kx + ky * ky + kz * kz, q[3] = {kx, ky, kz};
    double2 w[3] = {w0[l], w1[l], w2[l]}, r = {0, 0};
    double f = mask(g, i, j, k);
    for (int c2 = 0; c2 < 3; c2++) {
      w[c2].x *= f;
      w[c2].y *= f;
    }
    if (m == 2) {
      r.x = -g.nu * kk * g.s * w[c].x;
      r.y = -g.nu * kk * g.s * w[c].y;
    } else if (m == 3) {
      if (!nyq(g, i, j, k)) {
        r.x = -q[a] * w[c].y * g.s;
        r.y = q[a] * w[c].x * g.s;
      }
    } else if (kk > 0 && !nyq(g, i, j, k)) {
      int p = (c + 1) % 3, n = (c + 2) % 3;
      double2 x = {q[p] * w[n].x - q[n] * w[p].x, q[p] * w[n].y - q[n] * w[p].y};
      double2 u = {-x.y * g.s / kk, x.x * g.s / kk};
      if (m == 0)
        r = u;
      else {
        r.x = -q[a] * u.y;
        r.y = q[a] * u.x;
      }
    }
    o[l] = r;
  }
}
__global__ static void curl(G g, int c, const double2 *u0, const double2 *u1,
                            const double2 *u2, double2 *o) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.nf;
       l += (long)gridDim.x * blockDim.x) {
    long i, j, k;
    double kx, ky, kz;
    mode(g, l, &i, &j, &k, &kx, &ky, &kz);
    double q[3] = {kx, ky, kz};
    double2 u[3] = {u0[l], u1[l], u2[l]}, r = {0, 0};
    if (!nyq(g, i, j, k)) {
      int p = (c + 1) % 3, n = (c + 2) % 3;
      double2 x = {q[p] * u[n].x - q[n] * u[p].x, q[p] * u[n].y - q[n] * u[p].y};
      r.x = -x.y * g.s;
      r.y = x.x * g.s;
    }
    o[l] = r;
  }
}
__global__ static void proj(G g, int c, int pm, const double2 *w0, const double2 *w1,
                            const double2 *w2, double2 *o) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.nf;
       l += (long)gridDim.x * blockDim.x) {
    long i, j, k;
    double kx, ky, kz;
    mode(g, l, &i, &j, &k, &kx, &ky, &kz);
    double kk = kx * kx + ky * ky + kz * kz, q[3] = {kx, ky, kz};
    double2 w[3] = {w0[l], w1[l], w2[l]}, r = w[c];
    if (pm < 0) {
      double f = mask(g, i, j, k);
      r.x *= f * g.s;
      r.y *= f * g.s;
      o[l] = r;
      continue;
    }
    if (kk > 0) {
      double2 d = {kx * w[0].x + ky * w[1].x + kz * w[2].x,
                   kx * w[0].y + ky * w[1].y + kz * w[2].y};
      r.x -= q[c] * d.x / kk;
      r.y -= q[c] * d.y / kk;
    }
    double f = pm == 3 ? (g.fp == 1 ? mask(g, i, j, k) : pow(mask(g, i, j, k), g.fp)) : pm ? mask(g, i, j, k) : 1;
    r.x *= f * g.s;
    r.y *= f * g.s;
    o[l] = r;
  }
}
__global__ static void addc(long n, const double *a, double c, double *o) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < n;
       l += (long)gridDim.x * blockDim.x)
    o[l] = a[l] + c;
}
__global__ static void axpy(long n, double c, double *r, const double *a, const double *b) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < n;
       l += (long)gridDim.x * blockDim.x)
    r[l] += c * a[l] * b[l];
}
__global__ static void lin(long n, const double *a, double c, const double *b, double *o) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < n;
       l += (long)gridDim.x * blockDim.x)
    o[l] = a[l] + c * b[l];
}
__device__ static void line(G g, int d, long l, long *b, long *st, long *n, long *i) {
  long x = l / (g.ny * g.nz), y = l / g.nz % g.ny, z = l % g.nz;
  if (d == 0) {
    *n = g.nx, *st = g.ny * g.nz, *i = x, *b = l - x * *st;
  } else if (d == 1) {
    *n = g.ny, *st = g.nz, *i = y, *b = l - y * *st;
  } else {
    *n = g.nz, *st = 1, *i = z, *b = l - z;
  }
}
__device__ static double ip1(const double *f, long b, long st, long n, double x) {
  long c = (long)floor(x);
  double v = 0, t = 0;
  for (int o = 0; o < 2 * cS; o++) {
    long q = c - cS + 1 + o;
    double w = W(q - x);
    t += w;
    v += w * f[b + wrap(q, n) * st];
  }
  return v / t;
}
__global__ static void inv1(G g, int d, const double *e, const double *f0, const double *f1,
                            double *o0, double *o1, double *bad) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.N;
       l += (long)gridDim.x * blockDim.x) {
    long b, st, n, i;
    line(g, d, l, &b, &st, &n, &i);
    double off = -e[l], prev = 1e300;
    int it;
    for (it = 0; it < 60 && fabs(off - prev) > 1e-14 * (1 + fabs(off)); it++) {
      prev = off;
      off = -ip1(e, b, st, n, i + off);
    }
    bad[l] = it;
    o0[l] = ip1(f0, b, st, n, i + off);
    if (o1 != NULL)
      o1[l] = ip1(f1, b, st, n, i + off);
  }
}
__global__ static void push1(G g, int d, long D, const double *e, const double *a0,
                             const double *a1, const double *a2, double *o0, double *o1,
                             double *o2) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.N;
       l += (long)gridDim.x * blockDim.x) {
    long b, st, n, i;
    line(g, d, l, &b, &st, &n, &i);
    double s0 = 0, s1 = 0, s2 = 0;
    for (long p = -cS + 1 - D; p <= cS - 1 + D; p++) {
      long m = b + wrap(i - p, n) * st;
      double x = p - e[m], w = W(x);
      if (w == 0)
        continue;
      long c = (long)floor(e[m]);
      double t = 0;
      for (int o = 0; o < 2 * cS; o++)
        t += W(c - cS + 1 + o - e[m]);
      w /= t;
      s0 += w * a0[m];
      if (a1 != NULL) {
        s1 += w * a1[m];
        s2 += w * a2[m];
      }
    }
    o0[l] = s0;
    if (o1 != NULL) {
      o1[l] = s1;
      o2[l] = s2;
    }
  }
}
__global__ static void pull1(G g, int d, const double *e, const double *f, double *o) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.N;
       l += (long)gridDim.x * blockDim.x) {
    long b, st, n, i;
    line(g, d, l, &b, &st, &n, &i);
    o[l] = ip1(f, b, st, n, i + e[l]);
  }
}
__global__ static void tint(G g, const double *ex, const double *ey, const double *ez,
                            const double *f, double *o) {
  for (long m = blockIdx.x * (long)blockDim.x + threadIdx.x; m < g.N;
       m += (long)gridDim.x * blockDim.x) {
    long i = m / (g.ny * g.nz), j = m / g.nz % g.ny, k = m % g.nz;
    double e[3] = {ex[m], ey[m], ez[m]}, w[3][8], v = 0;
    long b[3];
    for (int d = 0; d < 3; d++) {
      double t = 0;
      b[d] = (long)floor(e[d]);
      for (int q = 0; q < 2 * cS; q++)
        t += w[d][q] = W(b[d] - cS + 1 + q - e[d]);
      for (int q = 0; q < 2 * cS; q++)
        w[d][q] /= t;
    }
    for (int p = 0; p < 2 * cS; p++)
      for (int q = 0; q < 2 * cS; q++)
        for (int r = 0; r < 2 * cS; r++)
          v += w[0][p] * w[1][q] * w[2][r] *
               f[(wrap(i + b[0] - cS + 1 + p, g.nx) * g.ny + wrap(j + b[1] - cS + 1 + q, g.ny)) *
                     g.nz +
                 wrap(k + b[2] - cS + 1 + r, g.nz)];
    o[m] = v;
  }
}
__global__ static void norm(G g, const double *ex, const double *ey, const double *ez,
                           double *f) {
  for (long m = blockIdx.x * (long)blockDim.x + threadIdx.x; m < g.N;
       m += (long)gridDim.x * blockDim.x) {
    double e[3] = {ex[m], ey[m], ez[m]}, q = 1;
    for (int d = 0; d < 3; d++) {
      long b = (long)floor(e[d]);
      double t = 0;
      for (int o = 0; o < 2 * cS; o++)
        t += W(b - cS + 1 + o - e[d]);
      q *= t;
    }
    f[m] = 1 / q;
  }
}
__global__ static void remesh(G g, long Dx, long Dy, long Dz, const double *nf,
                              const double *ex, const double *ey, const double *ez,
                              const double *a0,
                              const double *a1, const double *a2, double *o0,
                              double *o1, double *o2) {
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.N;
       l += (long)gridDim.x * blockDim.x) {
    long i = l / (g.ny * g.nz), j = l / g.nz % g.ny, k = l % g.nz;
    double s0 = 0, s1 = 0, s2 = 0;
    long kr[23];
    long h = cS - 1;
    for (long r = -h - Dz; r <= h + Dz; r++)
      kr[r + h + Dz] = wrap(k - r, g.nz);
    for (long p = -h - Dx; p <= h + Dx; p++) {
      long mi = wrap(i - p, g.nx) * g.ny;
      for (long q = -h - Dy; q <= h + Dy; q++) {
        long b = (mi + wrap(j - q, g.ny)) * g.nz;
        for (long r = -h - Dz; r <= h + Dz; r++) {
          long m = b + kr[r + h + Dz];
          double wx = W(p - ex[m]);
          if (wx == 0)
            continue;
          double wy = W(q - ey[m]);
          if (wy == 0)
            continue;
          double w = wx * wy * W(r - ez[m]) * nf[m];
          s0 += w * a0[m];
          s1 += w * a1[m];
          s2 += w * a2[m];
        }
      }
    }
    o0[l] = s0;
    o1[l] = s1;
    o2[l] = s2;
  }
}
__global__ static void scatter(G g, const double *nf, const double *ex, const double *ey,
                               const double *ez, const double *a0, const double *a1,
                               const double *a2, double *o0, double *o1, double *o2) {
  for (long m = blockIdx.x * (long)blockDim.x + threadIdx.x; m < g.N;
       m += (long)gridDim.x * blockDim.x) {
    long i = m / (g.ny * g.nz), j = m / g.nz % g.ny, k = m % g.nz;
    long bi = (long)floor(ex[m]), bj = (long)floor(ey[m]), bk = (long)floor(ez[m]);
    for (long p = bi - cS + 1; p <= bi + cS; p++)
      for (long q = bj - cS + 1; q <= bj + cS; q++)
        for (long r = bk - cS + 1; r <= bk + cS; r++) {
          double w = W(p - ex[m]) * W(q - ey[m]) * W(r - ez[m]) * nf[m];
          long n = (wrap(i + p, g.nx) * g.ny + wrap(j + q, g.ny)) * g.nz + wrap(k + r, g.nz);
          atomicAdd(&o0[n], w * a0[m]);
          atomicAdd(&o1[n], w * a1[m]);
          atomicAdd(&o2[n], w * a2[m]);
        }
  }
}
__global__ static void stage(G g, int first, int pre, double hb, double hc, double dt,
                             const double *f0, const double *f1, const double *f2,
                             const double *f3, const double *f4, const double *f5,
                             double *ex, double *ey, double *ez, double *a0,
                             double *a1, double *a2, const double *b0, const double *b1,
                             const double *b2, double *ae0, double *ae1, double *ae2,
                             double *aw0, double *aw1, double *aw2) {
  for (long m = blockIdx.x * (long)blockDim.x + threadIdx.x; m < g.N;
       m += (long)gridDim.x * blockDim.x) {
    double v[6];
    if (first || pre) {
      v[0] = f0[m];
      v[1] = f1[m];
      v[2] = f2[m];
      v[3] = f3[m];
      v[4] = f4[m];
      v[5] = f5[m];
    } else {
      long i = m / (g.ny * g.nz), j = m / g.nz % g.ny, k = m % g.nz;
      double e[3] = {ex[m], ey[m], ez[m]}, w[3][8];
      long b[3];
      for (int d = 0; d < 3; d++) {
        double t = 0;
        b[d] = (long)floor(e[d]);
        for (int o = 0; o < 2 * cS; o++)
          t += w[d][o] = W(b[d] - cS + 1 + o - e[d]);
        for (int o = 0; o < 2 * cS; o++)
          w[d][o] /= t;
      }
      for (int c = 0; c < 6; c++)
        v[c] = 0;
      for (int p = 0; p < 2 * cS; p++) {
        long ni = wrap(i + b[0] - cS + 1 + p, g.nx);
        for (int q = 0; q < 2 * cS; q++) {
          long nj = wrap(j + b[1] - cS + 1 + q, g.ny);
          double wpq = w[0][p] * w[1][q];
          for (int r = 0; r < 2 * cS; r++) {
            long n = (ni * g.ny + nj) * g.nz + wrap(k + b[2] - cS + 1 + r, g.nz);
            double x = wpq * w[2][r];
            v[0] += x * f0[n];
            v[1] += x * f1[n];
            v[2] += x * f2[n];
            v[3] += x * f3[n];
            v[4] += x * f4[n];
            v[5] += x * f5[n];
          }
        }
      }
    }
    double ke[3] = {v[0] / g.hx, v[1] / g.hy, v[2] / g.hz};
    if (first) {
      ae0[m] = hb * ke[0];
      ae1[m] = hb * ke[1];
      ae2[m] = hb * ke[2];
      aw0[m] = hb * v[3];
      aw1[m] = hb * v[4];
      aw2[m] = hb * v[5];
    } else {
      ae0[m] += hb * ke[0];
      ae1[m] += hb * ke[1];
      ae2[m] += hb * ke[2];
      aw0[m] += hb * v[3];
      aw1[m] += hb * v[4];
      aw2[m] += hb * v[5];
    }
    if (hc > 0) {
      ex[m] = hc * dt * ke[0];
      ey[m] = hc * dt * ke[1];
      ez[m] = hc * dt * ke[2];
      a0[m] = b0[m] + hc * dt * v[3];
      a1[m] = b1[m] + hc * dt * v[4];
      a2[m] = b2[m] + hc * dt * v[5];
    } else {
      ex[m] = dt * ae0[m];
      ey[m] = dt * ae1[m];
      ez[m] = dt * ae2[m];
      a0[m] = b0[m] + dt * aw0[m];
      a1[m] = b1[m] + dt * aw1[m];
      a2[m] = b2[m] + dt * aw2[m];
    }
  }
}
__host__ __device__ static void add2(double *s, double *c, double x) {
  double t = *s + x, z = t - *s;
  *c += (*s - (t - z)) + (x - z);
  *s = t;
}
__device__ static void blockreduce2(int K, const double *v, const double *w, double *s,
                                    double *out) {
  for (int c = 0; c < K; c++) {
    s[2 * c * nthread + threadIdx.x] = v[c];
    s[(2 * c + 1) * nthread + threadIdx.x] = w[c];
  }
  __syncthreads();
  for (int h = nthread / 2; h > 0; h /= 2) {
    if ((int)threadIdx.x < h)
      for (int c = 0; c < K; c++) {
        double *a = &s[2 * c * nthread + threadIdx.x], *b = &s[(2 * c + 1) * nthread + threadIdx.x];
        double x = a[h], y = b[h];
        add2(a, b, x);
        *b += y;
      }
    __syncthreads();
  }
  if (threadIdx.x == 0)
    for (int c = 0; c < K; c++) {
      out[2 * c * gridDim.x + blockIdx.x] = s[2 * c * nthread];
      out[(2 * c + 1) * gridDim.x + blockIdx.x] = s[(2 * c + 1) * nthread];
    }
}
__device__ static void blockreduce(int K, int op, const double *v, double *s,
                                   double *out) {
  for (int c = 0; c < K; c++)
    s[c * nthread + threadIdx.x] = v[c];
  __syncthreads();
  for (int h = nthread / 2; h > 0; h /= 2) {
    if ((int)threadIdx.x < h)
      for (int c = 0; c < K; c++) {
        double *a = &s[c * nthread + threadIdx.x], b = a[h];
        *a = op ? fmax(*a, b) : *a + b;
      }
    __syncthreads();
  }
  if (threadIdx.x == 0)
    for (int c = 0; c < K; c++)
      out[c * gridDim.x + blockIdx.x] = s[c * nthread];
}
__global__ static void rsum(long n, const double *a, const double *b, const double *c,
                            double *out) {
  __shared__ double s[8 * nthread];
  double e[4] = {0, 0, 0, 0}, k[4] = {0, 0, 0, 0};
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < n;
       l += (long)gridDim.x * blockDim.x) {
    add2(&e[0], &k[0], a[l] * a[l]);
    add2(&e[0], &k[0], b[l] * b[l]);
    add2(&e[0], &k[0], c[l] * c[l]);
    add2(&e[1], &k[1], a[l]);
    add2(&e[2], &k[2], b[l]);
    add2(&e[3], &k[3], c[l]);
  }
  blockreduce2(4, e, k, s, out);
}
__global__ static void rmax(long n, const double *a, const double *b, const double *c,
                           double *out) {
  __shared__ double s[4 * nthread];
  double e[4] = {0, 0, 0, 0};
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < n;
       l += (long)gridDim.x * blockDim.x) {
    e[0] = fmax(e[0], sqrt(a[l] * a[l] + b[l] * b[l] + c[l] * c[l]));
    e[1] = fmax(e[1], fabs(a[l]));
    e[2] = fmax(e[2], fabs(b[l]));
    e[3] = fmax(e[3], fabs(c[l]));
  }
  blockreduce(4, 1, e, s, out);
}
__global__ static void rimp(G g, const double *w0, const double *w1, const double *w2,
                            const double *u0, const double *u1, const double *u2,
                            double *out) {
  __shared__ double s[14 * nthread];
  double e[7] = {0, 0, 0, 0, 0, 0, 0}, k[7] = {0, 0, 0, 0, 0, 0, 0};
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.N;
       l += (long)gridDim.x * blockDim.x) {
    double x[3] = {(l / (g.ny * g.nz)) * g.hx, (l / g.nz % g.ny) * g.hy, (l % g.nz) * g.hz},
           w[3] = {w0[l], w1[l], w2[l]};
    double xw[3] = {x[1] * w[2] - x[2] * w[1], x[2] * w[0] - x[0] * w[2],
                    x[0] * w[1] - x[1] * w[0]};
    for (int c = 0; c < 3; c++) {
      int p = (c + 1) % 3, q = (c + 2) % 3;
      add2(&e[c], &k[c], 0.5 * xw[c]);
      add2(&e[3 + c], &k[3 + c], (x[p] * xw[q] - x[q] * xw[p]) / 3);
    }
    add2(&e[6], &k[6], u0[l] * w[0] + u1[l] * w[1] + u2[l] * w[2]);
  }
  blockreduce2(7, e, k, s, out);
}
__global__ static void rng(long n, const double *a, double *out) {
  __shared__ double s[2 * nthread];
  double e[2] = {0, -INFINITY};
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < n;
       l += (long)gridDim.x * blockDim.x) {
    double x = fabs(a[l]);
    e[0] = fmax(e[0], x);
    if (x > 0)
      e[1] = fmax(e[1], -x);
  }
  blockreduce(2, 1, e, s, out);
}
__global__ static void rspec(G g, const double2 *w0, const double2 *w1,
                             const double2 *w2, double *out) {
  __shared__ double s[8 * nthread];
  double e[4] = {0, 0, 0, 0}, kc[4] = {0, 0, 0, 0};
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < g.nf;
       l += (long)gridDim.x * blockDim.x) {
    long i, j, k;
    double kx, ky, kz;
    mode(g, l, &i, &j, &k, &kx, &ky, &kz);
    double kk = kx * kx + ky * ky + kz * kz, h = k == 0 || 2 * k == g.nz ? 0.5 : 1;
    double2 a = w0[l], b = w1[l], c = w2[l];
    double q = a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y + c.x * c.x + c.y * c.y;
    double dr = kx * a.x + ky * b.x + kz * c.x, di = kx * a.y + ky * b.y + kz * c.y;
    double dd = dr * dr + di * di;
    add2(&e[0], &kc[0], h * (kk > 0 ? q - dd / kk : 0));
    add2(&e[3], &kc[3], h * (kk > 0 ? (q - dd / kk) / kk : 0));
    add2(&e[1], &kc[1], h * dd);
    add2(&e[2], &kc[2], h * kk * q);
  }
  blockreduce2(4, e, kc, s, out);
}
static double *dP, *P;
static void collect(int K, int op, double *r) {
  cuda(cudaMemcpy(P, dP, K * nblock * sizeof(double), cudaMemcpyDeviceToHost),
       "cudaMemcpy");
  for (int c = 0; c < K; c++) {
    r[c] = 0;
    for (int b = 0; b < nblock; b++)
      r[c] = op ? fmax(r[c], P[c * nblock + b]) : r[c] + P[c * nblock + b];
  }
}
static void collect2(int K, double *r) {
  cuda(cudaMemcpy(P, dP, 2 * K * nblock * sizeof(double), cudaMemcpyDeviceToHost),
       "cudaMemcpy");
  for (int c = 0; c < K; c++) {
    double a = 0, b = 0;
    for (int i = 0; i < nblock; i++) {
      add2(&a, &b, P[2 * c * nblock + i]);
      b += P[(2 * c + 1) * nblock + i];
    }
    r[c] = a + b;
  }
}
static cufftHandle fplan, bplan;
static int verbose;
static void look(long n, const char *name, double *a) {
  double r[2];
  rng<<<nblock, nthread>>>(n, a, dP);
  launched("rng");
  cuda(cudaMemcpy(P, dP, 2 * nblock * sizeof(double), cudaMemcpyDeviceToHost),
       "cudaMemcpy");
  r[0] = 0;
  r[1] = -INFINITY;
  for (int b = 0; b < nblock; b++) {
    r[0] = fmax(r[0], P[b]);
    r[1] = fmax(r[1], P[nblock + b]);
  }
  fprintf(stderr, "vm3: look %-6s max|x| % .3e min|x|>0 % .3e log10 range %.1f\n", name,
          r[0], -r[1], r[0] > 0 && -r[1] > 0 ? log10(r[0] / -r[1]) : 0.0);
}
static void plan(cufftHandle *p, long nx, long ny, long nz, cufftType t) {
  long long n[3] = {nx, ny, nz};
  size_t w;
  cufft(cufftCreate(p), "cufftCreate");
  cufft(cufftMakePlanMany64(*p, 3, n, NULL, 1, 0, NULL, 1, 0, t, 1, &w),
        "cufftMakePlanMany64");
}
static double2 *oh[3], *sc;
static double *T;
static void fwd(double *a, double2 *o) {
  cufft(cufftExecD2Z(fplan, a, o), "cufftExecD2Z");
}
static void bwd(double2 *c, double *a) {
  cufft(cufftExecZ2D(bplan, c, a), "cufftExecZ2D");
}
static double gmax, *wu[3];
static int sform;
static void eval(G g, double **om, double **u, double **r, int diag, double **kE) {
  for (int c = 0; c < 3; c++)
    fwd(om[c], oh[c]);
  if (sform == 2 || g.f)
    for (int c = 0; c < 3; c++) {
      proj<<<grid(g.nf), nthread>>>(g, c, sform == 2 ? 1 : -1, oh[0], oh[1], oh[2], sc);
      launched("proj");
      bwd(sc, wu[c]);
    }
  for (int c = 0; c < 3; c++) {
    spec<<<grid(g.nf), nthread>>>(g, 0, 0, c, oh[0], oh[1], oh[2], sc);
    launched("spec");
    bwd(sc, u[c]);
    spec<<<grid(g.nf), nthread>>>(g, 2, 0, c, oh[0], oh[1], oh[2], sc);
    launched("spec");
    bwd(sc, r[c]);
  }
  gmax = 0;
  for (int a = 0; a < 3; a++)
    for (int c = 0; c < 3; c++) {
      spec<<<grid(g.nf), nthread>>>(g, 1, a, c, oh[0], oh[1], oh[2], sc);
      launched("spec");
      bwd(sc, T);
      double **v = sform == 2 || g.f ? wu : om;
      if (sform == 1 || sform == 3)
        axpy<<<grid(g.N), nthread>>>(g.N, sform == 3 ? 0.5 : 1, r[c], v[a], T);
      if (sform != 1)
        axpy<<<grid(g.N), nthread>>>(g.N, sform == 3 ? 0.5 : 1, r[a], v[c], T);
      launched("axpy");
      if (diag) {
        double m[4];
        rmax<<<nblock, nthread>>>(g.N, T, T, T, dP);
        launched("rmax");
        collect(4, 1, m);
        gmax = fmax(gmax, m[1]);
      }
    }
  if (g.f)
    for (int c = 0; c < 3; c++) {
      fwd(r[c], sc);
      filt<<<grid(g.nf), nthread>>>(g, sc);
      launched("filt");
      bwd(sc, r[c]);
    }
  if (kE == NULL)
    return;
  for (int c = 0; c < 3; c++) {
    cuda(cudaMemset(kE[c], 0, g.N * sizeof(double)), "cudaMemset");
    for (int a = 0; a < 3; a++) {
      spec<<<grid(g.nf), nthread>>>(g, 3, a, c, oh[0], oh[1], oh[2], sc);
      launched("spec");
      bwd(sc, T);
      axpy<<<grid(g.N), nthread>>>(g.N, -1, kE[c], u[a], T);
      launched("axpy");
    }
    if (g.f) {
      fwd(kE[c], sc);
      filt<<<grid(g.nf), nthread>>>(g, sc);
      launched("filt");
      bwd(sc, kE[c]);
    }
    lin<<<grid(g.N), nthread>>>(g.N, kE[c], 1, r[c], kE[c]);
    launched("lin");
  }
}
static void maxdisp(G g, double *e0, double *e1, double *e2, long cap, long *D) {
  double m[4];
  rmax<<<nblock, nthread>>>(g.N, e0, e1, e2, dP);
  launched("rmax");
  collect(4, 1, m);
  for (int d = 0; d < 3; d++) {
    D[d] = (long)ceil(m[d + 1]);
    if (D[d] > cap) {
      fprintf(stderr, "vm3: error: displacement %g cells > %ld (reduce -s or -C)\n",
              m[d + 1], cap);
      exit(1);
    }
  }
}
static double *nf;
static void doremesh(G g, double **e, double **a, double **o) {
  long D[3];
  maxdisp(g, e[0], e[1], e[2], 8, D);
  norm<<<grid(g.N), nthread>>>(g, e[0], e[1], e[2], nf);
  launched("norm");
  remesh<<<grid(g.N), nthread>>>(g, D[0], D[1], D[2], nf, e[0], e[1], e[2], a[0], a[1],
                                 a[2], o[0], o[1], o[2]);
  launched("remesh");
}
static double *fy, *fz1, *fz, *T2, *itn;
static int itmax;
static void factor(G g, double **e) {
  double m[4];
  inv1<<<grid(g.N), nthread>>>(g, 0, e[0], e[1], e[2], fy, fz1, itn);
  launched("inv1");
  rmax<<<nblock, nthread>>>(g.N, itn, itn, itn, dP);
  launched("rmax");
  collect(4, 1, m);
  itmax = m[1];
  inv1<<<grid(g.N), nthread>>>(g, 1, fy, fz1, NULL, fz, NULL, itn);
  launched("inv1");
  rmax<<<nblock, nthread>>>(g.N, itn, itn, itn, dP);
  launched("rmax");
  collect(4, 1, m);
  itmax = fmax(itmax, m[1]);
  if (itmax >= 60) {
    fprintf(stderr, "vm3: error: map inversion did not converge (dt |grad u| too large)\n");
    exit(1);
  }
}
static void splitremesh(G g, double **e, double **a, double **o) {
  long D[3];
  double *t[3] = {fz1, itn, T2};
  maxdisp(g, e[0], fy, fz, 64, D);
  push1<<<grid(g.N), nthread>>>(g, 0, D[0], e[0], a[0], a[1], a[2], o[0], o[1], o[2]);
  launched("push1");
  push1<<<grid(g.N), nthread>>>(g, 1, D[1], fy, o[0], o[1], o[2], t[0], t[1], t[2]);
  launched("push1");
  push1<<<grid(g.N), nthread>>>(g, 2, D[2], fz, t[0], t[1], t[2], o[0], o[1], o[2]);
  launched("push1");
}
static void splitinterp(G g, double **e, double **f, int n) {
  for (int c = 0; c < n; c++) {
    pull1<<<grid(g.N), nthread>>>(g, 2, fz, f[c], T2);
    launched("pull1");
    pull1<<<grid(g.N), nthread>>>(g, 1, fy, T2, itn);
    launched("pull1");
    pull1<<<grid(g.N), nthread>>>(g, 0, e[0], itn, f[c]);
    launched("pull1");
  }
}
static void project(G g, double **om) {
  for (int c = 0; c < 3; c++)
    fwd(om[c], oh[c]);
  for (int c = 0; c < 3; c++) {
    proj<<<grid(g.nf), nthread>>>(g, c, g.f ? 3 : 0, oh[0], oh[1], oh[2], sc);
    launched("proj");
    bwd(sc, om[c]);
  }
}
static double hrand(unsigned long *x) {
  *x ^= *x << 13;
  *x ^= *x >> 7;
  *x ^= *x << 17;
  return (*x >> 11) * (1.0 / 9007199254740992.0);
}
static double cpow(double x, int m) {
  double r = 1;
  while (m-- > 0)
    r *= x;
  return r;
}
static void tsplit(long n, double amp) {
  G g;
  g.nx = g.ny = g.nz = n;
  g.nzf = n / 2 + 1;
  g.nf = n * n * g.nzf;
  g.N = n * n * n;
  g.ax = g.ay = g.az = 1;
  g.hx = g.hy = g.hz = 2 * pi / n;
  g.s = 1.0 / g.N;
  g.nu = 0;
  g.f = 0;
  g.fp = 1;
  double *h[8], *d[16];
  for (int c = 0; c < 8; c++)
    h[c] = (double *)malloc(g.N * sizeof(double));
  for (int c = 0; c < 16; c++)
    d[c] = dalloc(g.N);
  for (long m = 0; m < g.N; m++) {
    long i = m / (n * n), j = m / n % n, k = m % n;
    double x = i * g.hx, y = j * g.hy, z = k * g.hz;
    h[0][m] = amp * sin(y + 0.3) * cos(z) / g.hx;
    h[1][m] = amp * sin(z + 0.7) * cos(x) / g.hy;
    h[2][m] = amp * sin(x + 1.1) * cos(y) / g.hz;
    h[3][m] = exp(sin(x)) * cos(y - z);
    h[4][m] = cos(x + 2 * y) * sin(z);
    h[5][m] = sin(x - z) + cos(2 * y);
    double X = x + amp * sin(y + 0.3) * cos(z), Y = y + amp * sin(z + 0.7) * cos(x),
           Z = z + amp * sin(x + 1.1) * cos(y);
    h[6][m] = exp(sin(X)) * cos(Y - Z);
  }
  for (int c = 0; c < 6; c++)
    cuda(cudaMemcpy(d[c], h[c], g.N * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
  double *e[3] = {d[0], d[1], d[2]}, *a[3] = {d[3], d[4], d[5]}, *o1[3] = {d[6], d[7], d[8]},
         *o2[3] = {d[9], d[10], d[11]};
  fy = d[12];
  fz1 = d[13];
  fz = d[14];
  itn = d[15];
  T2 = dalloc(g.N);
  double *nf0 = nf;
  nf = dalloc(g.N);
  doremesh(g, e, a, o1);
  factor(g, e);
  int it = itmax;
  splitremesh(g, e, a, o2);
  for (int c = 0; c < 3; c++) {
    cuda(cudaMemcpy(h[c], o1[c], g.N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    cuda(cudaMemcpy(h[3 + c], o2[c], g.N * sizeof(double), cudaMemcpyDeviceToHost),
         "cudaMemcpy");
  }
  double dm = 0, mm = 0, s1 = 0, k1 = 0, s2 = 0, k2 = 0;
  for (long m = 0; m < g.N; m++) {
    dm = fmax(dm, fabs(h[0][m] - h[3][m]));
    mm = fmax(mm, fabs(h[0][m]));
    add2(&s1, &k1, h[0][m]);
    add2(&s2, &k2, h[3][m]);
  }
  printf("split n %ld amp %g remesh max|split - tensor|/max %.3e sums tensor %.17g split %.17g it %d\n",
         n, amp, dm / mm, s1 + k1, s2 + k2, it);
  if (n == 32) {
    long D[3];
    double *hx = (double *)malloc(g.N * sizeof(double)), *ha = (double *)malloc(g.N * sizeof(double)),
           *hp = (double *)malloc(g.N * sizeof(double));
    unsigned long sd = 1234567UL;
    for (long m = 0; m < g.N; m++) {
      long i = m / (n * n);
      ha[m] = i >= 8 && i < n - 8 ? hrand(&sd) - 0.5 : 0;
    }
    cuda(cudaMemcpy(a[0], ha, g.N * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    maxdisp(g, e[0], e[0], e[0], 64, D);
    push1<<<grid(g.N), nthread>>>(g, 0, D[0], e[0], a[0], NULL, NULL, o2[0], NULL, NULL);
    launched("push1");
    cuda(cudaMemcpy(hx, e[0], g.N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    cuda(cudaMemcpy(hp, o2[0], g.N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
    for (int q = 0; q <= 7; q++) {
      double em = 0;
      for (long jk = 0; jk < n * n; jk++) {
        double P0 = 0, G0 = 0, S = 0;
        for (long i = 0; i < n; i++) {
          long m = i * n * n + jk;
          double x = i + hx[m] - n / 2;
          P0 += cpow(x, q) * ha[m];
          S += fabs(cpow(x, q) * ha[m]);
        }
        for (long i = 0; i < n; i++)
          G0 += cpow(i - n / 2.0, q) * hp[i * n * n + jk];
        em = fmax(em, fabs(P0 - G0) / S);
      }
      printf("split pass x line moment %d max rel %.3e\n", q, em);
    }
    free(hx);
    free(ha);
    free(hp);
  }
  for (long m = 0; m < g.N; m++) {
    long i = m / (n * n), j = m / n % n, k = m % n;
    double x = i * g.hx, y = j * g.hy, z = k * g.hz;
    h[3][m] = exp(sin(x)) * cos(y - z);
  }
  cuda(cudaMemcpy(d[3], h[3], g.N * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
  tint<<<grid(g.N), nthread>>>(g, e[0], e[1], e[2], d[3], d[6]);
  launched("tint");
  cuda(cudaMemcpy(h[0], d[6], g.N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
  double *f[1] = {d[3]};
  splitinterp(g, e, f, 1);
  cuda(cudaMemcpy(h[1], d[3], g.N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
  double et = 0, es = 0;
  for (long m = 0; m < g.N; m++) {
    et = fmax(et, fabs(h[0][m] - h[6][m]));
    es = fmax(es, fabs(h[1][m] - h[6][m]));
  }
  printf("split n %ld amp %g interp max err tensor %.3e split %.3e\n", n, amp, et, es);
  for (int c = 0; c < 16; c++)
    cudaFree(d[c]);
  cudaFree(T2);
  cudaFree(nf);
  nf = nf0;
  for (int c = 0; c < 8; c++)
    free(h[c]);
}
static void tests(void) {
  unsigned long seed = 88172645463325252UL;
  printf("test kernel m maxerr\n");
  for (int K = 0; K < 3; K++) {
    double mom[10];
    for (int m = 0; m < 10; m++)
      mom[m] = 0;
    for (int it = 0; it < 4001; it++) {
      double t = it / 4000.0;
      for (int m = 0; m < 10; m++) {
        double s = 0;
        for (int i = -4; i <= 4; i++)
          s += cpow(t - i, m) * kw(KT[K], KS[K], t - i);
        mom[m] = fmax(mom[m], fabs(s - (m == 0)));
      }
    }
    for (int m = 0; m <= 2 * KS[K]; m++)
      printf("moment m%d %d % .16e\n", 2 * KS[K], m, mom[m]);
  }
  for (int j = 1; j <= 4; j++)
    for (int q = 0; q <= 4; q++) {
      double a = 0, b = 0;
      for (int e = q; e < 8; e++) {
        double f = 1;
        for (int z = 0; z < q; z++)
          f *= e - z;
        a += f * K8[j - 1][e];
        if (j < 4 && e == q)
          b += f * K8[j][e];
      }
      printf("jump j %d deriv %d % .16e\n", j, q, a - b);
    }
  long n = 32;
  G g;
  g.nx = g.ny = g.nz = n;
  g.nzf = n / 2 + 1;
  g.nf = n * n * g.nzf;
  g.N = n * n * n;
  g.ax = g.ay = g.az = 1;
  g.hx = g.hy = g.hz = 2 * pi / n;
  g.s = 1.0 / g.N;
  g.nu = 0.01;
  g.f = 0;
  g.fp = 1;
  double *h[12], *d[12];
  for (int c = 0; c < 12; c++) {
    h[c] = (double *)malloc(g.N * sizeof(double));
    d[c] = dalloc(g.N);
  }
  for (long m = 0; m < g.N; m++) {
    long i = m / (n * n), j = m / n % n, k = m % n;
    int in = i >= 10 && i < 22 && j >= 10 && j < 22 && k >= 10 && k < 22;
    for (int c = 0; c < 3; c++)
      h[c][m] = in ? 2.4 * (hrand(&seed) - 0.5) : 0;
    for (int c = 3; c < 6; c++)
      h[c][m] = in ? hrand(&seed) - 0.5 : 0;
  }
  for (int c = 0; c < 6; c++)
    cuda(cudaMemcpy(d[c], h[c], g.N * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
  for (int c = 6; c < 12; c++)
    cuda(cudaMemset(d[c], 0, g.N * sizeof(double)), "cudaMemset");
  nf = dalloc(g.N);
  norm<<<grid(g.N), nthread>>>(g, d[0], d[1], d[2], nf);
  launched("norm");
  remesh<<<grid(g.N), nthread>>>(g, 2, 2, 2, nf, d[0], d[1], d[2], d[3], d[4], d[5], d[6],
                                 d[7], d[8]);
  launched("remesh");
  scatter<<<grid(g.N), nthread>>>(g, nf, d[0], d[1], d[2], d[3], d[4], d[5], d[9], d[10],
                                  d[11]);
  launched("scatter");
  for (int c = 6; c < 12; c++)
    cuda(cudaMemcpy(h[c], d[c], g.N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
  double gs = 0, gm = 0;
  for (long m = 0; m < g.N; m++)
    for (int c = 0; c < 3; c++) {
      gs = fmax(gs, fabs(h[6 + c][m] - h[9 + c][m]));
      gm = fmax(gm, fabs(h[6 + c][m]));
    }
  printf("gather-scatter maxabs % .16e max % .16e\n", gs, gm);
  double em = 0;
  for (int a = 0; a <= 7; a++)
    for (int b = 0; a + b <= 7; b++)
      for (int c = 0; a + b + c <= 7; c++) {
        double P0 = 0, G0 = 0, S = 0;
        for (long m = 0; m < g.N; m++) {
          long i = m / (n * n), j = m / n % n, k = m % n;
          double x = i + h[0][m] - 16, y = j + h[1][m] - 16, z = k + h[2][m] - 16;
          double f = cpow(x, a) * cpow(y, b) * cpow(z, c);
          P0 += f * h[3][m];
          S += fabs(f * h[3][m]);
          G0 += cpow(i - 16.0, a) * cpow(j - 16.0, b) * cpow(k - 16.0, c) * h[6][m];
        }
        double r = fabs(P0 - G0) / S;
        if (a <= 6 && b <= 6 && c <= 6 && a + b + c <= 6)
          em = fmax(em, r);
        if (a + b + c == 7 && (a == 7 || b == 7 || c == 7))
          printf("moment %d %d %d rel % .16e\n", a, b, c, r);
      }
  printf("remesh moments total degree <= 6 max rel % .16e\n", em);
  for (long m = 0; m < g.N; m++)
    for (int c = 0; c < 3; c++)
      h[c][m] = 1.6 * (hrand(&seed) - 0.5);
  double ie = 0, im = 0;
  for (long m = 0; m < g.N; m++) {
    long i = m / (n * n), j = m / n % n, k = m % n;
    if (i < 8 || i >= 24 || j < 8 || j >= 24 || k < 8 || k >= 24)
      continue;
    double e[3] = {h[0][m], h[1][m], h[2][m]}, v = 0, x = i + e[0] - 16,
           y = j + e[1] - 16, z = k + e[2] - 16;
    long b[3];
    for (int dd = 0; dd < 3; dd++)
      b[dd] = (long)floor(e[dd]);
    for (int p = 0; p < 8; p++)
      for (int q = 0; q < 8; q++)
        for (int r = 0; r < 8; r++) {
          double w = kw(K8, 4, b[0] - 3 + p - e[0]) * kw(K8, 4, b[1] - 3 + q - e[1]) *
                     kw(K8, 4, b[2] - 3 + r - e[2]);
          double X = i + b[0] - 3 + p - 16, Y = j + b[1] - 3 + q - 16,
                 Z = k + b[2] - 3 + r - 16;
          v += w * (cpow(X, 6) - 2 * cpow(X, 3) * cpow(Y, 2) * Z + cpow(Y, 4) * cpow(Z, 2) + 3);
        }
    double f = cpow(x, 6) - 2 * cpow(x, 3) * cpow(y, 2) * z + cpow(y, 4) * cpow(z, 2) + 3;
    ie = fmax(ie, fabs(v - f));
    im = fmax(im, fabs(f));
  }
  printf("interp degree-6 polynomial maxabs % .16e max % .16e\n", ie, im);
  for (long m = 0; m < g.N; m++) {
    long i = m / (n * n), j = m / n % n, k = m % n;
    double x = i * g.hx, y = j * g.hy, z = k * g.hz;
    h[0][m] = -cos(x) * sin(y) * sin(z);
    h[1][m] = -sin(x) * cos(y) * sin(z);
    h[2][m] = 2 * sin(x) * sin(y) * cos(z);
  }
  for (int c = 0; c < 3; c++)
    cuda(cudaMemcpy(d[c], h[c], g.N * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
  plan(&fplan, n, n, n, CUFFT_D2Z);
  plan(&bplan, n, n, n, CUFFT_Z2D);
  for (int c = 0; c < 3; c++)
    cuda(cudaMalloc(&oh[c], g.nf * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&sc, g.nf * sizeof(double2)), "cudaMalloc");
  T = dalloc(g.N);
  double *om[3] = {d[0], d[1], d[2]}, *u[3] = {d[3], d[4], d[5]}, *r[3] = {d[6], d[7], d[8]};
  eval(g, om, u, r, 1, NULL);
  for (int c = 3; c < 9; c++)
    cuda(cudaMemcpy(h[c], d[c], g.N * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
  double eu = 0, er = 0;
  for (long m = 0; m < g.N; m++) {
    long i = m / (n * n), j = m / n % n, k = m % n;
    double x = i * g.hx, y = j * g.hy, z = k * g.hz;
    double U[3] = {sin(x) * cos(y) * cos(z), -cos(x) * sin(y) * cos(z), 0};
    double gu[3][3] = {{cos(x) * cos(y) * cos(z), sin(x) * sin(y) * cos(z), 0},
                       {-sin(x) * sin(y) * cos(z), -cos(x) * cos(y) * cos(z), 0},
                       {-sin(x) * cos(y) * sin(z), cos(x) * sin(y) * sin(z), 0}};
    double w[3] = {h[0][m], h[1][m], h[2][m]};
    for (int c = 0; c < 3; c++) {
      eu = fmax(eu, fabs(h[3 + c][m] - U[c]));
      double R = -g.nu * 3 * w[c];
      for (int j2 = 0; j2 < 3; j2++)
        R += w[j2] * gu[c][j2];
      er = fmax(er, fabs(h[6 + c][m] - R));
    }
  }
  printf("poisson velocity maxabs % .16e\n", eu);
  printf("rhs (transpose stretching + diffusion) maxabs % .16e\n", er);
  printf("gmax % .16e (exact 1)\n", gmax);
  for (long n = 32; n <= 128; n *= 2)
    tsplit(n, 0.2);
  for (double amp = 0.1; amp > 0.006; amp /= 2)
    tsplit(64, amp);
  for (long n = 32; n <= 128; n *= 2) {
    tsplit(n, 0.1);
    tsplit(n, 0.025);
  }
}
int main(int argc, char **argv) {
  G g;
  FILE *file;
  g.f = 0;
  char path[FILENAME_MAX], *ipath, *end;
  double nu, dt, T1, t, x, X, Y, Z, sums[4], mx[4], sp[4], im[7], dd, F0 = 0, C = 0;
  long nx, ny, nz, nt, ns, nd, np, tstep;
  int test, kern = 2, R = 0, M = 0;
  double O = 0, tout = 0;
  (void)argc;
  ipath = NULL;
  nx = ny = nz = 0;
  X = Y = Z = 2;
  nu = dt = -1;
  T1 = dd = 0;
  ns = 10;
  np = 0;
  test = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    if (argv[0][1] == 'h') {
      fprintf(stderr, "Usage: vm3 -x <nx> -y <ny> -z <nz> [-X <Lx/pi>] [-Y <Ly/pi>] "
                      "[-Z <Lz/pi>] -i <u.raw> -n <viscosity> -t <end time> -s <time "
                      "step> [-o <diag every steps>] [-p <reproject every steps>] [-d "
                      "<dump interval>] [-S <stretching: 0 transpose, 1 classic, 2 transpose with curl u, 3 mixed (classic + transpose)/2>] [-f <0 none, 1 2/3 rule, 2 Hou-Li filter>] [-C <Lagrangian CFL: dt = C / max|grad u| each step, replaces -s>] [-O <output time interval with -C>] [-R 1 (Eulerian RK stages, one remesh per step)] [-M 1 (dimension-split remeshing and interpolation)] [-k <remeshing kernel 4 (M4'), 6, 8>] [-v (log magnitude ranges)] [-F <t0>: the filter applied at each projection is mask^(p dt / t0), the same damping per unit time for any dt]\n       vm3 -T (unit tests)\n");
      exit(1);
    }
    if (argv[0][1] == 'T') {
      test = 1;
      continue;
    }
    if (argv[0][1] == 'v') {
      verbose = 1;
      continue;
    }
    if (argv[1] == NULL) {
      fprintf(stderr, "vm3: error: %s needs an argument\n", argv[0]);
      exit(1);
    }
    if (argv[0][1] == 'i') {
      ipath = *++argv;
      continue;
    }
    x = strtod(argv[1], &end);
    if (*end != '\0') {
      fprintf(stderr, "vm3: error: '%s' is not a number\n", argv[1]);
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
      T1 = x;
      break;
    case 'o':
      ns = x;
      break;
    case 'p':
      np = x;
      break;
    case 'd':
      dd = x;
      break;
    case 'S':
      sform = x;
      break;
    case 'f':
      g.f = x;
      break;
    case 'F':
      F0 = x;
      break;
    case 'C':
      C = x;
      break;
    case 'R':
      R = x;
      break;
    case 'M':
      M = x;
      break;
    case 'O':
      O = x;
      break;
    case 'k':
      kern = x == 4 ? 0 : x == 6 ? 1 : x == 8 ? 2 : -1;
      break;
    default:
      fprintf(stderr, "vm3: error: unknown option '%s'\n", *argv);
      exit(1);
    }
    argv++;
  }
  if (kern < 0 || kern > 2) {
    fprintf(stderr, "vm3: error: -k must be 4, 6 or 8\n");
    exit(1);
  }
  cuda(cudaMemcpyToSymbol(cK, test ? KT[2] : KT[kern], sizeof KT[0]), "cudaMemcpyToSymbol");
  cuda(cudaMemcpyToSymbol(cS, test ? &KS[2] : &KS[kern], sizeof(int)), "cudaMemcpyToSymbol");
  cuda(cudaMalloc(&dP, 16 * nblock * sizeof(double)), "cudaMalloc");
  P = (double *)malloc(16 * nblock * sizeof(double));
  if (test) {
    tests();
    return 0;
  }
  if (C > 0) {
    if (dd > 0) {
      fprintf(stderr, "vm3: error: -d needs a fixed time step (no -C)\n");
      exit(1);
    }
    dt = 1;
  }
  if (nx < 8 || ny < 8 || nz < 8 || nz % 2 != 0 || nu < 0 || dt <= 0 || T1 <= 0 ||
      ipath == NULL || ns < 1) {
    fprintf(stderr, "vm3: error: need -x -y -z (>= 8, nz even), -n, -s, -t, -i\n");
    exit(1);
  }
  if (nx > (1L << 20) || ny > (1L << 20) || nz > (1L << 20) ||
      (double)nx * ny * nz * 8 * 32 > 9.0e18) {
    fprintf(stderr, "vm3: error: grid %ld x %ld x %ld too large for 64-bit indexing\n", nx,
            ny, nz);
    exit(1);
  }
  if (T1 / dt > 1e15) {
    fprintf(stderr, "vm3: error: too many steps\n");
    exit(1);
  }
  nt = lround(T1 / dt);
  if (C <= 0 && fabs(nt * dt - T1) > 1e-9 * T1) {
    fprintf(stderr, "vm3: error: -t %g is not a multiple of -s %g\n", T1, dt);
    exit(1);
  }
  nd = dd > 0 ? lround(dd / dt) : 0;
  if (dd > 0 && fabs(nd * dt - dd) > 1e-9 * dd) {
    fprintf(stderr, "vm3: error: -d %g is not a multiple of -s %g\n", dd, dt);
    exit(1);
  }
  g.nx = nx;
  g.ny = ny;
  g.nz = nz;
  g.nzf = nz / 2 + 1;
  g.nf = nx * ny * g.nzf;
  g.N = nx * ny * nz;
  g.ax = 2 / X;
  g.ay = 2 / Y;
  g.az = 2 / Z;
  g.hx = X * pi / nx;
  g.hy = Y * pi / ny;
  g.hz = Z * pi / nz;
  g.s = 1.0 / g.N;
  g.nu = nu;
  g.fp = F0 > 0 ? np * dt / F0 : 1;
  size_t mfree0, mfree1, mtot;
  cuda(cudaMemGetInfo(&mfree0, &mtot), "cudaMemGetInfo");
  plan(&fplan, nx, ny, nz, CUFFT_D2Z);
  plan(&bplan, nx, ny, nz, CUFFT_Z2D);
  for (int c = 0; c < 3; c++)
    cuda(cudaMalloc(&oh[c], g.nf * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&sc, g.nf * sizeof(double2)), "cudaMalloc");
  T = dalloc(g.N);
  if (M != 1 || verbose)
    nf = dalloc(g.N);
  if (sform < 0 || sform > 3) {
    fprintf(stderr, "vm3: error: -S must be 0, 1, 2 or 3\n");
    exit(1);
  }
  if (g.f < 0 || g.f > 2) {
    fprintf(stderr, "vm3: error: -f must be 0, 1 or 2\n");
    exit(1);
  }
  double *kE[3];
  for (int c = 0; c < 3 && R; c++)
    kE[c] = dalloc(g.N);
  if (M) {
    fy = dalloc(g.N);
    fz1 = dalloc(g.N);
    fz = dalloc(g.N);
    T2 = dalloc(g.N);
    itn = dalloc(g.N);
  }
  for (int c = 0; c < 3 && (sform == 2 || g.f); c++)
    wu[c] = dalloc(g.N);
  double *w0[3], *e[3], *a[3], *ae[3], *aw[3], *om[3], *u[3], *r[3], *host;
  for (int c = 0; c < 3; c++) {
    w0[c] = dalloc(g.N);
    e[c] = dalloc(g.N);
    a[c] = dalloc(g.N);
    ae[c] = dalloc(g.N);
    aw[c] = dalloc(g.N);
    om[c] = dalloc(g.N);
    u[c] = dalloc(g.N);
    r[c] = dalloc(g.N);
  }
  cuda(cudaMemGetInfo(&mfree1, &mtot), "cudaMemGetInfo");
  fprintf(stderr, "vm3: device memory %.6f GB (%.1f bytes per node)\n",
          (mfree0 - mfree1) / 1e9, (double)(mfree0 - mfree1) / g.N);
  if ((host = (double *)malloc(g.N * sizeof(double))) == NULL) {
    fprintf(stderr, "vm3: error: malloc failed\n");
    exit(1);
  }
  if ((file = fopen(ipath, "r")) == NULL) {
    fprintf(stderr, "vm3: error: fail to open '%s'\n", ipath);
    exit(1);
  }
  for (int c = 0; c < 3; c++) {
    if (fread(host, sizeof(double), g.N, file) != (size_t)g.N) {
      fprintf(stderr, "vm3: error: fail to read '%s' (wrong -x -y -z?)\n", ipath);
      exit(1);
    }
    cuda(cudaMemcpy(u[c], host, g.N * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    fwd(u[c], oh[c]);
  }
  if (fgetc(file) != EOF) {
    fprintf(stderr, "vm3: error: '%s' is too long (wrong -x -y -z?)\n", ipath);
    exit(1);
  }
  fclose(file);
  for (int c = 0; c < 3; c++) {
    curl<<<grid(g.nf), nthread>>>(g, c, oh[0], oh[1], oh[2], sc);
    launched("curl");
    bwd(sc, w0[c]);
  }
  printf("step t E Omega Omegau wmax gmax divr sumx sumy sumz H Epsi Ix Iy Iz Ax Ay Az dt wall\n");
  double t0 = now();
  const double bb[4] = {1. / 6, 1. / 3, 1. / 3, 1. / 6}, cc[4] = {0.5, 0.5, 1, 0};
  t = 0;
  for (tstep = 0;; tstep++) {
    int last;
    if (C > 0)
      last = t >= T1 * (1 - 1e-14);
    else {
      t = tstep * dt;
      last = tstep >= nt;
    }
    int diag = (C > 0 && O > 0) ? (last || t >= tout * (1 - 1e-14)) : (tstep % ns == 0 || last);
    eval(g, w0, u, r, diag || C > 0, R ? kE : NULL);
    if (diag && C > 0 && O > 0)
      tout += O;
    if (C > 0 && !last) {
      dt = C / gmax;
      if (O > 0 && t + dt > tout)
        dt = tout - t;
      if (t + dt > T1)
        dt = T1 - t;
      g.fp = F0 > 0 ? np * dt / F0 : 1;
    }
    if (diag) {
      cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
      rsum<<<nblock, nthread>>>(g.N, u[0], u[1], u[2], dP);
      launched("rsum");
      collect2(4, sums);
      double E = 0.5 * sums[0] * g.s;
      rsum<<<nblock, nthread>>>(g.N, w0[0], w0[1], w0[2], dP);
      launched("rsum");
      collect2(4, sums);
      double Om = 0.5 * sums[0] * g.s;
      rmax<<<nblock, nthread>>>(g.N, w0[0], w0[1], w0[2], dP);
      launched("rmax");
      collect(4, 1, mx);
      for (int c = 0; c < 3; c++)
        fwd(w0[c], oh[c]);
      rspec<<<nblock, nthread>>>(g, oh[0], oh[1], oh[2], dP);
      launched("rspec");
      collect2(4, sp);
      rimp<<<nblock, nthread>>>(g, w0[0], w0[1], w0[2], u[0], u[1], u[2], dP);
      launched("rimp");
      collect2(7, im);
      double q[] = {t,
                    E,
                    Om,
                    sp[0] * g.s * g.s,
                    mx[0],
                    gmax,
                    sp[2] > 0 ? sp[1] / sp[2] : 0,
                    sums[1] * g.s,
                    sums[2] * g.s,
                    sums[3] * g.s,
                    im[6] * g.s,
                    sp[3] * g.s * g.s,
                    im[0] * g.s,
                    im[1] * g.s,
                    im[2] * g.s,
                    im[3] * g.s,
                    im[4] * g.s,
                    im[5] * g.s,
                    dt,
                    now() - t0};
      if (verbose) {
        const char *nm[3] = {"x", "y", "z"};
        char b[32];
        for (int c = 0; c < 3; c++) {
          sprintf(b, "w%s", nm[c]);
          look(g.N, b, w0[c]);
          sprintf(b, "u%s", nm[c]);
          look(g.N, b, u[c]);
          sprintf(b, "r%s", nm[c]);
          look(g.N, b, r[c]);
          if (tstep > 0) {
            sprintf(b, "e%s", nm[c]);
            look(g.N, b, e[c]);
          }
        }
        if (tstep > 0) {
          norm<<<grid(g.N), nthread>>>(g, e[0], e[1], e[2], nf);
          launched("norm");
          addc<<<grid(g.N), nthread>>>(g.N, nf, -1, T);
          launched("addc");
          look(g.N, "1/sumW-1", T);
        }
      }
      printf("% 10ld", tstep);
      for (int i = 0; i < (int)(sizeof q / sizeof *q); i++)
        printf(" % .16e", q[i]);
      printf("\n");
      if (fflush(stdout) != 0 || ferror(stdout)) {
        fprintf(stderr, "vm3: error: fail to write stdout\n");
        exit(1);
      }
    }
    if (nd > 0 && tstep % nd == 0) {
      sprintf(path, "u.%08ld", tstep);
      if ((file = fopen(path, "w")) == NULL) {
        fprintf(stderr, "vm3: error: fail to open '%s'\n", path);
        exit(1);
      }
      for (int c = 0; c < 3; c++) {
        cuda(cudaMemcpy(host, u[c], g.N * sizeof(double), cudaMemcpyDeviceToHost),
             "cudaMemcpy");
        if (fwrite(host, sizeof(double), g.N, file) != (size_t)g.N) {
          fprintf(stderr, "vm3: error: fail to write '%s'\n", path);
          exit(1);
        }
      }
      if (ferror(file) || fclose(file) != 0) {
        fprintf(stderr, "vm3: error: fail to write '%s'\n", path);
        exit(1);
      }
    }
    if (last)
      break;
    for (int s = 0; s < 4; s++) {
      if (s > 0) {
        if (R) {
          for (int c = 0; c < 3; c++) {
            lin<<<grid(g.N), nthread>>>(g.N, w0[c], cc[s - 1] * dt, kE[c], om[c]);
            launched("lin");
          }
          eval(g, om, u, r, 0, s < 3 ? kE : NULL);
          if (M)
            factor(g, e);
        } else if (M) {
          factor(g, e);
          splitremesh(g, e, a, om);
          eval(g, om, u, r, 0, NULL);
        } else {
          doremesh(g, e, a, om);
          eval(g, om, u, r, 0, NULL);
        }
        if (M) {
          double *f[6] = {u[0], u[1], u[2], r[0], r[1], r[2]};
          splitinterp(g, e, f, 6);
        }
      }
      stage<<<grid(g.N), nthread>>>(g, s == 0, M && s > 0, bb[s], cc[s], dt, u[0], u[1],
                                    u[2], r[0], r[1], r[2], e[0], e[1], e[2], a[0], a[1],
                                    a[2], w0[0], w0[1], w0[2], ae[0], ae[1], ae[2], aw[0],
                                    aw[1], aw[2]);
      launched("stage");
    }
    if (M) {
      factor(g, e);
      splitremesh(g, e, a, w0);
    } else
      doremesh(g, e, a, w0);
    if (np > 0 && (tstep + 1) % np == 0)
      project(g, w0);
    if (C > 0)
      t += dt;
  }
  return 0;
}
