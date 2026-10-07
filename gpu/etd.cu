#include <cufft.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
enum { SIN, COS };
enum { EVEN, ODD };
enum { R00C, R00S, R10C, R10S, R01C, R01S };
enum { nthread = 256, nblock = 1024 };
struct V {
  double *a[2][3];
};
static long N, N1, n3, B;
static double *R, *T, *cs, *sn, *F[3], *tmp;
static double2 *C;
static cufftHandle plans[2][2];
static V w, U, W;
static const int tu[3][3] = {{SIN, COS, COS}, {COS, SIN, COS}, {COS, COS, SIN}};
static const int tw[3][3] = {{COS, SIN, SIN}, {SIN, COS, SIN}, {SIN, SIN, COS}};
static void cuda(cudaError_t e, const char *what) {
  if (e != cudaSuccess) {
    fprintf(stderr, "etd: error: %s: %s\n", what, cudaGetErrorString(e));
    exit(1);
  }
}
static void cufft(cufftResult r, const char *what) {
  if (r != CUFFT_SUCCESS) {
    fprintf(stderr, "etd: error: %s: cufft error %d\n", what, (int)r);
    exit(1);
  }
}
static void launched(const char *what) { cuda(cudaGetLastError(), what); }
static unsigned grid(long n) { return (n + nthread - 1) / nthread; }
__host__ __device__ static long wavenumber(int c, long n) {
  return c == EVEN ? 2 * n : 2 * n + 1;
}
__device__ static void line(int ax, long t, long L, long N1, long B, long *m,
                            long *r, long *base, long *s) {
  long b;
  if (ax < 2) {
    b = t % B;
    *m = t / B;
    *r = *m * B + b;
  } else {
    *m = t % L;
    b = t / L;
    *r = b * L + *m;
  }
  *base = ax == 0 ? b : ax == 1 ? b / N1 * N1 * N1 + b % N1 : b * N1;
  *s = ax == 0 ? N1 * N1 : ax == 1 ? N1 : 1;
}
__global__ static void trig(long N, double *cs, double *sn) {
  long k = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (k > N)
    return;
  sincospi((double)k / (2 * N), &sn[k], &cs[k]);
}
__global__ static void extend(int ax, int kd, long N, long N1, long B,
                              const double *a, double *R) {
  long t = blockIdx.x * (long)blockDim.x + threadIdx.x, n = 2 * N, m, r, o, s;
  double x;
  if (t >= n * B)
    return;
  line(ax, t, n, N1, B, &m, &r, &o, &s);
  if (kd == R00C)
    x = a[o + s * (m <= N ? m : n - m)];
  else if (kd == R00S)
    x = m == 0 || m == N ? 0 : m < N ? a[o + s * m] : -a[o + s * (n - m)];
  else if (kd == R10C)
    x = a[o + s * (m < N ? m : n - 1 - m)];
  else
    x = m < N ? a[o + s * m] : -a[o + s * (n - 1 - m)];
  R[r] = x;
}
__global__ static void fold(int ax, int kd, long N, long N1, long B,
                            const double *cs, const double *sn,
                            const double2 *C, double *b) {
  long t = blockIdx.x * (long)blockDim.x + threadIdx.x, k, r, o, s;
  double y;
  double2 c;
  if (t >= N1 * B)
    return;
  line(ax, t, N1, N1, B, &k, &r, &o, &s);
  c = C[r];
  if (kd == R00C)
    y = c.x;
  else if (kd == R00S)
    y = k == 0 || k == N ? 0 : -c.y;
  else if (kd == R10C)
    y = k < N ? c.x * cs[k] + c.y * sn[k] : 0;
  else
    y = k > 0 ? c.x * sn[k] - c.y * cs[k] : 0;
  b[o + s * k] = y;
}
__global__ static void twiddle(int ax, int kd, long N, long N1, long B,
                               const double *cs, const double *sn,
                               const double *a, double2 *C) {
  long t = blockIdx.x * (long)blockDim.x + threadIdx.x, j, r, o, s;
  double x;
  if (t >= N1 * B)
    return;
  line(ax, t, N1, N1, B, &j, &r, &o, &s);
  if (kd == R01C) {
    x = j < N ? a[o + s * j] : 0;
    C[r] = make_double2(x * cs[j], x * sn[j]);
  } else {
    x = j > 0 ? a[o + s * j] : 0;
    C[r] = j == N ? make_double2(x, 0) : make_double2(x * sn[j], -x * cs[j]);
  }
}
__global__ static void crop(int ax, long N, long N1, long B, const double *R,
                            double *b) {
  long t = blockIdx.x * (long)blockDim.x + threadIdx.x, k, r, o, s;
  if (t >= N1 * B)
    return;
  line(ax, t, N1, N1, B, &k, &r, &o, &s);
  b[o + s * k] = k < N ? R[ax < 2 ? r : t / N1 * 2 * N + k] : 0;
}
__device__ static double weight(int syn, int c, int type, long i, long N) {
  double s = c == EVEN && type == COS && (i == 0 || i == N) ? 1 : 0.5;
  return syn ? s : 1 / (2 * N * s);
}
__global__ static void scale(int syn, int c, int t0, int t1, int t2, long N,
                             long N1, long n3, const double *a, double *b) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  if (l >= n3)
    return;
  i = l / (N1 * N1);
  j = l / N1 % N1;
  k = l % N1;
  b[l] = a[l] * (weight(syn, c, t0, i, N) * weight(syn, c, t1, j, N) *
                 weight(syn, c, t2, k, N));
}
__global__ static void curl(int c, long N1, long n3, V u, V w) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n3)
    return;
  double m = wavenumber(c, l / (N1 * N1)), n = wavenumber(c, l / N1 % N1),
         p = wavenumber(c, l % N1);
  w.a[c][0][l] = -n * u.a[c][2][l] + p * u.a[c][1][l];
  w.a[c][1][l] = -p * u.a[c][0][l] + m * u.a[c][2][l];
  w.a[c][2][l] = -m * u.a[c][1][l] + n * u.a[c][0][l];
}
__device__ static void blocksum(int K, const double *v, double *s,
                                double *out) {
  int i = threadIdx.x;
  for (int q = 0; q < K; q++)
    s[q * nthread + i] = v[q];
  __syncthreads();
  for (int h = nthread / 2; h > 0; h /= 2) {
    if (i < h)
      for (int q = 0; q < K; q++)
        s[q * nthread + i] += s[q * nthread + i + h];
    __syncthreads();
  }
  if (i == 0)
    for (int q = 0; q < K; q++)
      out[blockIdx.x * K + q] = s[q * nthread];
}
__global__ static void norms(long N1, long n3, long kcut, V u, V w,
                             double *out) {
  __shared__ double s[5 * nthread];
  double v[5] = {0, 0, 0, 0, 0};
  for (long t = blockIdx.x * (long)blockDim.x + threadIdx.x; t < 2 * n3;
       t += (long)gridDim.x * blockDim.x) {
    int c = t / n3;
    long l = t % n3, g[3] = {wavenumber(c, l / (N1 * N1)),
                             wavenumber(c, l / N1 % N1), wavenumber(c, l % N1)};
    for (int d = 0; d < 3; d++) {
      double su = 1, sw = 1, x = u.a[c][d][l], y = w.a[c][d][l];
      for (int e = 0; e < 3; e++) {
        su *= e != d && g[e] == 0 ? 1 : 0.5;
        sw *= e == d && g[e] == 0 ? 1 : 0.5;
      }
      v[0] += su * x * x;
      v[1] += sw * y * y;
      v[2] += sw * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) * y * y;
      if (l / (N1 * N1) >= kcut - 1 || l / N1 % N1 >= kcut - 1 ||
          l % N1 >= kcut - 1)
        v[3] += su * x * x;
      if (g[0] != 0 || g[1] != 0 || g[2] != 0)
        v[4] += su * x * x / sqrt((double)(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]));
    }
  }
  blocksum(5, v, s, out);
}
__global__ static void peaks(long n3, V U, V W, double *out) {
  __shared__ double s[2 * nthread];
  double v[2] = {0, 0};
  int i = threadIdx.x;
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < n3;
       l += (long)gridDim.x * blockDim.x)
    for (int q = -1; q <= 1; q += 2) {
      double a = 0, b = 0, x;
      for (int d = 0; d < 3; d++) {
        a += fabs(U.a[EVEN][d][l] + q * U.a[ODD][d][l]);
        x = W.a[EVEN][d][l] + q * W.a[ODD][d][l];
        b += x * x;
      }
      v[0] = fmax(v[0], a);
      v[1] = fmax(v[1], b);
    }
  for (int q = 0; q < 2; q++)
    s[q * nthread + i] = v[q];
  __syncthreads();
  for (int h = nthread / 2; h > 0; h /= 2) {
    if (i < h)
      for (int q = 0; q < 2; q++)
        s[q * nthread + i] = fmax(s[q * nthread + i], s[q * nthread + i + h]);
    __syncthreads();
  }
  if (i == 0)
    for (int q = 0; q < 2; q++)
      out[blockIdx.x * 2 + q] = s[q * nthread];
}
__global__ static void deriv(int c, int q, long N1, long n3, const double *u,
                             double *du) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n3)
    return;
  double m = wavenumber(c, l / (N1 * N1));
  du[l] = q == 0 ? m * u[l] : -m * m * u[l];
}
__global__ static void moments(long N, long N1, long n3, const double *fe,
                               const double *fo, double *out) {
  __shared__ double s[9 * nthread];
  double m[9] = {0}, w, x, y, z;
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < n3;
       l += (long)gridDim.x * blockDim.x) {
    long i = l / (N1 * N1), j = l / N1 % N1, k = l % N1;
    w = (i == 0 || i == N ? 0.5 : 1) * (j == 0 || j == N ? 0.5 : 1) *
        (k == 0 || k == N ? 0.5 : 1);
    x = fe[l] + fo[l];
    z = fe[l] - fo[l];
    y = w / 2;
    for (int n = 0; n < 9; n++) {
      m[n] += y;
      y *= x;
    }
    y = w / 2;
    for (int n = 0; n < 9; n++) {
      m[n] += y;
      y *= z;
    }
  }
  blocksum(9, m, s, out);
}
__global__ static void spectrum(long N1, long n3, long nb, V u, V f,
                                double *E) {
  long t = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (t >= 2 * n3)
    return;
  int c = t / n3;
  long l = t % n3, g[3] = {wavenumber(c, l / (N1 * N1)),
                           wavenumber(c, l / N1 % N1), wavenumber(c, l % N1)};
  long b = (long)(2 * sqrt((double)(g[0] * g[0] + g[1] * g[1] + g[2] * g[2])));
  double e = 0, tr = 0;
  for (int d = 0; d < 3; d++) {
    double su = 1;
    for (int q = 0; q < 3; q++)
      su *= q != d && g[q] == 0 ? 1 : 0.5;
    e += su * u.a[c][d][l] * u.a[c][d][l] / 2;
    tr += su * u.a[c][d][l] * f.a[c][d][l];
  }
  if (b < nb && e != 0) {
    atomicAdd(&E[b], e);
    atomicAdd(&E[nb + b], tr);
  }
}
__global__ static void cross(int c, long n3, V U, V W, double *F0, double *F1,
                             double *F2) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n3)
    return;
  double *const *ue = U.a[EVEN], *const *uo = U.a[ODD], *const *wa = W.a[c],
                *const *wb = W.a[1 - c];
  F0[l] = ue[1][l] * wa[2][l] - ue[2][l] * wa[1][l] + uo[1][l] * wb[2][l] -
          uo[2][l] * wb[1][l];
  F1[l] = ue[2][l] * wa[0][l] - ue[0][l] * wa[2][l] + uo[2][l] * wb[0][l] -
          uo[0][l] * wb[2][l];
  F2[l] = ue[0][l] * wa[1][l] - ue[1][l] * wa[0][l] + uo[0][l] * wb[1][l] -
          uo[1][l] * wb[0][l];
}
__device__ static double2 cm(double2 a, double2 b) {
  return make_double2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
__device__ static double2 cd(double2 a, double2 b) {
  double d = b.x * b.x + b.y * b.y;
  return make_double2((a.x * b.x + a.y * b.y) / d, (a.y * b.x - a.x * b.y) / d);
}
__device__ static double2 cx(double2 a) {
  double e = exp(a.x);
  return make_double2(e * cos(a.y), e * sin(a.y));
}
__device__ static double2 poly(double a0, double a1, double a2, double2 L,
                               double2 L2) {
  return make_double2(a0 + a1 * L.x + a2 * L2.x, a1 * L.y + a2 * L2.y);
}
__device__ static double2 add(double2 a, double2 b) {
  return make_double2(a.x + b.x, a.y + b.y);
}
__global__ static void coef(long nk, double nu, double dt, double *tab) {
  long q = blockIdx.x * (long)blockDim.x + threadIdx.x;
  double z, m[4] = {0, 0, 0, 0}, sn, cs;
  double2 L, L2, L3, e, e2;
  if (q >= nk)
    return;
  z = -nu * dt * q;
  for (int j = 0; j < 16; j++) {
    sincospi((j + 0.5) / 16, &sn, &cs);
    L = make_double2(z + cs, sn);
    L2 = cm(L, L);
    L3 = cm(L2, L);
    e = cx(L);
    e2 = cx(make_double2(L.x / 2, L.y / 2));
    m[0] += cd(make_double2(e2.x - 1, e2.y), L).x;
    m[1] += cd(add(poly(-4, -1, 0, L, L2), cm(e, poly(4, -3, 1, L, L2))), L3).x;
    m[2] += cd(add(poly(2, 1, 0, L, L2), cm(e, poly(-2, 1, 0, L, L2))), L3).x;
    m[3] += cd(add(poly(-4, -3, -1, L, L2), cm(e, poly(4, -1, 0, L, L2))), L3).x;
  }
  tab[q] = exp(z);
  tab[nk + q] = exp(z / 2);
  for (int i = 0; i < 4; i++)
    tab[(2 + i) * nk + q] = dt * m[i] / 16;
}
__global__ static void project(int c, long N1, long n3, long kcut, V f) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n3)
    return;
  long i = l / (N1 * N1), j = l / N1 % N1, k = l % N1;
  double m = wavenumber(c, i), n = wavenumber(c, j), p = wavenumber(c, k),
         kk = m * m + n * n + p * p, P, g[3];
  for (int d = 0; d < 3; d++)
    g[d] = f.a[c][d][l];
  if (i > kcut || j > kcut || k > kcut)
    g[0] = g[1] = g[2] = 0;
  P = kk > 0 ? (m * g[0] + n * g[1] + p * g[2]) / kk : 0;
  f.a[c][0][l] = g[0] - m * P;
  f.a[c][1][l] = g[1] - n * P;
  f.a[c][2][l] = g[2] - p * P;
}
__global__ static void stage(int q, int c, long N1, long n3, long nk,
                             const double *tab, V u, V X, V Nv, V acc, V s) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n3)
    return;
  long m = wavenumber(c, l / (N1 * N1)), n = wavenumber(c, l / N1 % N1),
       p = wavenumber(c, l % N1), kk = m * m + n * n + p * p;
  double E = tab[kk], E2 = tab[nk + kk], Q = tab[2 * nk + kk],
         f1 = tab[3 * nk + kk], f2 = tab[4 * nk + kk], f3 = tab[5 * nk + kk];
  for (int d = 0; d < 3; d++) {
    double v = u.a[c][d][l], x = X.a[c][d][l], y = Nv.a[c][d][l];
    if (q == 0) {
      acc.a[c][d][l] = E * v + f1 * y;
      s.a[c][d][l] = E2 * v + Q * y;
    } else if (q == 1) {
      acc.a[c][d][l] += 2 * f2 * x;
      s.a[c][d][l] = E2 * v + Q * x;
    } else if (q == 2) {
      acc.a[c][d][l] += 2 * f2 * x;
      s.a[c][d][l] = E2 * (E2 * v + Q * y) + Q * (2 * x - y);
    } else
      u.a[c][d][l] = acc.a[c][d][l] + f3 * x;
  }
}
static int kind(int dir, int c, int type) {
  if (c == EVEN)
    return type == COS ? R00C : R00S;
  if (type == COS)
    return dir == 0 ? R10C : R01C;
  return dir == 0 ? R10S : R01S;
}
static void pass(int ax, int kd, const double *a, double *b) {
  if (kd < R01C) {
    extend<<<grid(2 * N * B), nthread>>>(ax, kd, N, N1, B, a, R);
    launched("extend");
    cufft(cufftExecD2Z(plans[ax == 2][0], R, C), "cufftExecD2Z");
    fold<<<grid(N1 * B), nthread>>>(ax, kd, N, N1, B, cs, sn, C, b);
    launched("fold");
  } else {
    twiddle<<<grid(N1 * B), nthread>>>(ax, kd, N, N1, B, cs, sn, a, C);
    launched("twiddle");
    cufft(cufftExecZ2D(plans[ax == 2][1], C, R), "cufftExecZ2D");
    crop<<<grid(N1 * B), nthread>>>(ax, N, N1, B, R, b);
    launched("crop");
  }
}
static void transform(int dir, int c, const int *t, const double *a,
                      double *b) {
  pass(0, kind(dir, c, t[0]), a, b);
  pass(1, kind(dir, c, t[1]), b, T);
  pass(2, kind(dir, c, t[2]), T, b);
}
static void synthesis(int c, const int *t, const double *u, double *tmp,
                      double *f) {
  scale<<<grid(n3), nthread>>>(1, c, t[0], t[1], t[2], N, N1, n3, u, tmp);
  launched("scale");
  transform(0, c, t, tmp, f);
}
static void analysis(int c, const int *t, double *f, double *u) {
  transform(1, c, t, f, u);
  scale<<<grid(n3), nthread>>>(0, c, t[0], t[1], t[2], N, N1, n3, u, u);
  launched("scale");
}
static void physical(V v) {
  int c, d;
  for (c = 0; c < 2; c++) {
    curl<<<grid(n3), nthread>>>(c, N1, n3, v, w);
    launched("curl");
  }
  for (c = 0; c < 2; c++)
    for (d = 0; d < 3; d++) {
      synthesis(c, tu[d], v.a[c][d], tmp, U.a[c][d]);
      synthesis(c, tw[d], w.a[c][d], tmp, W.a[c][d]);
    }
}
static void rhs(long kcut, V v, V f) {
  int c, d;
  physical(v);
  for (c = 0; c < 2; c++) {
    cross<<<grid(n3), nthread>>>(c, n3, U, W, F[0], F[1], F[2]);
    launched("cross");
    for (d = 0; d < 3; d++)
      analysis(c, tu[d], F[d], f.a[c][d]);
  }
  for (c = 0; c < 2; c++) {
    project<<<grid(n3), nthread>>>(c, N1, n3, kcut, f);
    launched("project");
  }
}
static void plan(cufftHandle *p, cufftType type, long long n, long long st,
                 long long dist, long long ost, long long odist, size_t *ws) {
  long long ne = type == CUFFT_D2Z ? n : n / 2 + 1,
            no = type == CUFFT_D2Z ? n / 2 + 1 : n;
  size_t size;
  cufft(cufftCreate(p), "cufftCreate");
  cufft(cufftSetAutoAllocation(*p, 0), "cufftSetAutoAllocation");
  cufft(cufftMakePlanMany64(*p, 1, &n, &ne, st, dist, &no, ost, odist, type, B,
                            &size),
        "cufftMakePlanMany64");
  if (size > *ws)
    *ws = size;
}
static void collect(int K, const double *dP, double *P, double *sum) {
  cuda(cudaMemcpy(P, dP, nblock * K * sizeof(double), cudaMemcpyDeviceToHost),
       "cudaMemcpy");
  for (int q = 0; q < K; q++)
    sum[q] = 0;
  for (int b = 0; b < nblock; b++)
    for (int q = 0; q < K; q++)
      sum[q] += P[b * K + q];
}
int main(int argc, char **argv) {
  long M, tstep, ne, nd, kcut, nk;
  int c, d;
  double nu, dt, T_end, t, x, e, dd, *dP, *P, *dE, *host, *tab, one;
  size_t ws;
  void *work;
  V u, du, Nv, acc, s;
  char *end;
  (void)argc;
  M = 0;
  nu = -1;
  dt = -1;
  T_end = 0;
  e = 0;
  dd = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    if (argv[0][1] == 'h') {
      fprintf(stderr, "Usage: etd -M <modes> -n <viscosity> -t <end time> -s "
                      "<time step> [-e <spectrum interval>] [-d <dump interval>]\n");
      exit(1);
    }
    if (argv[1] == NULL) {
      fprintf(stderr, "etd: error: %s needs an argument\n", argv[0]);
      exit(1);
    }
    x = strtod(argv[1], &end);
    if (*end != '\0') {
      fprintf(stderr, "etd: error: '%s' is not a number\n", argv[1]);
      exit(1);
    }
    switch (argv[0][1]) {
    case 'M':
      M = x;
      break;
    case 'n':
      nu = x;
      break;
    case 's':
      dt = x;
      break;
    case 't':
      T_end = x;
      break;
    case 'e':
      e = x;
      break;
    case 'd':
      dd = x;
      break;
    default:
      fprintf(stderr, "etd: error: unknown option '%s'\n", *argv);
      exit(1);
    }
    argv++;
  }
  if (M < 4 || M % 2 != 0 || nu < 0 || dt <= 0 || T_end <= 0) {
    fprintf(stderr, "etd: error: need -M (even, >= 4), -n, -s, -t\n");
    exit(1);
  }
  N = M / 2;
  N1 = N + 1;
  n3 = N1 * N1 * N1;
  B = N1 * N1;
  kcut = M / 3;
  ne = e > 0 ? lround(e / dt) : 0;
  nd = dd > 0 ? lround(dd / dt) : 0;
  for (c = 0; c < 2; c++)
    for (d = 0; d < 3; d++) {
      cuda(cudaMalloc(&u.a[c][d], n3 * sizeof(double)), "cudaMalloc");
      cuda(cudaMalloc(&Nv.a[c][d], n3 * sizeof(double)), "cudaMalloc");
      cuda(cudaMalloc(&acc.a[c][d], n3 * sizeof(double)), "cudaMalloc");
      cuda(cudaMalloc(&s.a[c][d], n3 * sizeof(double)), "cudaMalloc");
      cuda(cudaMalloc(&w.a[c][d], n3 * sizeof(double)), "cudaMalloc");
      cuda(cudaMalloc(&du.a[c][d], n3 * sizeof(double)), "cudaMalloc");
      cuda(cudaMalloc(&U.a[c][d], n3 * sizeof(double)), "cudaMalloc");
      cuda(cudaMalloc(&W.a[c][d], n3 * sizeof(double)), "cudaMalloc");
      cuda(cudaMemset(u.a[c][d], 0, n3 * sizeof(double)), "cudaMemset");
    }
  for (d = 0; d < 3; d++)
    cuda(cudaMalloc(&F[d], n3 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&tmp, n3 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&T, n3 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&R, 2 * N * B * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&C, N1 * B * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&cs, N1 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&sn, N1 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&dP, 9 * nblock * sizeof(double)), "cudaMalloc");
  if ((P = (double *)malloc(9 * nblock * sizeof(double))) == NULL) {
    fprintf(stderr, "etd: error: malloc failed\n");
    exit(1);
  }
  nk = 3 * (2 * N + 1) * (2 * N + 1) + 1;
  cuda(cudaMalloc(&tab, 6 * nk * sizeof(double)), "cudaMalloc");
  coef<<<grid(nk), nthread>>>(nk, nu, dt, tab);
  launched("coef");
  trig<<<grid(N1), nthread>>>(N, cs, sn);
  launched("trig");
  ws = 0;
  plan(&plans[0][0], CUFFT_D2Z, 2 * N, B, 1, B, 1, &ws);
  plan(&plans[0][1], CUFFT_Z2D, 2 * N, B, 1, B, 1, &ws);
  plan(&plans[1][0], CUFFT_D2Z, 2 * N, 1, 2 * N, 1, N1, &ws);
  plan(&plans[1][1], CUFFT_Z2D, 2 * N, 1, N1, 1, 2 * N, &ws);
  cuda(cudaMalloc(&work, ws > 0 ? ws : 1), "cudaMalloc");
  for (int a = 0; a < 2; a++)
    for (int b = 0; b < 2; b++)
      cufft(cufftSetWorkArea(plans[a][b], work), "cufftSetWorkArea");
  one = 1;
  cuda(cudaMemcpy(u.a[ODD][0], &one, sizeof(double), cudaMemcpyHostToDevice),
       "cudaMemcpy");
  one = -1;
  cuda(cudaMemcpy(u.a[ODD][1], &one, sizeof(double), cudaMemcpyHostToDevice),
       "cudaMemcpy");
  host = NULL;
  t = 0;
  tstep = 0;
  for (;;) {
    for (c = 0; c < 2; c++) {
      curl<<<grid(n3), nthread>>>(c, N1, n3, u, w);
      launched("curl");
    }
    if (tstep % 10 == 0) {
      double sums[9], S[9], Sb[9];
      double umax = 0, wmax = 0, Eface;
      norms<<<nblock, nthread>>>(N1, n3, kcut, u, w, dP);
      launched("norms");
      collect(5, dP, P, sums);
      double Ek = sums[0] / 2, eps = nu * sums[1], u2 = 2 * Ek / 3;
      double Lint = M_PI / (2 * u2) * sums[4] / 2, lam = sqrt(10 * nu * Ek / eps);
      Eface = sums[3] / sums[0];
      physical(u);
      peaks<<<nblock, nthread>>>(n3, U, W, dP);
      launched("peaks");
      cuda(cudaMemcpy(P, dP, 2 * nblock * sizeof(double),
                      cudaMemcpyDeviceToHost),
           "cudaMemcpy");
      for (int b = 0; b < nblock; b++) {
        umax = fmax(umax, P[2 * b]);
        wmax = fmax(wmax, P[2 * b + 1]);
      }
      double energy = sums[0], Omega = sums[1], Pal = sums[2];
      for (int q = 0; q < 2; q++) {
        for (c = 0; c < 2; c++) {
          static const int tc[3] = {COS, COS, COS}, ts[3] = {SIN, COS, COS};
          deriv<<<grid(n3), nthread>>>(c, q, N1, n3, u.a[c][0], du.a[c][0]);
          launched("deriv");
          synthesis(c, q == 0 ? tc : ts, du.a[c][0], tmp, F[c]);
        }
        moments<<<nblock, nthread>>>(N, N1, n3, F[EVEN], F[ODD], dP);
        launched("moments");
        collect(9, dP, P, sums);
        for (int n = 3; n < 9; n++)
          (q == 0 ? S : Sb)[n] = (n % 2 ? -1 : 1) * (sums[n] / sums[0]) /
                                 pow(sums[2] / sums[0], n / 2.0);
      }
      printf("% 10ld % .16e % .16e % .16e", tstep, t, energy / 2, Omega / 2);
      for (int n = 3; n < 9; n++)
        printf(" % .6e", S[n]);
      printf(" % .6e % .6e % .16e", Sb[4], Sb[6], Pal / 2);
      printf(" % .6e % .6e % .6e % .6e", dt * umax * (2 * kcut + 1),
             sqrt(wmax), (2 * kcut + 1) * pow(nu * nu * nu / (nu * Omega), 0.25),
             Eface);
      printf(" % .6e % .6e % .6e", u2 * sqrt(15 / (nu * eps)), lam, Lint);
      printf("\n");
      fflush(stdout);
    }
    if (ne > 0 && tstep % ne == 0) {
      long nb = 2 * (long)(sqrt(3.0) * (M + 1)) + 2;
      double *E, Pi;
      char path[FILENAME_MAX];
      FILE *file;
      if ((E = (double *)malloc(2 * nb * sizeof(double))) == NULL) {
        fprintf(stderr, "etd: error: malloc failed\n");
        exit(1);
      }
      cuda(cudaMalloc(&dE, 2 * nb * sizeof(double)), "cudaMalloc");
      cuda(cudaMemset(dE, 0, 2 * nb * sizeof(double)), "cudaMemset");
      rhs(kcut, u, du);
      spectrum<<<grid(2 * n3), nthread>>>(N1, n3, nb, u, du, dE);
      launched("spectrum");
      cuda(cudaMemcpy(E, dE, 2 * nb * sizeof(double), cudaMemcpyDeviceToHost),
           "cudaMemcpy");
      cuda(cudaFree(dE), "cudaFree");
      sprintf(path, "e.%08ld", tstep);
      if ((file = fopen(path, "w")) == NULL) {
        fprintf(stderr, "etd: error: fail to open '%s'\n", path);
        exit(1);
      }
      fprintf(file, "# t = %.16e\n", t);
      Pi = 0;
      for (long b = 0; b < nb; b++) {
        Pi -= E[nb + b];
        fprintf(file, "%.1f %.16e % .16e % .16e\n", b / 2.0, E[b], E[nb + b],
                Pi);
      }
      if (fclose(file) != 0) {
        fprintf(stderr, "etd: error: fail to close '%s'\n", path);
        exit(1);
      }
      free(E);
    }
    if (nd > 0 && tstep % nd == 0) {
      char path[FILENAME_MAX];
      FILE *file;
      if (host == NULL &&
          (host = (double *)malloc(n3 * sizeof(double))) == NULL) {
        fprintf(stderr, "etd: error: malloc failed\n");
        exit(1);
      }
      sprintf(path, "u.%08ld", tstep);
      if ((file = fopen(path, "w")) == NULL) {
        fprintf(stderr, "etd: error: fail to open '%s'\n", path);
        exit(1);
      }
      for (c = 0; c < 2; c++)
        for (d = 0; d < 3; d++) {
          cuda(cudaMemcpy(host, u.a[c][d], n3 * sizeof(double),
                          cudaMemcpyDeviceToHost),
               "cudaMemcpy");
          if (fwrite(host, sizeof(double), n3, file) != (size_t)n3) {
            fprintf(stderr, "etd: error: fail to write '%s'\n", path);
            exit(1);
          }
        }
      if (fclose(file) != 0) {
        fprintf(stderr, "etd: error: fail to close '%s'\n", path);
        exit(1);
      }
    }
    if (t > T_end)
      break;
    for (int q = 0; q < 4; q++) {
      rhs(kcut, q == 0 ? u : s, q == 0 ? Nv : du);
      for (c = 0; c < 2; c++) {
        stage<<<grid(n3), nthread>>>(q, c, N1, n3, nk, tab, u, q == 0 ? Nv : du,
                                     Nv, acc, s);
        launched("stage");
      }
    }
    t += dt;
    tstep++;
  }
  cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
}
