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
static double *R, *T, *cs, *sn;
static double2 *C;
static cufftHandle plans[2][2];
static void cuda(cudaError_t e, const char *what) {
  if (e != cudaSuccess) {
    fprintf(stderr, "tg: error: %s: %s\n", what, cudaGetErrorString(e));
    exit(1);
  }
}
static void cufft(cufftResult r, const char *what) {
  if (r != CUFFT_SUCCESS) {
    fprintf(stderr, "tg: error: %s: cufft error %d\n", what, (int)r);
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
__global__ static void norms(long N1, long n3, V u, V w, double *out) {
  __shared__ double s[3 * nthread];
  double v[3] = {0, 0, 0};
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
    }
  }
  blocksum(3, v, s, out);
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
__global__ static void spectrum(long N1, long n3, long nb, V u, double *E) {
  long t = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (t >= 2 * n3)
    return;
  int c = t / n3;
  long l = t % n3, g[3] = {wavenumber(c, l / (N1 * N1)),
                           wavenumber(c, l / N1 % N1), wavenumber(c, l % N1)};
  long b = (long)(2 * sqrt((double)(g[0] * g[0] + g[1] * g[1] + g[2] * g[2])));
  double e = 0;
  for (int d = 0; d < 3; d++) {
    double su = 1;
    for (int q = 0; q < 3; q++)
      su *= q != d && g[q] == 0 ? 1 : 0.5;
    e += su * u.a[c][d][l] * u.a[c][d][l] / 2;
  }
  if (b < nb && e != 0)
    atomicAdd(&E[b], e);
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
__global__ static void advance(int c, int first, long N1, long n3, long kcut,
                               double nu, double dt, V du, V u, V um) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n3)
    return;
  long i = l / (N1 * N1), j = l / N1 % N1, k = l % N1;
  double m = wavenumber(c, i), n = wavenumber(c, j), p = wavenumber(c, k),
         kk = m * m + n * n + p * p, P, cv, v, f[3];
  for (int d = 0; d < 3; d++)
    f[d] = du.a[c][d][l];
  if (i > kcut || j > kcut || k > kcut)
    f[0] = f[1] = f[2] = 0;
  P = kk > 0 ? (m * f[0] + n * f[1] + p * f[2]) / kk : 0;
  f[0] -= m * P;
  f[1] -= n * P;
  f[2] -= p * P;
  for (int d = 0; d < 3; d++) {
    v = u.a[c][d][l];
    if (first) {
      cv = nu * dt * kk / 2;
      u.a[c][d][l] = ((1 - cv) * v + dt * f[d]) / (1 + cv);
    } else {
      cv = nu * dt * kk;
      u.a[c][d][l] = ((1 - cv) * um.a[c][d][l] + 2 * dt * f[d]) / (1 + cv);
    }
    um.a[c][d][l] = v;
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
static long steps(const char *opt, double x, double dt) {
  long n = lround(x / dt);
  if (fabs(n * dt - x) > 1e-9 * x) {
    fprintf(stderr, "tg: error: %s %g is not a multiple of -s %g\n", opt, x,
            dt);
    exit(1);
  }
  return n;
}
int main(int argc, char **argv) {
  long M, tstep, ne, nd, kcut, nt;
  int c, d;
  double nu, dt, T_end, t, x, e, dd, *F[3], *tmp, *dP, *P, *dE, *host, one;
  size_t ws;
  void *work;
  V u, um, w, du, U, W;
  char *end;
  static const int tu[3][3] = {
      {SIN, COS, COS}, {COS, SIN, COS}, {COS, COS, SIN}};
  static const int tw[3][3] = {
      {COS, SIN, SIN}, {SIN, COS, SIN}, {SIN, SIN, COS}};
  (void)argc;
  M = 0;
  nu = -1;
  dt = -1;
  T_end = 0;
  e = 0;
  dd = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    if (argv[0][1] == 'h') {
      fprintf(stderr, "Usage: tg -M <modes> -n <viscosity> -t <end time> -s "
                      "<time step> [-e <spectrum interval>] [-d <dump interval>]\n");
      exit(1);
    }
    if (argv[1] == NULL) {
      fprintf(stderr, "tg: error: %s needs an argument\n", argv[0]);
      exit(1);
    }
    x = strtod(argv[1], &end);
    if (*end != '\0') {
      fprintf(stderr, "tg: error: '%s' is not a number\n", argv[1]);
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
      fprintf(stderr, "tg: error: unknown option '%s'\n", *argv);
      exit(1);
    }
    argv++;
  }
  if (M < 4 || M % 2 != 0 || nu < 0 || dt <= 0 || T_end <= 0) {
    fprintf(stderr, "tg: error: need -M (even, >= 4), -n, -s, -t\n");
    exit(1);
  }
  N = M / 2;
  N1 = N + 1;
  n3 = N1 * N1 * N1;
  B = N1 * N1;
  kcut = M / 3;
  ne = e > 0 ? steps("-e", e, dt) : 0;
  nd = dd > 0 ? steps("-d", dd, dt) : 0;
  nt = steps("-t", T_end, dt);
  for (c = 0; c < 2; c++)
    for (d = 0; d < 3; d++) {
      cuda(cudaMalloc(&u.a[c][d], n3 * sizeof(double)), "cudaMalloc");
      cuda(cudaMalloc(&um.a[c][d], n3 * sizeof(double)), "cudaMalloc");
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
    fprintf(stderr, "tg: error: malloc failed\n");
    exit(1);
  }
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
      norms<<<nblock, nthread>>>(N1, n3, u, w, dP);
      launched("norms");
      collect(3, dP, P, sums);
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
      printf("\n");
      if (fflush(stdout) != 0 || ferror(stdout)) {
        fprintf(stderr, "tg: error: fail to write stdout\n");
        exit(1);
      }
    }
    if (ne > 0 && tstep % ne == 0) {
      long nb = 2 * (long)(sqrt(3.0) * (M + 1)) + 2;
      double *E;
      char path[FILENAME_MAX];
      FILE *file;
      if ((E = (double *)malloc(nb * sizeof(double))) == NULL) {
        fprintf(stderr, "tg: error: malloc failed\n");
        exit(1);
      }
      cuda(cudaMalloc(&dE, nb * sizeof(double)), "cudaMalloc");
      cuda(cudaMemset(dE, 0, nb * sizeof(double)), "cudaMemset");
      spectrum<<<grid(2 * n3), nthread>>>(N1, n3, nb, u, dE);
      launched("spectrum");
      cuda(cudaMemcpy(E, dE, nb * sizeof(double), cudaMemcpyDeviceToHost),
           "cudaMemcpy");
      cuda(cudaFree(dE), "cudaFree");
      sprintf(path, "e.%08ld", tstep);
      if ((file = fopen(path, "w")) == NULL) {
        fprintf(stderr, "tg: error: fail to open '%s'\n", path);
        exit(1);
      }
      fprintf(file, "# t = %.16e\n", t);
      for (long b = 0; b < nb; b++)
        fprintf(file, "%.1f %.16e\n", b / 2.0, E[b]);
      if (ferror(file) || fclose(file) != 0) {
        fprintf(stderr, "tg: error: fail to write '%s'\n", path);
        exit(1);
      }
      free(E);
    }
    if (nd > 0 && tstep % nd == 0) {
      char path[FILENAME_MAX];
      FILE *file;
      if (host == NULL &&
          (host = (double *)malloc(n3 * sizeof(double))) == NULL) {
        fprintf(stderr, "tg: error: malloc failed\n");
        exit(1);
      }
      sprintf(path, "u.%08ld", tstep);
      if ((file = fopen(path, "w")) == NULL) {
        fprintf(stderr, "tg: error: fail to open '%s'\n", path);
        exit(1);
      }
      for (c = 0; c < 2; c++)
        for (d = 0; d < 3; d++) {
          cuda(cudaMemcpy(host, u.a[c][d], n3 * sizeof(double),
                          cudaMemcpyDeviceToHost),
               "cudaMemcpy");
          if (fwrite(host, sizeof(double), n3, file) != (size_t)n3) {
            fprintf(stderr, "tg: error: fail to write '%s'\n", path);
            exit(1);
          }
        }
      if (ferror(file) || fclose(file) != 0) {
        fprintf(stderr, "tg: error: fail to write '%s'\n", path);
        exit(1);
      }
    }
    if (tstep >= nt)
      break;
    for (c = 0; c < 2; c++)
      for (d = 0; d < 3; d++) {
        synthesis(c, tu[d], u.a[c][d], tmp, U.a[c][d]);
        synthesis(c, tw[d], w.a[c][d], tmp, W.a[c][d]);
      }
    for (c = 0; c < 2; c++) {
      cross<<<grid(n3), nthread>>>(c, n3, U, W, F[0], F[1], F[2]);
      launched("cross");
      for (d = 0; d < 3; d++)
        analysis(c, tu[d], F[d], du.a[c][d]);
    }
    for (c = 0; c < 2; c++) {
      advance<<<grid(n3), nthread>>>(c, tstep == 0, N1, n3, kcut, nu, dt, du,
                                     u, um);
      launched("advance");
    }
    tstep++;
    t = tstep * dt;
  }
  cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
}
