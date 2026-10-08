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
struct G {
  long nx, ny, nz, nzf, nf, nr;
  double ax, ay, az;
};
static void cuda(cudaError_t e, const char *what) {
  if (e != cudaSuccess) {
    fprintf(stderr, "box: error: %s: %s\n", what, cudaGetErrorString(e));
    exit(1);
  }
}
static void cufft(cufftResult r, const char *what) {
  if (r != CUFFT_SUCCESS) {
    fprintf(stderr, "box: error: %s: cufft error %d\n", what, (int)r);
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
__device__ static int kept(G g, long i, long j, long k) {
  return 3 * labs(wave(i, g.nx)) < g.nx && 3 * labs(wave(j, g.ny)) < g.ny &&
         3 * k < g.nz;
}
__global__ static void spread(G g, const double2 *u, const double2 *v,
                              const double2 *w, double2 *b0, double2 *b1,
                              double2 *b2, double2 *b3, double2 *b4,
                              double2 *b5) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double2 a = u[l], b = v[l], c = w[l];
  b0[l] = a;
  b1[l] = b;
  b2[l] = c;
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
__global__ static void update(G g, double nu, double dt, double a, double b,
                              const double2 *f0, const double2 *f1,
                              const double2 *f2, double2 *u, double2 *v,
                              double2 *w, double2 *du, double2 *dv,
                              double2 *dw) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double kk = kx * kx + ky * ky + kz * kz, m = kept(g, i, j, k), visc = nu * kk;
  double2 x = f0[l], y = f1[l], z = f2[l], p;
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
  u[l] = make_double2(U.x + b * dx.x, U.y + b * dx.y);
  v[l] = make_double2(V.x + b * dy.x, V.y + b * dy.y);
  w[l] = make_double2(W.x + b * dz.x, W.y + b * dz.y);
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
__global__ static void norms(G g, const double2 *u, const double2 *v,
                             const double2 *w, double *out) {
  __shared__ double s[3 * nthread];
  double e[3] = {0, 0, 0};
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
  }
  blockreduce(3, 0, e, s, out);
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
    fprintf(stderr, "box: error: %s %g is not a multiple of -s %g\n", opt, x,
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
  double nu, dt, T, t, dd, X, Y, Z, x, N, s, cx, cy, cz, kc, sums[3], mx[2],
      *host, *dP, *P;
  double2 *u[3], *du[3], *b[6];
  long tstep, start, nd, nx, ny, nz, nt;
  size_t wf, wb;
  void *work;
  (void)argc;
  ipath = rpath = NULL;
  nx = ny = nz = 0;
  X = Y = Z = 2;
  nu = dt = -1;
  T = dd = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    if (argv[0][1] == 'h') {
      fprintf(stderr,
              "Usage: box -x <nx> -y <ny> -z <nz> [-X <Lx/pi>] [-Y <Ly/pi>] "
              "[-Z <Lz/pi>] (-i <u.raw> | -r <dump u.STEP>) -n <viscosity> "
              "-t <end time> -s <time step> [-d <dump interval>]\n");
      exit(1);
    }
    if (argv[1] == NULL) {
      fprintf(stderr, "box: error: %s needs an argument\n", argv[0]);
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
    x = strtod(argv[1], &end);
    if (*end != '\0') {
      fprintf(stderr, "box: error: '%s' is not a number\n", argv[1]);
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
    default:
      fprintf(stderr, "box: error: unknown option '%s'\n", *argv);
      exit(1);
    }
    argv++;
  }
  if (nx < 4 || ny < 4 || nz < 4 || nz % 2 != 0 || nu < 0 || dt <= 0 ||
      T <= 0 || X <= 0 || Y <= 0 || Z <= 0 || (ipath == NULL) == (rpath == NULL)) {
    fprintf(stderr, "box: error: need -x -y -z (nz even), -n, -s, -t and one "
                    "of -i, -r\n");
    exit(1);
  }
  g.nx = nx;
  g.ny = ny;
  g.nz = nz;
  g.nzf = nz / 2 + 1;
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
  cuda(cudaMalloc(&dP, 3 * nblock * sizeof(double)), "cudaMalloc");
  if ((P = (double *)malloc(3 * nblock * sizeof(double))) == NULL ||
      (host = (double *)malloc(nx * ny * nz * sizeof(double))) == NULL) {
    fprintf(stderr, "box: error: malloc failed\n");
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
  if (rpath != NULL) {
    const char *p = strrchr(rpath, '.');
    if (p == NULL || (tstep = strtol(p + 1, &end, 10), *end != '\0' || end == p + 1)) {
      fprintf(stderr, "box: error: '%s' is not named u.STEP\n", rpath);
      exit(1);
    }
    name = rpath;
  }
  if ((file = fopen(name, "r")) == NULL) {
    fprintf(stderr, "box: error: fail to open '%s'\n", name);
    exit(1);
  }
  for (int c = 0; c < 3; c++) {
    if (fread(host, sizeof(double), nx * ny * nz, file) != (size_t)(nx * ny * nz)) {
      fprintf(stderr, "box: error: fail to read '%s' (wrong -x -y -z?)\n", name);
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
    fprintf(stderr, "box: error: '%s' is too long (wrong -x -y -z?)\n", name);
    exit(1);
  }
  fclose(file);
  start = tstep;
  t = tstep * dt;
  for (;;) {
    int diag = tstep % 10 == 0, dump = nd > 0 && tstep % nd == 0 && (rpath == NULL || tstep != start);
    if (diag || dump) {
      spread<<<grid(g.nf), nthread>>>(g, u[0], u[1], u[2], b[0], b[1], b[2], b[3],
                                      b[4], b[5]);
      launched("spread");
      for (int c = 0; c < 6; c++)
        cufft(cufftExecZ2D(bplan, b[c], (double *)b[c]), "cufftExecZ2D");
    }
    if (diag) {
      norms<<<nblock, nthread>>>(g, u[0], u[1], u[2], dP);
      launched("norms");
      collect(3, 0, dP, P, sums);
      peaks<<<nblock, nthread>>>(g, s, cx, cy, cz, (double *)b[0], (double *)b[1],
                                 (double *)b[2], (double *)b[3], (double *)b[4],
                                 (double *)b[5], dP);
      launched("peaks");
      collect(2, 1, dP, P, mx);
      double E = sums[0] * s * s, Om = sums[1] * s * s, eps = 2 * nu * Om,
             u2 = 2 * E / 3;
      double q[] = {t,
                    E,
                    Om,
                    dt * mx[0],
                    sqrt(mx[1]),
                    kc * pow(nu * nu * nu / eps, 0.25),
                    u2 * sqrt(15 / (nu * eps)),
                    sqrt(10 * nu * E / eps),
                    pi / (2 * u2) * sums[2] * s * s};
      if (tstep == start)
        printf("step t E Omega C wmax keta Rlambda lambda L\n");
      printf("% 10ld", tstep);
      for (int i = 0; i < (int)(sizeof q / sizeof *q); i++)
        printf(" % .16e", q[i]);
      printf("\n");
      if (fflush(stdout) != 0 || ferror(stdout)) {
        fprintf(stderr, "box: error: fail to write stdout\n");
        exit(1);
      }
    }
    if (dump) {
      sprintf(path, "u.%08ld", tstep);
      if ((file = fopen(path, "w")) == NULL) {
        fprintf(stderr, "box: error: fail to open '%s'\n", path);
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
          fprintf(stderr, "box: error: fail to write '%s'\n", path);
          exit(1);
        }
      }
      if (ferror(file) || fclose(file) != 0) {
        fprintf(stderr, "box: error: fail to write '%s'\n", path);
        exit(1);
      }
    }
    if (tstep >= nt)
      break;
    for (int r = 0; r < 5; r++) {
      spread<<<grid(g.nf), nthread>>>(g, u[0], u[1], u[2], b[0], b[1], b[2], b[3],
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
      update<<<grid(g.nf), nthread>>>(g, nu, dt, A[r], B[r], b[0], b[1], b[2],
                                      u[0], u[1], u[2], du[0], du[1], du[2]);
      launched("update");
    }
    tstep++;
    t = tstep * dt;
  }
  cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
}
