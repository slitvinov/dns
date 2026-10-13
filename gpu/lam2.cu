// full-domain lambda2 (float32) of a box dump u.STEP (3 x nx x ny x nz doubles), as lam2g.py:
// gradients by FFT (cuFFT, double), eigenvalues of S^2 + O^2 in closed form (double);
// writes out.npy, out_xyz.npy (x0 y0 z0 dx dy dz) and out_ctr.npy (pair centre x y: circular
// mean of wz^2 over x and y). Domain 6pi x 4pi x 2pi. Build: nvcc -O3 -arch=sm_90 -o lam2 lam2.cu -lcufft
#include <cufft.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { nthread = 256 };
static const double pi = 3.141592653589793238;
static void cuda(cudaError_t e, const char *what) {
  if (e != cudaSuccess) {
    fprintf(stderr, "lam2: error: %s: %s\n", what, cudaGetErrorString(e));
    exit(1);
  }
}
static void cufft(cufftResult r, const char *what) {
  if (r != CUFFT_SUCCESS) {
    fprintf(stderr, "lam2: error: %s: cufft error %d\n", what, (int)r);
    exit(1);
  }
}
static unsigned grid(long n) { return (n + nthread - 1) / nthread; }

// h is the D2Z spectrum (nx, ny, nzf); d = i k_j h / N (torch.fft.fftfreq signs: Nyquist negative in x, y)
__global__ static void deriv(long nx, long ny, long nzf, double scale, double kx0, double ky0,
                             double kz0, int j, const double2 *h, double2 *d) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, n = nx * ny * nzf;
  if (l >= n)
    return;
  long k = l % nzf, y = l / nzf % ny, x = l / nzf / ny;
  double kk = j == 0 ? (x < nx / 2 ? x : x - nx) * kx0
              : j == 1 ? (y < ny / 2 ? y : y - ny) * ky0
                       : k * kz0;
  double2 a = h[l];
  d[l] = make_double2(-kk * a.y * scale, kk * a.x * scale);
}
// padded real (nx, ny, 2 nzf) -> float (nx, ny, nz)
__global__ static void tofloat(long nx, long ny, long nz, long nzf, const double *r, float *f) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= nx * ny * nz)
    return;
  long z = l % nz, xy = l / nz;
  f[l] = (float)r[xy * 2 * nzf + z];
}
struct Grad {
  const float *g[3][3];
};
__global__ static void lam2(long n, Grad G, float *lam, double *wz2) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n)
    return;
  double g[3][3], s[3][3], o[3][3], M[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g[i][j] = G.g[i][j][l];
  wz2[l] = (g[1][0] - g[0][1]) * (g[1][0] - g[0][1]);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      s[i][j] = 0.5 * (g[i][j] + g[j][i]);
      o[i][j] = 0.5 * (g[i][j] - g[j][i]);
    }
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      M[i][j] = 0;
      for (int k = 0; k < 3; k++)
        M[i][j] += s[i][k] * s[k][j] + o[i][k] * o[k][j];
    }
  for (int sw = 0; sw < 8; sw++)
    for (int pq = 0; pq < 3; pq++) {
      int a = pq == 2 ? 1 : 0, b = pq == 0 ? 1 : 2, c = 3 - a - b;
      double apq = M[a][b];
      if (apq == 0)
        continue;
      double th = (M[b][b] - M[a][a]) / (2 * apq);
      double t = fabs(th) > 1e150 ? 0.5 / th : copysign(1.0, th) / (fabs(th) + sqrt(th * th + 1));
      double co = 1 / sqrt(t * t + 1), si = t * co;
      double mac = M[a][c], mbc = M[b][c];
      M[a][a] -= t * apq;
      M[b][b] += t * apq;
      M[a][b] = M[b][a] = 0;
      M[a][c] = M[c][a] = co * mac - si * mbc;
      M[b][c] = M[c][b] = si * mac + co * mbc;
    }
  double e0 = M[0][0], e1 = M[1][1], e2 = M[2][2];
  lam[l] = (float)fmax(fmin(e0, e1), fmin(fmax(e0, e1), e2));
}
// w[a] = sum of wz2 over the other two indices; one block per a (deterministic)
__global__ static void sumx(long nx, long m, const double *wz2, double *w) {
  __shared__ double sh[nthread];
  long a = blockIdx.x;
  double acc = 0;
  for (long l = threadIdx.x; l < m; l += blockDim.x)
    acc += wz2[a * m + l];
  sh[threadIdx.x] = acc;
  __syncthreads();
  for (int s = blockDim.x / 2; s > 0; s /= 2) {
    if (threadIdx.x < s)
      sh[threadIdx.x] += sh[threadIdx.x + s];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    w[a] = sh[0];
}
__global__ static void sumy(long nx, long ny, long nz, const double *wz2, double *w) {
  __shared__ double sh[nthread];
  long b = blockIdx.x;
  double acc = 0;
  for (long l = threadIdx.x; l < nx * nz; l += blockDim.x)
    acc += wz2[((l / nz) * ny + b) * nz + l % nz];
  sh[threadIdx.x] = acc;
  __syncthreads();
  for (int s = blockDim.x / 2; s > 0; s /= 2) {
    if (threadIdx.x < s)
      sh[threadIdx.x] += sh[threadIdx.x + s];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    w[b] = sh[0];
}
static void npy(const char *path, const char *descr, const char *shape, const void *data,
                size_t size) {
  char h[128];
  FILE *f;
  int n = snprintf(h, sizeof h, "{'descr': '%s', 'fortran_order': False, 'shape': %s, }",
                   descr, shape);
  int tot = (10 + n + 1 + 63) / 64 * 64;
  memset(h + n, ' ', tot - 10 - n - 1);
  h[tot - 10 - 1] = '\n';
  unsigned short hl = tot - 10;
  if ((f = fopen(path, "w")) == NULL) {
    fprintf(stderr, "lam2: error: fail to open '%s'\n", path);
    exit(1);
  }
  fwrite("\x93NUMPY\x01\x00", 1, 8, f);
  fwrite(&hl, 2, 1, f);
  fwrite(h, 1, hl, f);
  if (fwrite(data, 1, size, f) != size || ferror(f) || fclose(f) != 0) {
    fprintf(stderr, "lam2: error: fail to write '%s'\n", path);
    exit(1);
  }
}
static double centre(const double *w, long m, double L) {
  double sn = 0, cs = 0;
  for (long i = 0; i < m; i++) {
    double a = 2 * pi * i / m;
    sn += w[i] * sin(a);
    cs += w[i] * cos(a);
  }
  double f = fmod(atan2(sn, cs), 2 * pi);
  if (f < 0)
    f += 2 * pi;
  long i = (long)nearbyint(f / (2 * pi) * m) % m;
  return -L / 2 + L / m * i;
}
int main(int argc, char **argv) {
  if (argc != 6) {
    fprintf(stderr, "usage: lam2 u.STEP out.npy nx ny nz\n");
    return 1;
  }
  const char *in = argv[1], *out = argv[2];
  long nx = atol(argv[3]), ny = atol(argv[4]), nz = atol(argv[5]), nzf = nz / 2 + 1;
  long n = nx * ny * nz, nf = nx * ny * nzf;
  double L[3] = {6 * pi, 4 * pi, 2 * pi};
  double *host, *u, *wz2, *wx, *wy;
  double2 *h, *r;
  float *G[3][3], *lam;
  cufftHandle fplan, bplan;
  FILE *f;
  if (nx <= 0 || ny <= 0 || nz <= 0 || strlen(out) < 5 || strcmp(out + strlen(out) - 4, ".npy")) {
    fprintf(stderr, "lam2: error: bad arguments\n");
    return 1;
  }
  cuda(cudaMallocHost(&host, n * sizeof(double)), "cudaMallocHost");
  cuda(cudaMalloc(&u, n * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&r, nf * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&h, nf * sizeof(double2)), "cudaMalloc");
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      cuda(cudaMalloc(&G[i][j], n * sizeof(float)), "cudaMalloc");
  cufft(cufftPlan3d(&fplan, nx, ny, nz, CUFFT_D2Z), "cufftPlan3d");
  cufft(cufftPlan3d(&bplan, nx, ny, nz, CUFFT_Z2D), "cufftPlan3d");
  if ((f = fopen(in, "r")) == NULL) {
    fprintf(stderr, "lam2: error: fail to open '%s'\n", in);
    return 1;
  }
  for (int i = 0; i < 3; i++) {
    if (fread(host, sizeof(double), n, f) != (size_t)n) {
      fprintf(stderr, "lam2: error: fail to read '%s' (wrong nx ny nz?)\n", in);
      return 1;
    }
    cuda(cudaMemcpy(u, host, n * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    cufft(cufftExecD2Z(fplan, u, h), "cufftExecD2Z");
    for (int j = 0; j < 3; j++) {
      deriv<<<grid(nf), nthread>>>(nx, ny, nzf, 1.0 / n, 2 * pi / L[0], 2 * pi / L[1],
                                   2 * pi / L[2], j, h, r);
      cufft(cufftExecZ2D(bplan, r, (double *)r), "cufftExecZ2D");
      tofloat<<<grid(n), nthread>>>(nx, ny, nz, nzf, (double *)r, G[i][j]);
    }
  }
  if (fgetc(f) != EOF) {
    fprintf(stderr, "lam2: error: '%s' is too long (wrong nx ny nz?)\n", in);
    return 1;
  }
  fclose(f);
  cufftDestroy(fplan);
  cufftDestroy(bplan);
  cudaFree(h);
  wz2 = u;
  lam = (float *)r;
  Grad g;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g.g[i][j] = G[i][j];
  lam2<<<grid(n), nthread>>>(n, g, lam, wz2);
  cuda(cudaGetLastError(), "lam2");
  cuda(cudaMalloc(&wx, (nx + ny) * sizeof(double)), "cudaMalloc");
  wy = wx + nx;
  sumx<<<nx, nthread>>>(nx, ny * nz, wz2, wx);
  sumy<<<ny, nthread>>>(nx, ny, nz, wz2, wy);
  cuda(cudaGetLastError(), "sum");
  double w[nx + ny], c[2], xyz[6] = {-L[0] / 2, -L[1] / 2, 0, L[0] / nx, L[1] / ny, L[2] / nz};
  cuda(cudaMemcpy(w, wx, (nx + ny) * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
  c[0] = centre(w, nx, L[0]);
  c[1] = centre(w + nx, ny, L[1]);
  float *hl = (float *)host;
  cuda(cudaMemcpy(hl, lam, n * sizeof(float), cudaMemcpyDeviceToHost), "cudaMemcpy");
  float mn = hl[0], mx = hl[0];
  for (long l = 1; l < n; l++) {
    mn = fminf(mn, hl[l]);
    mx = fmaxf(mx, hl[l]);
  }
  printf("lam2: %s lambda2 min %.6g max %.6g ctr %.6g %.6g\n", out, mn, mx, c[0], c[1]);
  char shape[64], path[4096];
  snprintf(shape, sizeof shape, "(%ld, %ld, %ld)", nx, ny, nz);
  npy(out, "<f4", shape, hl, n * sizeof(float));
  size_t b = strlen(out) - 4;
  snprintf(path, sizeof path, "%.*s_ctr.npy", (int)b, out);
  npy(path, "<f8", "(2,)", c, sizeof c);
  snprintf(path, sizeof path, "%.*s_xyz.npy", (int)b, out);
  npy(path, "<f8", "(6,)", xyz, sizeof xyz);
  return 0;
}
