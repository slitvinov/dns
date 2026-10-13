/* boxsymmp: MPI multi-GPU version of gpu/boxsym.cu (box restricted to the Sx, Sz symmetric fields of the
   colliding tubes: u(-x,y,z) = (-u,v,w)(x,y,z), u(x,y,-z) = (u,v,-w)(x,y,z)). Same options and out lines
   as boxsym; reads and writes only q.STEP (double [3][nx/2+1][nz/2+1][ny]; gpu/symq makes one from a full
   box u.STEP, gpu/symfull expands one back).
   Decomposition (X1 = nx/2+1, Z1 = nz/2+1, Q = ny/2+1, P ranks, one GPU each): physical fields and the
   x-stage in z-slabs [i][kloc][j], spectral fields and the z-stage in x-slabs [ploc][r][q]; one
   all-to-all per 3-D transform between the x- and z-stages. y is R2C/C2R, x and z are full-length
   complex FFTs of parity-extended lines (as in boxsym, tiled transposes). The all-to-all is NCCL grouped
   send/recv on its own stream, pipelined against the compute of the neighbouring fields (MPI only
   bootstraps NCCL and reduces the diagnostics). */
#include <cufft.h>
#include <fcntl.h>
#include <math.h>
#include <mpi.h>
#include <nccl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

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
static const int ux[] = {1, 0, 0}, uz[] = {0, 0, 1}, wx[] = {0, 1, 1}, wz[] = {1, 1, 0};
/* parities of the six synthesized fields u, omega */
static const int ox6[] = {1, 0, 0, 0, 1, 1}, oz6[] = {0, 0, 1, 1, 1, 0};
struct G {
  long nx, ny, nz, Nx, Nz, X1, Z1, Q, nf, nr;
  long P, p0, pc, k0, kc; /* ranks, my x-slab p0+pc (spectral), my z-slab k0+kc (physical) */
  double ax, ay, az;
  int f;
  double fp;
};
static G g;
static int rank, size;
static cufftHandle px, pz, py[2];
static double2 *R, *T[2];
static int *scz, *sdz, *rcz, *rdz; /* alltoallv: z-slab sections <-> x-slab chunks */
static ncclComm_t nccl;
static int usempi; /* BOXSYM_MPI=1: CUDA-aware MPI_Alltoallv instead of NCCL (e.g. several ranks per GPU) */
static MPI_Datatype cplx;
static cudaStream_t st, sc; /* compute, all-to-all */
static cudaEvent_t ev[4][8];   /* pipeline events: [kind][field] */
static void fail(const char *fmt, const char *a, const char *b) {
  fprintf(stderr, "boxsymmp: error: rank %d: ", rank);
  fprintf(stderr, fmt, a, b);
  fprintf(stderr, "\n");
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}
static void cuda(cudaError_t e, const char *what) {
  if (e != cudaSuccess)
    fail("%s: %s", what, cudaGetErrorString(e));
}
static void cufft(cufftResult r, const char *what) {
  char c[32];
  if (r != CUFFT_SUCCESS) {
    snprintf(c, sizeof c, "cufft error %d", (int)r);
    fail("%s: %s", what, c);
  }
}
static void launched(const char *what) { cuda(cudaGetLastError(), what); }
static void nccheck(ncclResult_t r, const char *what) {
  if (r != ncclSuccess)
    fail("%s: %s", what, ncclGetErrorString(r));
}
static void devsync(void) { cuda(cudaStreamSynchronize(st), "cudaStreamSynchronize"); }
/* all-to-all of complex chunks: to rank q send sc[q] at sd[q], receive rc[q] at rd[q] */
static void a2a(const double2 *sb, const int *nsend, const int *sd, double2 *rb, const int *rc,
                const int *rd) {
  if (usempi) {
    devsync();
    MPI_Alltoallv(sb, nsend, sd, cplx, rb, rc, rd, cplx, MPI_COMM_WORLD);
    return;
  }
  if (rc[rank] > 0)
    cuda(cudaMemcpyAsync(rb + rd[rank], sb + sd[rank], rc[rank] * sizeof(double2),
                         cudaMemcpyDeviceToDevice, sc),
         "cudaMemcpyAsync");
  nccheck(ncclGroupStart(), "ncclGroupStart");
  for (int q = 0; q < size; q++) {
    if (q == rank)
      continue;
    if (nsend[q] > 0)
      nccheck(ncclSend(sb + sd[q], 2 * (size_t)nsend[q], ncclDouble, q, nccl, sc), "ncclSend");
    if (rc[q] > 0)
      nccheck(ncclRecv(rb + rd[q], 2 * (size_t)rc[q], ncclDouble, q, nccl, sc), "ncclRecv");
  }
  nccheck(ncclGroupEnd(), "ncclGroupEnd");
}
static unsigned grid(long n) { return (n + nthread - 1) / nthread; }
/* balanced split of n into P parts: part s starts at n s / P; owner of h: */
__host__ __device__ static long first(long n, long P, long s) { return n * s / P; }
__host__ __device__ static long owner(long n, long P, long h) { return ((h + 1) * P - 1) / n; }
__device__ static void mode(G g, long l, long *p, long *q, long *r, double *kx,
                            double *ky, double *kz) {
  *p = g.p0 + l / (g.Z1 * g.Q);
  *r = l / g.Q % g.Z1;
  *q = l % g.Q;
  *kx = g.ax * *p;
  *ky = g.ay * *q;
  *kz = g.az * *r;
}
__device__ static double mask(G g, long p, long q, long r) {
  if (g.f) {
    double x = 2.0 * p / g.nx, y = 2.0 * q / g.ny, z = 2.0 * r / g.nz;
    return exp(-36 * (pow(x, 36) + pow(y, 36) + pow(z, 36)));
  }
  return 3 * p < g.nx && 3 * q < g.ny && 3 * r < g.nz;
}
/* parity extension of n/2+1 values to a full period n, the same map for spectral (F(-p) = +-F(p)) and
   physical (u(-x) = +-u(x)) lines (odd lines vanish at 0 and n/2), fused with a transpose so that the
   FFT lines are contiguous: 32x32 tiles through shared memory, blocks of 32x8 threads.
   x: a[h][b] (b < M) <-> R[b][m] (m < n);  z: slab pl of a [pl][r][q] layout (S sections, see zrow)
   <-> R[pl][q][m]. fold* keep m <= n/2. */
enum { TL = 32, TR = 8 };
__device__ static long zrow(long Z1, long S, long pc, long Q, long pl, long h) {
  long s = ((h + 1) * S - 1) / Z1, r0 = Z1 * s / S, cnt = Z1 * (s + 1) / S - r0;
  return pc * r0 * Q + (pl * cnt + h - r0) * Q;
}
__global__ static void extx(long n, int odd, long M, const double2 *a, double2 *R) {
  __shared__ double2 t[TL][TL + 1];
  long b0 = blockIdx.x * (long)TL, m0 = blockIdx.y * (long)TL, N = n / 2;
  int x = threadIdx.x, y = threadIdx.y;
  for (int j = y; j < TL; j += TR) {
    long m = m0 + j, b = b0 + x, h = m <= N ? m : n - m;
    if (m < n && b < M) {
      double2 v = odd && (h == 0 || h == N) ? make_double2(0, 0) : a[h * M + b];
      t[j][x] = odd && m > N ? make_double2(-v.x, -v.y) : v;
    }
  }
  __syncthreads();
  for (int j = y; j < TL; j += TR) {
    long b = b0 + j, m = m0 + x;
    if (m < n && b < M)
      R[b * n + m] = t[x][j];
  }
}
__global__ static void foldx(long n, long M, const double2 *R, double2 *a) {
  __shared__ double2 t[TL][TL + 1];
  long b0 = blockIdx.x * (long)TL, p0 = blockIdx.y * (long)TL, N1 = n / 2 + 1;
  int x = threadIdx.x, y = threadIdx.y;
  for (int j = y; j < TL; j += TR) {
    long b = b0 + j, p = p0 + x;
    if (p < N1 && b < M)
      t[j][x] = R[b * n + p];
  }
  __syncthreads();
  for (int j = y; j < TL; j += TR) {
    long p = p0 + j, b = b0 + x;
    if (p < N1 && b < M)
      a[p * M + b] = t[x][j];
  }
}
__global__ static void extz(long n, int odd, long Q, long Z1, long pc, long S, const double2 *a,
                            double2 *R) {
  __shared__ double2 t[TL][TL + 1];
  long q0 = blockIdx.x * (long)TL, m0 = blockIdx.y * (long)TL, pl = blockIdx.z, N = n / 2;
  int x = threadIdx.x, y = threadIdx.y;
  for (int j = y; j < TL; j += TR) {
    long m = m0 + j, q = q0 + x, h = m <= N ? m : n - m;
    if (m < n && q < Q) {
      double2 v = odd && (h == 0 || h == N) ? make_double2(0, 0) : a[zrow(Z1, S, pc, Q, pl, h) + q];
      t[j][x] = odd && m > N ? make_double2(-v.x, -v.y) : v;
    }
  }
  __syncthreads();
  for (int j = y; j < TL; j += TR) {
    long q = q0 + j, m = m0 + x;
    if (m < n && q < Q)
      R[(pl * Q + q) * n + m] = t[x][j];
  }
}
__global__ static void foldz(long n, long Q, long Z1, long pc, long S, const double2 *R,
                             double2 *a) {
  __shared__ double2 t[TL][TL + 1];
  long q0 = blockIdx.x * (long)TL, k0 = blockIdx.y * (long)TL, pl = blockIdx.z;
  int x = threadIdx.x, y = threadIdx.y;
  for (int j = y; j < TL; j += TR) {
    long q = q0 + j, k = k0 + x;
    if (k < Z1 && q < Q)
      t[j][x] = R[(pl * Q + q) * n + k];
  }
  __syncthreads();
  for (int j = y; j < TL; j += TR) {
    long k = k0 + j, q = q0 + x;
    if (k < Z1 && q < Q)
      a[zrow(Z1, S, pc, Q, pl, k) + q] = t[x][j];
  }
}
static dim3 tiles(long a, long b, long c) { return dim3((a + TL - 1) / TL, (b + TL - 1) / TL, c); }
static void rec(int k, int c, cudaStream_t s) { cuda(cudaEventRecord(ev[k][c], s), "cudaEventRecord"); }
static void waitev(cudaStream_t s, int k, int c) {
  cuda(cudaStreamWaitEvent(s, ev[k][c], 0), "cudaStreamWaitEvent");
}
/* n fields a[c] with parities ox[c], oz[c]. Pipelined: the all-to-all of field c (stream sc, buffer
   T[c % 2]) overlaps the compute of fields c - 1 and c + 1 (stream st). Events: 0 compute side done
   with T[c % 2] / a[c] ready to send, 1 all-to-all c done, 2 T[c % 2] free again. */
static void synthesis(int n, double2 **a, const int *ox, const int *oz) {
  for (int c = 0; c <= n; c++) {
    if (c < n) {
      extz<<<tiles(g.Q, g.nz, g.pc), dim3(TL, TR), 0, st>>>(g.nz, oz[c], g.Q, g.Z1, g.pc, 1, a[c], R);
      launched("extz");
      cufft(cufftExecZ2Z(pz, R, R, CUFFT_INVERSE), "cufftExecZ2Z");
      if (c >= 2)
        waitev(st, 1, c - 2);
      foldz<<<tiles(g.Q, g.Z1, g.pc), dim3(TL, TR), 0, st>>>(g.nz, g.Q, g.Z1, g.pc, size, R,
                                                             T[c % 2]);
      launched("foldz");
      rec(0, c, st);
      waitev(sc, 0, c);
      a2a(T[c % 2], scz, sdz, a[c], rcz, rdz);
      rec(1, c, sc);
    }
    if (c >= 1) {
      int d = c - 1;
      waitev(st, 1, d);
      extx<<<tiles(g.kc * g.Q, g.nx, 1), dim3(TL, TR), 0, st>>>(g.nx, ox[d], g.kc * g.Q, a[d], R);
      launched("extx");
      cufft(cufftExecZ2Z(px, R, R, CUFFT_INVERSE), "cufftExecZ2Z");
      foldx<<<tiles(g.kc * g.Q, g.X1, 1), dim3(TL, TR), 0, st>>>(g.nx, g.kc * g.Q, R, a[d]);
      launched("foldx");
      cufft(cufftExecZ2D(py[1], a[d], (double *)a[d]), "cufftExecZ2D");
    }
  }
}
static void analysis(int n, double2 **a, const int *ox, const int *oz) {
  for (int c = 0; c <= n; c++) {
    if (c < n) {
      cufft(cufftExecD2Z(py[0], (double *)a[c], a[c]), "cufftExecD2Z");
      extx<<<tiles(g.kc * g.Q, g.nx, 1), dim3(TL, TR), 0, st>>>(g.nx, ox[c], g.kc * g.Q, a[c], R);
      launched("extx");
      cufft(cufftExecZ2Z(px, R, R, CUFFT_FORWARD), "cufftExecZ2Z");
      foldx<<<tiles(g.kc * g.Q, g.X1, 1), dim3(TL, TR), 0, st>>>(g.nx, g.kc * g.Q, R, a[c]);
      launched("foldx");
      rec(0, c, st);
      waitev(sc, 0, c);
      if (c >= 2)
        waitev(sc, 2, c - 2);
      a2a(a[c], rcz, rdz, T[c % 2], scz, sdz);
      rec(1, c, sc);
    }
    if (c >= 1) {
      int d = c - 1;
      waitev(st, 1, d);
      extz<<<tiles(g.Q, g.nz, g.pc), dim3(TL, TR), 0, st>>>(g.nz, oz[d], g.Q, g.Z1, g.pc, size,
                                                            T[d % 2], R);
      launched("extz");
      rec(2, d, st);
      cufft(cufftExecZ2Z(pz, R, R, CUFFT_FORWARD), "cufftExecZ2Z");
      foldz<<<tiles(g.Q, g.Z1, g.pc), dim3(TL, TR), 0, st>>>(g.nz, g.Q, g.Z1, g.pc, 1, R, a[d]);
      launched("foldz");
    }
  }
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
  if (l >= g.nr || l % (2 * g.Q) >= g.ny)
    return;
  double u = b0[l] * s, v = b1[l] * s, w = b2[l] * s;
  double p = b3[l] * s, q = b4[l] * s, r = b5[l] * s;
  b0[l] = v * r - w * q;
  b1[l] = w * p - u * r;
  b2[l] = u * q - v * p;
}
__global__ static void update(G g, int last, double nu, double dt, double a, double b,
                              const double2 *f0, const double2 *f1,
                              const double2 *f2, double2 *u, double2 *v,
                              double2 *w, double2 *du, double2 *dv,
                              double2 *dw) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double kk = kx * kx + ky * ky + kz * kz, m = mask(g, i, j, k), visc = nu * kk;
  /* modes that do not exist for an odd component (sin 0, sin N) */
  double ex = i == 0 || i == g.Nx ? 0 : 1, ez = k == 0 || k == g.Nz ? 0 : 1;
  double2 x = f0[l], y = f1[l], z = f2[l], p;
  x = make_double2(ex * m * x.x, ex * m * x.y);
  y = make_double2(m * y.x, m * y.y);
  z = make_double2(ez * m * z.x, ez * m * z.y);
  p = kk > 0 ? make_double2((kx * x.x + ky * y.x + kz * z.x) / kk,
                            (kx * x.y + ky * y.y + kz * z.y) / kk)
             : make_double2(0, 0);
  double2 U = u[l], V = v[l], W = w[l];
  double2 fx = make_double2(ex * (x.x - kx * p.x - visc * U.x),
                            ex * (x.y - kx * p.y - visc * U.y));
  double2 fy = make_double2(y.x - ky * p.x - visc * V.x,
                            y.y - ky * p.y - visc * V.y);
  double2 fz = make_double2(ez * (z.x - kz * p.x - visc * W.x),
                            ez * (z.y - kz * p.y - visc * W.y));
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
/* box sums over all its modes: p and -p (p = 1..Nx-1) and q, -q (q = 1..ny/2-1) appear twice, its r weights stay */
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
           h = (k == 0 || k == g.Nz ? 0.5 : 1) * (i == 0 || i == g.Nx ? 1 : 2) *
               (j == 0 || 2 * j == g.ny ? 1 : 2);
    double2 a = u[l], b = v[l], c = w[l];
    double q = h * (a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y + c.x * c.x +
                    c.y * c.y);
    e[0] += q;
    e[1] += kk * q;
    if (kk > 0)
      e[2] += q / sqrt(kk);
    if (i > 0.9 * fx)
      e[3] += q;
    if (j > 0.9 * fy)
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
    if (l % (2 * g.Q) >= g.ny)
      continue;
    double c = s * (cx * fabs(b0[l]) + cy * fabs(b1[l]) + cz * fabs(b2[l]));
    double o = s * s * (b3[l] * b3[l] + b4[l] * b4[l] + b5[l] * b5[l]);
    m[0] = fmax(m[0], c);
    m[1] = fmax(m[1], o);
  }
  blockreduce(2, 1, m, sh, out);
}
/* shell spectrum (bin width dk = min(ax, ay, az), nearest integer of |k|/dk) and 1D spectra in kx, ky, kz
   (bins p, q, r); weights as in norms, so the sums are box's over its full mode set */
__global__ static void spectrum(G g, double dk, long ns, const double2 *u, const double2 *v,
                                const double2 *w, double *E) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double h = (k == 0 || k == g.Nz ? 0.5 : 1) * (i == 0 || i == g.Nx ? 1 : 2) *
             (j == 0 || 2 * j == g.ny ? 1 : 2);
  double2 a = u[l], b = v[l], c = w[l];
  double q = h * (a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y + c.x * c.x + c.y * c.y);
  if (q == 0)
    return;
  long sb = (long)(sqrt(kx * kx + ky * ky + kz * kz) / dk + 0.5);
  if (sb < ns)
    atomicAdd(&E[sb], q);
  atomicAdd(&E[ns + i], q);
  atomicAdd(&E[ns + g.nx / 2 + 1 + j], q);
  atomicAdd(&E[ns + g.nx / 2 + 1 + g.ny / 2 + 1 + k], q);
}
/* d = scale * i k_j h (spectral derivative along j) */
__global__ static void lderiv(G g, int j, const double2 *h, double2 *d, double scale) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, p, q, r;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &p, &q, &r, &kx, &ky, &kz);
  double kk = j == 0 ? kx : j == 1 ? ky : kz;
  double2 a = h[l];
  d[l] = make_double2(-kk * a.y * scale, kk * a.x * scale);
}
__global__ static void ltofloat(long m, long ny, long Q, const double *r, float *f) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= m * ny)
    return;
  f[l] = (float)r[l / ny * 2 * Q + l % ny];
}
struct Grad {
  const float *g[3][3];
};
/* omega_z^2 = (d_x v - d_y u)^2 on the slab [i][kloc][j] (padded rows 2Q), weighted by the z multiplicity
   of plane k in the full box (1 at k = 0, Nz, else 2) */
__global__ static void lwz2(G g, const double *dv, const double *du, double *wz2) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= g.X1 * g.kc * g.ny)
    return;
  long row = l / g.ny, k = g.k0 + row % g.kc, r = row * 2 * g.Q + l % g.ny;
  double w = dv[r] - du[r];
  wz2[l] = (k == 0 || k == g.Nz ? 1 : 2) * w * w;
}
/* lambda2 (as boxmp/lam2) */
__global__ static void llam2(long n, Grad G_, float *lam) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n)
    return;
  double g_[3][3], s[3][3], o[3][3], M[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g_[i][j] = G_.g[i][j][l];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      s[i][j] = 0.5 * (g_[i][j] + g_[j][i]);
      o[i][j] = 0.5 * (g_[i][j] - g_[j][i]);
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
/* x profile: w[i] = sum over the slab rows of plane i (m values each) */
__global__ static void lsumx(long m, const double *wz2, double *w) {
  __shared__ double sh[nthread];
  long a = blockIdx.x;
  double acc = 0;
  for (long l = threadIdx.x; l < m; l += blockDim.x)
    acc += wz2[a * m + l];
  sh[threadIdx.x] = acc;
  __syncthreads();
  for (int q = blockDim.x / 2; q > 0; q /= 2) {
    if (threadIdx.x < q)
      sh[threadIdx.x] += sh[threadIdx.x + q];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    w[a] = sh[0];
}
/* y profile: w[j] = sum over i (x multiplicity 1 at i = 0, Nx, else 2) and the slab k of wz2[i][k][j] */
__global__ static void lsumy(G g, const double *wz2, double *w) {
  __shared__ double sh[nthread];
  long b = blockIdx.x;
  double acc = 0;
  for (long l = threadIdx.x; l < g.X1 * g.kc; l += blockDim.x) {
    long i = l / g.kc;
    acc += (i == 0 || i == g.Nx ? 1 : 2) * wz2[l * g.ny + b];
  }
  sh[threadIdx.x] = acc;
  __syncthreads();
  for (int q = blockDim.x / 2; q > 0; q /= 2) {
    if (threadIdx.x < q)
      sh[threadIdx.x] += sh[threadIdx.x + q];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    w[b] = sh[0];
}
static double lcentre(const double *w, long m, double L) {
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
static int npyhead(char *h, const char *descr, const char *shape) {
  int n = 10 + snprintf(h + 10, 118, "{'descr': '%s', 'fortran_order': False, 'shape': %s, }",
                        descr, shape);
  int tot = (n + 1 + 63) / 64 * 64;
  memcpy(h, "\x93NUMPY\x01\x00", 8);
  unsigned short hl = tot - 10;
  memcpy(h + 8, &hl, 2);
  memset(h + n, ' ', tot - n - 1);
  h[tot - 1] = '\n';
  return tot;
}
static int npysmall(const char *path, const double *v, int n) {
  char h[128], shape[16];
  FILE *f;
  snprintf(shape, sizeof shape, "(%d,)", n);
  int hn = npyhead(h, "<f8", shape);
  if ((f = fopen(path, "w")) == NULL)
    return -1;
  if (fwrite(h, 1, hn, f) != (size_t)hn || fwrite(v, sizeof(double), n, f) != (size_t)n ||
      ferror(f) || fclose(f) != 0)
    return -1;
  return 0;
}
static void collect(int K, int op, const double *dP, double *P, double *r) {
  devsync();
  cuda(cudaMemcpy(P, dP, nblock * K * sizeof(double), cudaMemcpyDeviceToHost),
       "cudaMemcpy");
  for (int q = 0; q < K; q++)
    r[q] = 0;
  for (int b = 0; b < nblock; b++)
    for (int q = 0; q < K; q++)
      r[q] = op ? fmax(r[q], P[b * K + q]) : r[q] + P[b * K + q];
  MPI_Allreduce(MPI_IN_PLACE, r, K, MPI_DOUBLE, op ? MPI_MAX : MPI_SUM, MPI_COMM_WORLD);
}
static long steps(const char *opt, double x, double dt) {
  long n = lround(x / dt);
  char s[64];
  if (fabs(n * dt - x) > 1e-9 * x) {
    snprintf(s, sizeof s, "%g is not a multiple of -s %g", x, dt);
    fail("%s %s", opt, s);
  }
  return n;
}
static void plan(cufftHandle *p, long long n, long long ine, long long idist,
                 long long one, long long odist, long long stride, long long batch,
                 cufftType type, size_t *ws) {
  size_t sz;
  cufft(cufftCreate(p), "cufftCreate");
  cufft(cufftSetAutoAllocation(*p, 0), "cufftSetAutoAllocation");
  cufft(cufftMakePlanMany64(*p, 1, &n, &ine, stride, idist, &one, stride, odist, type,
                            batch, &sz),
        "cufftMakePlanMany64");
  if (sz > *ws)
    *ws = sz;
}
static int rw(int w, int fd, void *p, size_t n, off_t o) {
  char *c = (char *)p;
  while (n > 0) {
    ssize_t k = w ? pwrite(fd, c, n, o) : pread(fd, c, n, o);
    if (k <= 0)
      return 1;
    c += k;
    o += k;
    n -= k;
  }
  return 0;
}
static int any(int bad) {
  MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return bad;
}
int main(int argc, char **argv) {
  cufftHandle *all[] = {&px, &pz, &py[0], &py[1]};
  char path[FILENAME_MAX], *ipath, *rpath, *end, *name;
  double ee = 0, LL = 0, L0 = 0, HH = 3.2;
  int L0start = 0; /* -0: lambda2 crop also at the restart step */
  char *tlist = NULL;
  long ne, nL, nT = 0, *Ts = NULL;
  double nu, dt, T_, t, dd, F, X, Y, Z, x, N, s, cx, cy, cz, kc, sums[6], mx[2],
      *host, *dP, *P, t0, t1, tio;
  double2 *u[3], *du[3], *b[6];
  long tstep, start, nd, nx, ny, nz, nt, nloc, nb, nR;
  size_t ws, fr, tot;
  void *work;
  int fd, ndev, local, bad;
  MPI_Comm shm;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &shm);
  MPI_Comm_rank(shm, &local);
  cuda(cudaGetDeviceCount(&ndev), "cudaGetDeviceCount");
  cuda(cudaSetDevice(local % ndev), "cudaSetDevice");
  cuda(cudaStreamCreate(&st), "cudaStreamCreate");
  cuda(cudaStreamCreate(&sc), "cudaStreamCreate");
  for (int k = 0; k < 4; k++)
    for (int c = 0; c < 8; c++)
      cuda(cudaEventCreateWithFlags(&ev[k][c], cudaEventDisableTiming), "cudaEventCreate");
  usempi = getenv("BOXSYM_MPI") != NULL && atoi(getenv("BOXSYM_MPI")) != 0;
  MPI_Type_contiguous(2, MPI_DOUBLE, &cplx);
  MPI_Type_commit(&cplx);
  if (!usempi) {
    ncclUniqueId id;
    if (rank == 0)
      nccheck(ncclGetUniqueId(&id), "ncclGetUniqueId");
    MPI_Bcast(&id, sizeof id, MPI_BYTE, 0, MPI_COMM_WORLD);
    nccheck(ncclCommInitRank(&nccl, size, id, rank), "ncclCommInitRank");
  }
  ipath = rpath = NULL;
  nx = ny = nz = 0;
  X = Y = Z = 2;
  nu = dt = -1;
  T_ = dd = 0;
  g.f = 0;
  F = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    if (argv[0][1] == 'h') {
      if (rank == 0)
        fprintf(stderr,
                "Usage: mpirun -np P boxsymmp -x <nx> -y <ny> -z <nz> [-X <Lx/pi>] [-Y <Ly/pi>] "
                "[-Z <Lz/pi>] (-i <q.raw> | -r <q.STEP>) -n <viscosity> -t <end time> "
                "-s <time step> [-d <dump interval>] [-T t1,t2,... (extra dump times)] [-f] [-F <dt0>]\n"
                "  [-e <spectrum interval> (e.STEP: shell and 1D spectra, as boxmp)]\n"
                "  [-L <lambda2 interval> (l.STEP.npy: full-box float32 crop +-H around the pair centre, as boxmp)]\n"
                "  [-l <first lambda2 time>] [-H <half width, default 3.2>] [-0 (crop also at the restart step)]\n"
                "  = gpu/boxsym on P GPUs; input and dumps in q. format only (gpu/symq, gpu/symfull)\n"
                "  all-to-all by NCCL; BOXSYM_MPI=1: CUDA-aware MPI_Alltoallv (allows several ranks per GPU)\n");
      MPI_Finalize();
      exit(1);
    }
    if (argv[0][1] == 'f') {
      g.f = 1;
      continue;
    }
    if (argv[0][1] == '0') {
      L0start = 1;
      continue;
    }
    if (argv[1] == NULL)
      fail("%s needs an argument%s", argv[0], "");
    if (argv[0][1] == 'i') {
      ipath = *++argv;
      continue;
    }
    if (argv[0][1] == 'r') {
      rpath = *++argv;
      continue;
    }
    if (argv[0][1] == 'T') {
      tlist = *++argv;
      continue;
    }
    x = strtod(argv[1], &end);
    if (*end != '\0')
      fail("'%s' is not a number%s", argv[1], "");
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
      T_ = x;
      break;
    case 'd':
      dd = x;
      break;
    case 'F':
      F = x;
      break;
    case 'e':
      ee = x;
      break;
    case 'L':
      LL = x;
      break;
    case 'l':
      L0 = x;
      break;
    case 'H':
      HH = x;
      break;
    default:
      fail("unknown option '%s'%s", *argv, "");
    }
    argv++;
  }
  if (nx < 4 || ny < 4 || nz < 4 || nx % 2 || ny % 2 || nz % 2 || nu < 0 || dt <= 0 ||
      T_ <= 0 || X <= 0 || Y <= 0 || Z <= 0 || (ipath == NULL) == (rpath == NULL))
    fail("need -x -y -z (all even), -n, -s, -t and one of -i, -r%s%s", "", "");
  g.nx = nx;
  g.ny = ny;
  g.nz = nz;
  g.Nx = nx / 2;
  g.Nz = nz / 2;
  g.X1 = g.Nx + 1;
  g.Z1 = g.Nz + 1;
  g.Q = ny / 2 + 1;
  if (size > g.Z1)
    fail("need nz/2+1 >= number of ranks%s%s", "", "");
  g.P = size;
  g.p0 = first(g.X1, size, rank);
  g.pc = first(g.X1, size, rank + 1) - g.p0;
  g.k0 = first(g.Z1, size, rank);
  g.kc = first(g.Z1, size, rank + 1) - g.k0;
  g.fp = F > 0 ? dt / F : 1;
  g.nf = g.pc * g.Z1 * g.Q;
  g.nr = 2 * g.X1 * g.kc * g.Q;
  g.ax = 2 / X;
  g.ay = 2 / Y;
  g.az = 2 / Z;
  nloc = g.X1 * g.kc * ny;
  nb = g.nf > g.X1 * g.kc * g.Q ? g.nf : g.X1 * g.kc * g.Q;
  nR = nx * g.kc > nz * g.pc ? nx * g.kc : nz * g.pc;
  nR *= g.Q;
  if ((scz = (int *)malloc(4 * size * sizeof(int))) == NULL)
    fail("malloc failed%s%s", "", "");
  sdz = scz + size;
  rcz = sdz + size;
  rdz = rcz + size;
  for (int q = 0; q < size; q++) {
    long k0 = first(g.Z1, size, q), k1 = first(g.Z1, size, q + 1), p0 = first(g.X1, size, q),
         p1 = first(g.X1, size, q + 1);
    if (g.pc * k1 * g.Q > 2147483647L || p1 * g.kc * g.Q > 2147483647L)
      fail("slab too large for MPI int counts: use more ranks%s%s", "", "");
    scz[q] = g.pc * (k1 - k0) * g.Q;
    sdz[q] = g.pc * k0 * g.Q;
    rcz[q] = (p1 - p0) * g.kc * g.Q;
    rdz[q] = p0 * g.kc * g.Q;
  }
  N = (double)nx * ny * nz;
  s = 1 / N;
  cx = g.ax * ((nx - 1) / 3);
  cy = g.ay * ((ny - 1) / 3);
  cz = g.az * ((nz - 1) / 3);
  kc = fmin(cx, fmin(cy, cz));
  nd = dd > 0 ? steps("-d", dd, dt) : 0;
  ne = ee > 0 ? steps("-e", ee, dt) : 0;
  nL = LL > 0 ? steps("-L", LL, dt) : 0;
  if (tlist != NULL) {
    char *c = tlist;
    if ((Ts = (long *)malloc((strlen(tlist) + 1) * sizeof(long))) == NULL)
      fail("malloc failed%s%s", "", "");
    for (;;) {
      double v = strtod(c, &end);
      if (end == c || (*end != ',' && *end != '\0'))
        fail("-T '%s' is not a list of times%s", tlist, "");
      Ts[nT++] = steps("-T", v, dt);
      if (*end == '\0')
        break;
      c = end + 1;
    }
  }
  nt = steps("-t", T_, dt);
  for (int c = 0; c < 3; c++) {
    cuda(cudaMalloc(&u[c], g.nf * sizeof(double2)), "cudaMalloc");
    cuda(cudaMalloc(&du[c], g.nf * sizeof(double2)), "cudaMalloc");
    cuda(cudaMemset(du[c], 0, g.nf * sizeof(double2)), "cudaMemset");
  }
  for (int c = 0; c < 6; c++)
    cuda(cudaMalloc(&b[c], nb * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&R, nR * sizeof(double2)), "cudaMalloc");
  for (int q = 0; q < 2; q++)
    cuda(cudaMalloc(&T[q], g.nf * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&dP, 6 * nblock * sizeof(double)), "cudaMalloc");
  if ((P = (double *)malloc(6 * nblock * sizeof(double))) == NULL ||
      (host = (double *)malloc(nloc * sizeof(double))) == NULL)
    fail("malloc failed%s%s", "", "");
  ws = 0;
  plan(&px, nx, nx, nx, nx, nx, 1, g.kc * g.Q, CUFFT_Z2Z, &ws);
  plan(&pz, nz, nz, nz, nz, nz, 1, g.pc * g.Q, CUFFT_Z2Z, &ws);
  plan(&py[0], ny, 2 * g.Q, 2 * g.Q, g.Q, g.Q, 1, g.X1 * g.kc, CUFFT_D2Z, &ws);
  plan(&py[1], ny, g.Q, g.Q, 2 * g.Q, 2 * g.Q, 1, g.X1 * g.kc, CUFFT_Z2D, &ws);
  cuda(cudaMalloc(&work, ws > 0 ? ws : 1), "cudaMalloc");
  for (int q = 0; q < 4; q++) {
    cufft(cufftSetWorkArea(*all[q], work), "cufftSetWorkArea");
    cufft(cufftSetStream(*all[q], st), "cufftSetStream");
  }
  cuda(cudaMemGetInfo(&fr, &tot), "cudaMemGetInfo");
  char host_name[256] = "?";
  gethostname(host_name, sizeof host_name);
  fprintf(stderr,
          "boxsymmp: rank %d/%d %s gpu %d: p %ld+%ld k %ld+%ld, fields %.3f GB, R+T %.3f GB, "
          "work %.3f GB, gpu used %.3f GB\n",
          rank, size, host_name, local % ndev, g.p0, g.pc, g.k0, g.kc,
          (6 * nb + 6 * g.nf) * sizeof(double2) / 1e9, (nR + 2 * g.nf) * sizeof(double2) / 1e9,
          ws / 1e9, (tot - fr) / 1e9);
  tstep = 0;
  name = ipath;
  if (rpath != NULL) {
    const char *p = strrchr(rpath, '.');
    if (p == NULL || (tstep = strtol(p + 1, &end, 10), *end != '\0' || end == p + 1))
      fail("'%s' is not named q.STEP%s", rpath, "");
    name = rpath;
  }
  /* file [c][i][k][j]: my rows are k0..k0+kc-1 of every i: X1 contiguous pieces of kc ny doubles */
  struct stat sb;
  tio = MPI_Wtime();
  if ((fd = open(name, O_RDONLY)) < 0)
    fail("fail to open '%s'%s", name, "");
  if (fstat(fd, &sb) != 0 || sb.st_size != (off_t)(3 * g.X1 * g.Z1 * ny * sizeof(double)))
    fail("'%s' has wrong size (not a q. file of this -x -y -z?)%s", name, "");
  for (int c = 0; c < 3; c++) {
    for (long i = 0; i < g.X1; i++)
      if (rw(0, fd, host + i * g.kc * ny, g.kc * ny * sizeof(double),
             (off_t)((((c * g.X1 + i) * g.Z1 + g.k0) * ny) * sizeof(double))) != 0)
        fail("fail to read '%s'%s", name, "");
    cuda(cudaMemcpy2D(b[0], 2 * g.Q * sizeof(double), host, ny * sizeof(double),
                      ny * sizeof(double), g.X1 * g.kc, cudaMemcpyHostToDevice),
         "cudaMemcpy2D");
    analysis(1, b, ux + c, uz + c);
    cuda(cudaMemcpy(u[c], b[0], g.nf * sizeof(double2), cudaMemcpyDeviceToDevice), "cudaMemcpy");
  }
  close(fd);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    fprintf(stderr, "boxsymmp: read '%s' in %.1f s\n", name, MPI_Wtime() - tio);
  start = tstep;
  t = tstep * dt;
  devsync();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for (;;) {
    int diag = tstep % 10 == 0, dump = nd > 0 && tstep % nd == 0;
    for (long q = 0; q < nT; q++)
      if (Ts[q] == tstep)
        dump = 1;
    dump = dump && (rpath == NULL || tstep != start);
    if (diag || dump) {
      spread<<<grid(g.nf), nthread, 0, st>>>(g, u[0], u[1], u[2], b[0], b[1], b[2], b[3], b[4], b[5]);
      launched("spread");
      synthesis(6, b, ox6, oz6);
    }
    if (diag) {
      norms<<<nblock, nthread, 0, st>>>(g, u[0], u[1], u[2], dP);
      launched("norms");
      collect(6, 0, dP, P, sums);
      peaks<<<nblock, nthread, 0, st>>>(g, s, cx, cy, cz, (double *)b[0], (double *)b[1],
                                 (double *)b[2], (double *)b[3], (double *)b[4],
                                 (double *)b[5], dP);
      launched("peaks");
      collect(2, 1, dP, P, mx);
      double E = sums[0] * s * s, Om = sums[1] * s * s, eps = 2 * nu * Om, u2 = 2 * E / 3;
      double q[] = {t,
                    E,
                    Om,
                    dt * mx[0],
                    sqrt(mx[1]),
                    kc * pow(nu * nu * nu / eps, 0.25),
                    u2 * sqrt(15 / (nu * eps)),
                    sqrt(10 * nu * E / eps),
                    pi / (2 * u2) * sums[2] * s * s,
                    sums[3] / sums[0],
                    sums[4] / sums[0],
                    sums[5] / sums[0]};
      char sv[64];
      if (!(isfinite(E) && isfinite(Om) && isfinite(mx[0]) && isfinite(mx[1]))) {
        snprintf(sv, sizeof sv, "%ld", tstep);
        fail("non-finite E, Omega, C or wmax at step %s%s", sv, "");
      }
      if ((g.f ? 1.5 : 1.0) * dt * mx[0] > 3.3) {
        snprintf(sv, sizeof sv, "%ld (C = %g)", tstep, dt * mx[0]);
        fail("Courant number too large at step %s%s: reduce -s", sv, "");
      }
      if (rank == 0) {
        if (tstep == start)
          printf("step t E Omega C wmax keta Rlambda lambda L tailx taily tailz\n");
        printf("% 10ld", tstep);
        for (int i = 0; i < (int)(sizeof q / sizeof *q); i++)
          printf(" % .16e", q[i]);
        printf("\n");
        if (fflush(stdout) != 0 || ferror(stdout))
          fail("fail to write stdout%s%s", "", "");
      }
    }
    if (dump) {
      sprintf(path, "q.%08ld.tmp", tstep);
      tio = MPI_Wtime();
      bad = 0;
      if (rank == 0 &&
          ((fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644)) < 0 ||
           ftruncate(fd, (off_t)(3 * g.X1 * g.Z1 * ny * sizeof(double))) != 0 || close(fd) != 0))
        bad = 1;
      if (any(bad))
        fail("fail to open '%s'%s", path, "");
      if ((fd = open(path, O_WRONLY)) < 0)
        bad = 1;
      for (int c = 0; c < 3 && !bad; c++) {
        cuda(cudaMemcpy2D(host, ny * sizeof(double), b[c], 2 * g.Q * sizeof(double),
                          ny * sizeof(double), g.X1 * g.kc, cudaMemcpyDeviceToHost),
             "cudaMemcpy2D");
        for (long l = 0; l < nloc; l++)
          host[l] *= s;
        for (long i = 0; i < g.X1 && !bad; i++)
          if (rw(1, fd, host + i * g.kc * ny, g.kc * ny * sizeof(double),
                 (off_t)((((c * g.X1 + i) * g.Z1 + g.k0) * ny) * sizeof(double))) != 0)
            bad = 1;
      }
      if (fd >= 0 && (fsync(fd) != 0 || close(fd) != 0))
        bad = 1;
      if (any(bad))
        fail("fail to write '%s'%s", path, "");
      char final[FILENAME_MAX];
      sprintf(final, "q.%08ld", tstep);
      if (rank == 0 && rename(path, final) != 0)
        bad = 1;
      if (any(bad))
        fail("fail to rename '%s'%s", path, "");
      if (rank == 0)
        fprintf(stderr, "boxsymmp: wrote '%s' in %.1f s\n", final, MPI_Wtime() - tio);
    }
    if (ne > 0 && tstep % ne == 0) {
      double dk = fmin(g.ax, fmin(g.ay, g.az));
      long ns = (long)(sqrt(pow(g.ax * (nx / 2), 2) + pow(g.ay * (ny / 2), 2) + pow(g.az * (nz / 2), 2)) / dk) + 2;
      long nbin = ns + nx / 2 + 1 + ny / 2 + 1 + nz / 2 + 1;
      double *dE, *hE;
      cuda(cudaMalloc(&dE, nbin * sizeof(double)), "cudaMalloc");
      if ((hE = (double *)malloc(nbin * sizeof(double))) == NULL)
        fail("malloc failed%s%s", "", "");
      cuda(cudaMemsetAsync(dE, 0, nbin * sizeof(double), st), "cudaMemsetAsync");
      spectrum<<<grid(g.nf), nthread, 0, st>>>(g, dk, ns, u[0], u[1], u[2], dE);
      launched("spectrum");
      devsync();
      cuda(cudaMemcpy(hE, dE, nbin * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
      cuda(cudaFree(dE), "cudaFree");
      MPI_Allreduce(MPI_IN_PLACE, hE, nbin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (rank == 0) {
        FILE *fe;
        sprintf(path, "e.%08ld", tstep);
        if ((fe = fopen(path, "w")) == NULL)
          fail("fail to open '%s'%s", path, "");
        fprintf(fe, "dir k E\n# t = %.16e\n# shell bin width %.16e; 1D bins in kx, ky, kz\n", t, dk);
        for (long q = 0; q < nbin; q++) {
          char d = q < ns ? 's' : q < ns + nx / 2 + 1 ? 'x' : q < ns + nx / 2 + 1 + ny / 2 + 1 ? 'y' : 'z';
          long m = d == 's' ? q : d == 'x' ? q - ns : d == 'y' ? q - ns - nx / 2 - 1
                                                              : q - ns - nx / 2 - 1 - ny / 2 - 1;
          double k = d == 's' ? m * dk : d == 'x' ? m * g.ax : d == 'y' ? m * g.ay : m * g.az;
          fprintf(fe, "%c %.10g %.16e\n", d, k, hE[q] * s * s);
        }
        if (ferror(fe) || fclose(fe) != 0)
          fail("fail to write '%s'%s", path, "");
      }
      free(hE);
    }
    if (nL > 0 && tstep % nL == 0 && t >= L0 - 1e-9 && (rpath == NULL || tstep != start || L0start)) {
      /* pair centre from the omega_z^2 profiles of the full box (omega_z = d_x v - d_y u synthesized on the
         slab), then the velocity gradient (9 syntheses through b) and lambda2 only for the planes
         i <= imax the crop needs; crop gathered on rank 0 and mirrored to the full box */
      double tl = MPI_Wtime(), Lx = X * pi, Ly = Y * pi, Lz = Z * pi;
      static double *lwd;
      if (lwd == NULL)
        cuda(cudaMalloc(&lwd, (g.X1 + ny) * sizeof(double)), "cudaMalloc");
      {
        int ox2[2] = {ux[1] ^ 1, ux[0]}, oz2[2] = {uz[1], uz[0]};
        lderiv<<<grid(g.nf), nthread, 0, st>>>(g, 0, u[1], b[0], s);
        lderiv<<<grid(g.nf), nthread, 0, st>>>(g, 1, u[0], b[1], s);
        launched("lderiv");
        synthesis(2, b, ox2, oz2);
        lwz2<<<grid(nloc), nthread, 0, st>>>(g, (const double *)b[0], (const double *)b[1], (double *)b[2]);
        launched("lwz2");
      }
      double *wz2 = (double *)b[2];
      lsumx<<<g.X1, nthread, 0, st>>>(g.kc * ny, wz2, lwd);
      lsumy<<<ny, nthread, 0, st>>>(g, wz2, lwd + g.X1);
      launched("lsum");
      devsync();
      double *lw = (double *)malloc((g.X1 + ny + nx) * sizeof(double));
      if (lw == NULL)
        fail("malloc failed%s%s", "", "");
      cuda(cudaMemcpy(lw, lwd, (g.X1 + ny) * sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
      MPI_Allreduce(MPI_IN_PLACE, lw, g.X1 + ny, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      /* full x profile: x_I = (I - Nx) dx, reduced index |x| / dx */
      double *wfx = lw + g.X1 + ny;
      for (long I = 0; I < nx; I++)
        wfx[I] = lw[labs(I - g.Nx)];
      double c[2] = {lcentre(wfx, nx, Lx), lcentre(lw + g.X1, ny, Ly)};
      double ldx = Lx / nx, ldy = Ly / ny;
      long ic = ((lround((c[0] + Lx / 2) / ldx) % nx) + nx) % nx,
           jc = ((lround((c[1] + Ly / 2) / ldy) % ny) + ny) % ny;
      long hx = (long)(HH / ldx), hy = (long)(HH / ldy);
      if (hx > nx / 2 - 1)
        hx = nx / 2 - 1;
      if (hy > ny / 2 - 1)
        hy = ny / 2 - 1;
      long n0 = 2 * hx + 1, n1 = 2 * hy + 1, imax = 0;
      for (long a = 0; a < n0; a++) {
        long I = ((ic - hx + a) % nx + nx) % nx, d = I - g.Nx;
        d = d > g.Nx ? d - nx : d < -g.Nx ? d + nx : d;
        if (labs(d) > imax)
          imax = labs(d);
      }
      long nlam = (imax + 1) * g.kc * ny;
      float *LG[3][3];
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          cuda(cudaMalloc(&LG[i][j], nlam * sizeof(float)), "cudaMalloc");
      for (int h = 0; h < 2; h++) {
        int m = h == 0 ? 6 : 3, ox9[6], oz9[6];
        for (int q = 0; q < m; q++) {
          int ij = 6 * h + q, i = ij / 3, j = ij % 3;
          lderiv<<<grid(g.nf), nthread, 0, st>>>(g, j, u[i], b[q], s);
          launched("lderiv");
          ox9[q] = ux[i] ^ (j == 0);
          oz9[q] = uz[i] ^ (j == 2);
        }
        synthesis(m, b, ox9, oz9);
        for (int q = 0; q < m; q++) {
          int ij = 6 * h + q;
          ltofloat<<<grid(nlam), nthread, 0, st>>>((imax + 1) * g.kc, ny, g.Q, (double *)b[q], LG[ij / 3][ij % 3]);
          launched("ltofloat");
        }
      }
      Grad gg;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          gg.g[i][j] = LG[i][j];
      float *lam = (float *)b[0];
      llam2<<<grid(nlam), nthread, 0, st>>>(nlam, gg, lam);
      launched("llam2");
      devsync();
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          cuda(cudaFree(LG[i][j]), "cudaFree");
      /* my block [i <= imax][kloc][jwin] */
      long nblk = (imax + 1) * g.kc * n1;
      float *hl = (float *)malloc(nlam * sizeof(float)), *blk = (float *)malloc(nblk * sizeof(float));
      if (hl == NULL || blk == NULL)
        fail("malloc failed%s%s", "", "");
      cuda(cudaMemcpy(hl, lam, nlam * sizeof(float), cudaMemcpyDeviceToHost), "cudaMemcpy");
      for (long i = 0; i <= imax; i++)
        for (long k = 0; k < g.kc; k++)
          for (long jw = 0; jw < n1; jw++) {
            long yg = ((jc - hy + jw) % ny + ny) % ny;
            blk[(i * g.kc + k) * n1 + jw] = hl[(i * g.kc + k) * ny + yg];
          }
      free(hl);
      int *cnt = NULL, *dsp = NULL;
      float *all = NULL;
      if (rank == 0) {
        cnt = (int *)malloc(2 * size * sizeof(int));
        dsp = cnt + size;
        long o = 0;
        for (int q = 0; q < size; q++) {
          cnt[q] = (imax + 1) * (first(g.Z1, size, q + 1) - first(g.Z1, size, q)) * n1;
          dsp[q] = o;
          o += cnt[q];
        }
        if ((all = (float *)malloc(o * sizeof(float))) == NULL)
          fail("malloc failed%s%s", "", "");
      }
      MPI_Gatherv(blk, nblk, MPI_FLOAT, all, cnt, dsp, MPI_FLOAT, 0, MPI_COMM_WORLD);
      free(blk);
      char lp[FILENAME_MAX], lt[FILENAME_MAX], shape[64], hd[128];
      sprintf(lp, "l.%08ld.npy", tstep);
      sprintf(lt, "l.%08ld.npy.tmp", tstep);
      bad = 0;
      if (rank == 0) {
        snprintf(shape, sizeof shape, "(%ld, %ld, %ld)", n0, n1, nz);
        int hn = npyhead(hd, "<f4", shape);
        double xyz[6] = {-Lx / 2 + ldx * (ic - hx), -Ly / 2 + ldy * (jc - hy), 0, ldx, ldy, Lz / nz};
        char sp[FILENAME_MAX];
        sprintf(sp, "l.%08ld_xyz.npy", tstep);
        if (npysmall(sp, xyz, 6) != 0)
          bad = 1;
        sprintf(sp, "l.%08ld_ctr.npy", tstep);
        if (npysmall(sp, c, 2) != 0)
          bad = 1;
        float *row = (float *)malloc(n1 * nz * sizeof(float));
        FILE *fl = fopen(lt, "w");
        if (row == NULL || fl == NULL || fwrite(hd, 1, hn, fl) != (size_t)hn)
          bad = 1;
        for (long a = 0; a < n0 && !bad; a++) {
          long I = ((ic - hx + a) % nx + nx) % nx, d = I - g.Nx;
          d = d > g.Nx ? d - nx : d < -g.Nx ? d + nx : d;
          long i = labs(d);
          for (long K = 0; K < nz; K++) {
            long k = K <= g.Nz ? K : nz - K, q = owner(g.Z1, size, k), k0 = first(g.Z1, size, q),
                 kq = first(g.Z1, size, q + 1) - k0;
            const float *src = all + dsp[q] + (i * kq + k - k0) * n1;
            for (long jw = 0; jw < n1; jw++)
              row[jw * nz + K] = src[jw];
          }
          if (fwrite(row, sizeof(float), n1 * nz, fl) != (size_t)(n1 * nz))
            bad = 1;
        }
        if (fl == NULL || ferror(fl) || fclose(fl) != 0 || bad || rename(lt, lp) != 0)
          bad = 1;
        free(row);
        free(all);
        free(cnt);
        fprintf(stderr, "boxsymmp: lambda2 '%s' (%ld x %ld x %ld, ctr %.6g %.6g) in %.1f s\n", lp, n0, n1,
                nz, c[0], c[1], MPI_Wtime() - tl);
      }
      free(lw);
      if (any(bad))
        fail("fail to write '%s'%s", lt, "");
    }
    if (tstep >= nt)
      break;
    for (int r = 0; r < 5; r++) {
      spread<<<grid(g.nf), nthread, 0, st>>>(g, u[0], u[1], u[2], b[0], b[1], b[2], b[3], b[4], b[5]);
      launched("spread");
      synthesis(6, b, ox6, oz6);
      cross<<<grid(g.nr), nthread, 0, st>>>(g, s, (double *)b[0], (double *)b[1], (double *)b[2],
                                     (double *)b[3], (double *)b[4], (double *)b[5]);
      launched("cross");
      analysis(3, b, ux, uz);
      update<<<grid(g.nf), nthread, 0, st>>>(g, r == 4, nu, dt, A[r], B[r], b[0], b[1], b[2], u[0],
                                      u[1], u[2], du[0], du[1], du[2]);
      launched("update");
    }
    tstep++;
    t = tstep * dt;
  }
  devsync();
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if (rank == 0)
    fprintf(stderr, "boxsymmp: %ld steps on %d ranks in %.2f s: %.4f s/step\n", nt - start, size,
            t1 - t0, nt > start ? (t1 - t0) / (nt - start) : 0.0);
  if (!usempi)
    ncclCommDestroy(nccl);
  MPI_Finalize();
}
