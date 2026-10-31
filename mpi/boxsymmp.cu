/* boxsymmp: MPI multi-GPU version of gpu/boxsym.cu (box restricted to the Sx, Sz symmetric fields of the
   colliding tubes: u(-x,y,z) = (-u,v,w)(x,y,z), u(x,y,-z) = (u,v,-w)(x,y,z)). Same options and out lines
   as boxsym; reads and writes only q.STEP (double [3][nx/2+1][nz/2+1][ny]; gpu/symq makes one from a full
   box u.STEP, gpu/symfull expands one back).
   Decomposition (X1 = nx/2+1, Z1 = nz/2+1, Q = ny/2+1, P ranks, one GPU each): physical fields and the
   x-stage in z-slabs [i][kloc][j], spectral fields and the z-stage in x-slabs [ploc][r][q]; one
   all-to-all per 3-D transform between the x- and z-stages. y is R2C/C2R, x and z are full-length
   complex FFTs of parity-extended lines (as in boxsym, tiled transposes). The all-to-all is NCCL grouped
   send/recv on its own stream, pipelined against the compute of the neighbouring fields (MPI only
   bootstraps NCCL and reduces the diagnostics).
   yao-spec additions: -I yao (Yao & Hussain 2020a tube pair generated on the GPU, no IC file), Gz column
   (circulation of omega_z through z = 0, x < 0, by Stokes from v on x = 0 and x = -pi), -c C (adaptive
   time step: dt = C / Courant rate at every diagnostic step, -s is the largest dt; -d -e -L -T are then times
   and q.STEP.t holds t), -S (all-to-all in float: x-, z-stage results rounded to fp32 before the exchange). */
#include <cuda_fp16.h>
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
  const double *mfx, *mfy, *mfz; /* mask (filter or 2/3 rule) = mfx[p] mfy[q] mfz[r] */
};
static G g;
static int rank, size;
static cufftHandle px, pz, py[2];
static double2 *R, *T[2];
static int *scz, *sdz, *rcz, *rdz; /* alltoallv: z-slab sections <-> x-slab chunks */
static ncclComm_t nccl;
static int usempi; /* BOXSYM_MPI=1: CUDA-aware MPI_Alltoallv instead of NCCL (e.g. several ranks per GPU) */
static MPI_Datatype cplx, cplxf;
static int lowp; /* -S: the all-to-all buffers hold float2 (set only around the RK stages) */
static int syn6; /* BOXSYM_SYN6=1: RK stages synthesize u and omega by six exchanges (the old path) */
static int nocomm; /* BOXSYM_NOCOMM=1 (timing only): the exchanges are skipped, results are garbage */
static int nofuse; /* BOXSYM_NOFUSE=1: diagnostics from a synthesis of their own (the old way) */
static int nl; /* nl_a/nl_b stage (needs b[6]) */
/* BOXSYM_HC=nch > 1: between the RK stages of a step the last z-stage, the update, spread5 and the first exchange of
   the next stage run in nch chunks of p rows (plans pzc, exchange counts per chunk); headdone: field 0 of this stage
   was sent that way */
static int nch = 1, headdone, nrb = 1;
static cufftHandle pxb[4], pzb[4];
static cufftHandle pzc[4];
static int *cnt_s[4], *dsp_s[4], *cnt_r[4], *dsp_r[4];
static int selfx = 1; /* the folds write the own part of an exchange directly (BOXSYM_SELFCOPY=1: copy instead) */
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
static void a2a(const void *sb_, const int *nsend, const int *sd, void *rb_, const int *rc,
                const int *rd) {
  size_t es = lowp ? sizeof(float2) : sizeof(double2);
  const char *sb = (const char *)sb_;
  char *rb = (char *)rb_;
  ncclDataType_t ty = lowp ? ncclFloat : ncclDouble;
  if (usempi) {
    devsync();
    MPI_Alltoallv(sb, nsend, sd, lowp ? cplxf : cplx, rb, rc, rd, lowp ? cplxf : cplx, MPI_COMM_WORLD);
    return;
  }
  if (nocomm)
    return;
  if (rc[rank] > 0 && !selfx)
    cuda(cudaMemcpyAsync(rb + rd[rank] * es, sb + sd[rank] * es, rc[rank] * es,
                         cudaMemcpyDeviceToDevice, sc),
         "cudaMemcpyAsync");
  nccheck(ncclGroupStart(), "ncclGroupStart");
  for (int q = 0; q < size; q++) {
    if (q == rank)
      continue;
    if (nsend[q] > 0)
      nccheck(ncclSend(sb + sd[q] * es, 2 * (size_t)nsend[q], ty, q, nccl, sc), "ncclSend");
    if (rc[q] > 0)
      nccheck(ncclRecv(rb + rd[q] * es, 2 * (size_t)rc[q], ty, q, nccl, sc), "ncclRecv");
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
__device__ static double2 ld2(const double2 *a, long i) { return a[i]; }
__device__ static double2 ld2(const float2 *a, long i) {
  float2 v = a[i];
  return make_double2(v.x, v.y);
}
__device__ static void st2(double2 *a, long i, double2 v) { a[i] = v; }
__device__ static void st2(float2 *a, long i, double2 v) { a[i] = make_float2((float)v.x, (float)v.y); }
/* x-stage kernels on the lines bb .. bb + Mb - 1 of M (R holds only those: R[(b - bb) n + m]) */
template <class In> __global__ static void extx(long n, int odd, long M, const In *a, double2 *R, long bb, long Mb) {
  __shared__ double2 t[TL][TL + 1];
  long b0 = bb + blockIdx.x * (long)TL, m0 = blockIdx.y * (long)TL, N = n / 2;
  int x = threadIdx.x, y = threadIdx.y;
  for (int j = y; j < TL; j += TR) {
    long m = m0 + j, b = b0 + x, h = m <= N ? m : n - m;
    if (m < n && b < bb + Mb) {
      double2 v = odd && (h == 0 || h == N) ? make_double2(0, 0) : ld2(a, h * M + b);
      t[j][x] = odd && m > N ? make_double2(-v.x, -v.y) : v;
    }
  }
  __syncthreads();
  for (int j = y; j < TL; j += TR) {
    long b = b0 + j, m = m0 + x;
    if (m < n && b < bb + Mb)
      R[(b - bb) * n + m] = t[x][j];
  }
}
/* fold* write the own rank's part straight into the receive buffer (self, rows sp0 .. sp0 + spc - 1 of the
   fold output for foldx, section srank for foldz), so the exchange skips its local copy */
template <class Out>
__global__ static void foldx(long n, long M, const double2 *R, Out *a, Out *self, long sp0, long spc, long bb,
                             long Mb) {
  __shared__ double2 t[TL][TL + 1];
  long b0 = bb + blockIdx.x * (long)TL, p0 = blockIdx.y * (long)TL, N1 = n / 2 + 1;
  int x = threadIdx.x, y = threadIdx.y;
  for (int j = y; j < TL; j += TR) {
    long b = b0 + j, p = p0 + x;
    if (p < N1 && b < bb + Mb)
      t[j][x] = R[(b - bb) * n + p];
  }
  __syncthreads();
  for (int j = y; j < TL; j += TR) {
    long p = p0 + j, b = b0 + x;
    if (p < N1 && b < bb + Mb) {
      if (self != NULL && p >= sp0 && p < sp0 + spc)
        st2(self, (p - sp0) * M + b, t[x][j]);
      else
        st2(a, p * M + b, t[x][j]);
    }
  }
}
template <class In>
__global__ static void extz(long n, int odd, long Q, long Z1, long pc, long S, const In *a, double2 *R, long pl0) {
  __shared__ double2 t[TL][TL + 1];
  __shared__ long ro[TL]; /* zrow of the tile's rows: the 64-bit divisions once per row, not per element */
  long q0 = blockIdx.x * (long)TL, m0 = blockIdx.y * (long)TL, pl = pl0 + blockIdx.z, N = n / 2;
  int x = threadIdx.x, y = threadIdx.y;
  if (y == 0) {
    long m = m0 + x, h = m <= N ? m : n - m;
    ro[x] = m < n ? zrow(Z1, S, pc, Q, pl, h) : 0;
  }
  __syncthreads();
  for (int j = y; j < TL; j += TR) {
    long m = m0 + j, q = q0 + x, h = m <= N ? m : n - m;
    if (m < n && q < Q) {
      double2 v = odd && (h == 0 || h == N) ? make_double2(0, 0) : ld2(a, ro[j] + q);
      t[j][x] = odd && m > N ? make_double2(-v.x, -v.y) : v;
    }
  }
  __syncthreads();
  for (int j = y; j < TL; j += TR) {
    long q = q0 + j, m = m0 + x;
    if (m < n && q < Q)
      R[((pl - pl0) * Q + q) * n + m] = t[x][j];
  }
}
template <class Out>
__global__ static void foldz(long n, long Q, long Z1, long pc, long S, const double2 *R, Out *a, Out *self,
                             long srank, long pl0) {
  __shared__ double2 t[TL][TL + 1];
  __shared__ long ro[TL];
  long q0 = blockIdx.x * (long)TL, k0 = blockIdx.y * (long)TL, pl = pl0 + blockIdx.z;
  int x = threadIdx.x, y = threadIdx.y;
  for (int j = y; j < TL; j += TR) {
    long q = q0 + j, k = k0 + x;
    if (k < Z1 && q < Q)
      t[j][x] = R[((pl - pl0) * Q + q) * n + k];
  }
  if (y == 0) {
    /* row offsets (64-bit divisions once per row); own section (self) encoded as -1 - offset */
    long k = k0 + x, sct = ((k + 1) * S - 1) / Z1, r0 = Z1 * sct / S, cnt = Z1 * (sct + 1) / S - r0;
    ro[x] = self != NULL && sct == srank ? -1 - (pl * cnt + k - r0) * Q : pc * r0 * Q + (pl * cnt + k - r0) * Q;
  }
  __syncthreads();
  for (int j = y; j < TL; j += TR) {
    long k = k0 + j, q = q0 + x;
    if (k < Z1 && q < Q) {
      long o = ro[j];
      if (o < 0)
        st2(self, -1 - o + q, t[x][j]);
      else
        st2(a, o + q, t[x][j]);
    }
  }
}
static dim3 tiles(long a, long b, long c) { return dim3((a + TL - 1) / TL, (b + TL - 1) / TL, c); }
static void rec(int k, int c, cudaStream_t s) { cuda(cudaEventRecord(ev[k][c], s), "cudaEventRecord"); }
static void waitev(cudaStream_t s, int k, int c) {
  cuda(cudaStreamWaitEvent(s, ev[k][c], 0), "cudaStreamWaitEvent");
}
/* exchange buffers: a ring of 2 (double, T[0], T[1]) or, with -S, of 4 (float halves of T[0], T[1]) */
static int nring(void) { return lowp ? 4 : 2; }
static void *tb(int s) {
  return lowp ? (void *)((float2 *)T[s % 4 / 2] + s % 2 * g.nf) : (void *)T[s % 2];
}
/* z-stage of field a (z parity oz) and its all-to-all; s = sequence number of the exchange in this synthesis
   (events [k][s], send buffer tb(s)). Events: 0 compute side done with tb(s), 1 all-to-all s done. */
static void waitrecv(cudaStream_t q, int s);
/* chunk j of the p rows (local rows pa .. pb - 1; j = -1: all) */
static long cha(int j) { return j < 0 ? 0 : first(g.pc, nch, j); }
static long chb(int j) { return j < 0 ? g.pc : first(g.pc, nch, j + 1); }
/* BOXSYM_RB = nrb > 1: x and z stages in nrb batches of lines, R holds one batch (plans pxb, pzb) */
static long xba(int k) { return first(g.kc * g.Q, nrb, k); }
static long zba(int k) { return first(g.pc, nrb, k); }
/* src: the spectral field (x-slab), a: its receive buffer (src == a allowed: the exchange waits for the z-stage) */
static void zsendc(int s, const double2 *src, double2 *a, int oz, int j) {
  int nk = j < 0 ? nrb : 1;
  for (int k = 0; k < nk; k++) {
    long pa = j < 0 ? zba(k) : cha(j), pb = j < 0 ? zba(k + 1) : chb(j);
    cufftHandle p = j < 0 ? pzb[k] : pzc[j];
    extz<<<tiles(g.Q, g.nz, pb - pa), dim3(TL, TR), 0, st>>>(g.nz, oz, g.Q, g.Z1, g.pc, 1, src, R, pa);
    launched("extz");
    cufft(cufftExecZ2Z(p, R, R, CUFFT_INVERSE), "cufftExecZ2Z");
    if (j <= 0 && k == 0 && s >= nring())
      waitrecv(st, s - nring());
    if (lowp)
      foldz<<<tiles(g.Q, g.Z1, pb - pa), dim3(TL, TR), 0, st>>>(g.nz, g.Q, g.Z1, g.pc, size, R, (float2 *)tb(s),
                                                                 selfx ? (float2 *)a + rdz[rank] : NULL, rank, pa);
    else
      foldz<<<tiles(g.Q, g.Z1, pb - pa), dim3(TL, TR), 0, st>>>(g.nz, g.Q, g.Z1, g.pc, size, R, (double2 *)tb(s),
                                                                 selfx ? a + rdz[rank] : NULL, rank, pa);
    launched("foldz");
  }
  if (j < 0) {
    rec(0, s, st);
    waitev(sc, 0, s);
    a2a(tb(s), scz, sdz, a, rcz, rdz);
    rec(1, s, sc);
  } else {
    /* chunk exchange: events [3][j] (compute done), [3][4 + j] (exchange done) */
    rec(3, j, st);
    waitev(sc, 3, j);
    a2a(tb(s), cnt_s[j], dsp_s[j], a, cnt_r[j], dsp_r[j]);
    rec(3, 4 + j, sc);
  }
}
static void zsend(int s, double2 *a, int oz) { zsendc(s, a, a, oz, -1); }
/* field s = 0 may have been sent in chunks at the end of the previous stage (headdone) */
static void waitrecv(cudaStream_t q, int s) {
  if (s == 0 && headdone)
    for (int j = 0; j < nch; j++)
      waitev(q, 3, 4 + j);
  else
    waitev(q, 1, s);
}
/* x-stage and y C2R of a received field (x parity ox) */
static void xstage(double2 *a, int ox) {
  for (int k = 0; k < nrb; k++) {
    long bb = xba(k), Mb = xba(k + 1) - bb;
    if (lowp)
      extx<<<tiles(Mb, g.nx, 1), dim3(TL, TR), 0, st>>>(g.nx, ox, g.kc * g.Q, (const float2 *)a, R, bb, Mb);
    else
      extx<<<tiles(Mb, g.nx, 1), dim3(TL, TR), 0, st>>>(g.nx, ox, g.kc * g.Q, (const double2 *)a, R, bb, Mb);
    launched("extx");
    cufft(cufftExecZ2Z(pxb[k], R, R, CUFFT_INVERSE), "cufftExecZ2Z");
    foldx<<<tiles(Mb, g.X1, 1), dim3(TL, TR), 0, st>>>(g.nx, g.kc * g.Q, R, a, (double2 *)NULL, 0, 0, bb, Mb);
    launched("foldx");
  }
  cufft(cufftExecZ2D(py[1], a, (double *)a), "cufftExecZ2D");
}
/* n fields a[c] with parities ox[c], oz[c]. Pipelined: the all-to-all of field c (stream sc, buffer
   T[c % 2]) overlaps the compute of fields c - 1 and c + 1 (stream st). */
static void synthesis(int n, double2 **a, const int *ox, const int *oz) {
  for (int c = 0; c <= n; c++) {
    if (c < n)
      zsend(c, a[c], oz[c]);
    if (c >= 1) {
      waitev(st, 1, c - 1);
      xstage(a[c - 1], ox[c - 1]);
    }
  }
}
/* omega_z = i kx v - i ky u in the mixed layout after the exchange ([p][kloc][q]: x, y spectral, z physical) */
template <class E> __global__ static void vortm(G g, const E *a, const E *b, E *c) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= g.X1 * g.kc * g.Q)
    return;
  double kx = g.ax * (l / (g.kc * g.Q)), ky = g.ay * (l % g.Q);
  double2 A = ld2(a, l), B = ld2(b, l);
  st2(c, l, make_double2(-(kx * B.y - ky * A.y), kx * B.x - ky * A.x));
}
static void vortz(double2 *a, double2 *b, double2 *c) {
  if (lowp)
    vortm<<<grid(g.X1 * g.kc * g.Q), nthread, 0, st>>>(g, (const float2 *)a, (const float2 *)b, (float2 *)c);
  else
    vortm<<<grid(g.X1 * g.kc * g.Q), nthread, 0, st>>>(g, (const double2 *)a, (const double2 *)b, c);
  launched("vortm");
}
/* u, omega on the physical slab from five exchanges: b0..b4 = u, v, w, omega_x, omega_y (spectral, spread5);
   omega_z = d_x v - d_y u is formed in the mixed layout once u and v have arrived (b5). The last field to
   arrive needs only its own x-stage. */
static void synth5(double2 **b) {
  for (int c = 0; c < 4; c++)
    zsend(c, b[c], oz6[c]);
  waitev(st, 1, 0);
  waitev(st, 1, 1);
  vortz(b[0], b[1], b[5]);
  xstage(b[0], ox6[0]);
  xstage(b[1], ox6[1]);
  zsend(4, b[4], oz6[4]);
  xstage(b[5], ox6[5]);
  for (int c = 2; c < 5; c++) {
    waitev(st, 1, c);
    xstage(b[c], ox6[c]);
  }
}
/* analysis pieces. asend: y R2C, x-stage of field a (x parity ox), its exchange into tb(s) (s = sequence number,
   events [k][s]); azstage: z-stage of the received field s into a (z parity oz). Event 2: tb(s) read, free again. */
static void asend(int s, double2 *a, int ox) {
  cufft(cufftExecD2Z(py[0], (double *)a, a), "cufftExecD2Z");
  for (int k = 0; k < nrb; k++) {
    long bb = xba(k), Mb = xba(k + 1) - bb;
    extx<<<tiles(Mb, g.nx, 1), dim3(TL, TR), 0, st>>>(g.nx, ox, g.kc * g.Q, (const double2 *)a, R, bb, Mb);
    launched("extx");
    cufft(cufftExecZ2Z(pxb[k], R, R, CUFFT_FORWARD), "cufftExecZ2Z");
    if (lowp)
      foldx<<<tiles(Mb, g.X1, 1), dim3(TL, TR), 0, st>>>(g.nx, g.kc * g.Q, R, (float2 *)a,
                                                         selfx ? (float2 *)tb(s) + sdz[rank] : NULL, g.p0, g.pc, bb, Mb);
    else
      foldx<<<tiles(Mb, g.X1, 1), dim3(TL, TR), 0, st>>>(g.nx, g.kc * g.Q, R, a,
                                                         selfx ? (double2 *)tb(s) + sdz[rank] : NULL, g.p0, g.pc, bb, Mb);
    launched("foldx");
  }
  rec(0, s, st);
  waitev(sc, 0, s);
  a2a(a, rcz, rdz, tb(s), scz, sdz);
  rec(1, s, sc);
}
static void azstagec(int s, double2 *a, int oz, int j) {
  int nk = j < 0 ? nrb : 1;
  if (j <= 0)
    waitev(st, 1, s);
  for (int k = 0; k < nk; k++) {
    long pa = j < 0 ? zba(k) : cha(j), pb = j < 0 ? zba(k + 1) : chb(j);
    cufftHandle p = j < 0 ? pzb[k] : pzc[j];
    if (lowp)
      extz<<<tiles(g.Q, g.nz, pb - pa), dim3(TL, TR), 0, st>>>(g.nz, oz, g.Q, g.Z1, g.pc, size, (const float2 *)tb(s), R, pa);
    else
      extz<<<tiles(g.Q, g.nz, pb - pa), dim3(TL, TR), 0, st>>>(g.nz, oz, g.Q, g.Z1, g.pc, size, (const double2 *)tb(s), R, pa);
    launched("extz");
    if ((j < 0 && k == nk - 1) || j == nch - 1)
      rec(2, s, st);
    cufft(cufftExecZ2Z(p, R, R, CUFFT_FORWARD), "cufftExecZ2Z");
    foldz<<<tiles(g.Q, g.Z1, pb - pa), dim3(TL, TR), 0, st>>>(g.nz, g.Q, g.Z1, g.pc, 1, R, a, (double2 *)NULL, 0, pa);
    launched("foldz");
  }
}
static void azstage(int s, double2 *a, int oz) { azstagec(s, a, oz, -1); }
static void analysis(int n, double2 **a, const int *ox, const int *oz) {
  for (int c = 0; c <= n; c++) {
    if (c < n) {
      if (c >= nring())
        waitev(sc, 2, c - nring());
      asend(c, a[c], ox[c]);
    }
    if (c >= 1)
      azstage(c - 1, a[c - 1], oz[c - 1]);
  }
}
/* N_x = s^2 (v omega_z - w omega_y) into n0 (physical slab) */
__global__ static void cross0(G g, double s, const double *v, const double *w, const double *oy, const double *oz,
                              double *n0) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= g.nr || l % (2 * g.Q) >= g.ny)
    return;
  n0[l] = s * s * (v[l] * oz[l] - w[l] * oy[l]);
}
/* N_y = s^2 (w omega_x - u omega_z) into b0, N_z = s^2 (u omega_y - v omega_x) into b1 */
__global__ static void cross12(G g, double s, double *b0, double *b1, const double *b2, const double *b3,
                               const double *b4, const double *b5) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= g.nr || l % (2 * g.Q) >= g.ny)
    return;
  double u = b0[l], v = b1[l], w = b2[l], p = b3[l], q = b4[l], r = b5[l];
  b0[l] = s * s * (w * p - u * r);
  b1[l] = s * s * (u * q - v * p);
}
/* the nonlinear term of an RK stage with N_x started before omega_x has arrived. Roles U V W OX OY OZ N0 live in
   b[bm[role]] (the map rotates at chunked stage boundaries). u, v, w are sent from the state u[] (receive into U, V,
   W), omega_y and omega_x from OY, OX (spectral, spreadw); omega_z is formed after u, v arrive. N_x into N0 and
   its exchange (seq 5) as soon as v, w, omega_y, omega_z are physical; nl_a ends with the six fields physical
   (diagnostics may run there); nl_b leaves N spectral in N0, U, V (N_x, N_y, N_z). */
enum { U_, V_, W_, OX_, OY_, OZ_, N0_ };
static int bm[7] = {0, 1, 2, 3, 4, 5, 6};
static void nl_a(double2 **b, double2 **u, double s) {
#define B_(r) b[bm[r]]
  for (int c = headdone; c < 3; c++)
    zsendc(c, u[c], B_(c), oz6[c], -1);
  zsend(3, B_(OY_), oz6[4]);
  waitrecv(st, 0);
  waitev(st, 1, 1);
  vortz(B_(U_), B_(V_), B_(OZ_));
  xstage(B_(U_), ox6[0]);
  xstage(B_(V_), ox6[1]);
  zsend(4, B_(OX_), oz6[3]);
  xstage(B_(OZ_), ox6[5]);
  waitev(st, 1, 2);
  xstage(B_(W_), ox6[2]);
  waitev(st, 1, 3);
  xstage(B_(OY_), ox6[4]);
  cross0<<<grid(g.nr), nthread, 0, st>>>(g, s, (double *)B_(V_), (double *)B_(W_), (double *)B_(OY_),
                                         (double *)B_(OZ_), (double *)B_(N0_));
  launched("cross0");
  asend(5, B_(N0_), ux[0]);
  waitev(st, 1, 4);
  xstage(B_(OX_), ox6[3]);
}
static void nl_b(double2 **b, double s, int head) {
  cross12<<<grid(g.nr), nthread, 0, st>>>(g, s, (double *)B_(U_), (double *)B_(V_), (double *)B_(W_),
                                          (double *)B_(OX_), (double *)B_(OY_), (double *)B_(OZ_));
  launched("cross12");
  asend(6, B_(U_), ux[1]);
  azstage(5, B_(N0_), uz[0]);
  if (7 - 5 >= nring())
    waitev(sc, 2, 5);
  asend(7, B_(V_), ux[2]);
  azstage(6, B_(U_), uz[1]);
  if (!head)
    azstage(7, B_(V_), uz[2]);
}
/* omega_x, omega_y spectral (modes l0 .. l1 - 1) for nl_a */
__global__ static void spreadw(G g, const double2 *u, const double2 *v, const double2 *w, double2 *ox, double2 *oy,
                               long l0, long l1) {
  long l = l0 + blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= l1)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double2 a = u[l], b = v[l], c = w[l];
  ox[l] = make_double2(-(ky * c.y - kz * b.y), ky * c.x - kz * b.x);
  oy[l] = make_double2(-(kz * a.y - kx * c.y), kz * a.x - kx * c.x);
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
/* spectral inputs of synth5: u, v, w, omega_x, omega_y (as spread) */
__global__ static void spread5(G g, const double2 *u, const double2 *v, const double2 *w, double2 *b0, double2 *b1,
                               double2 *b2, double2 *b3, double2 *b4, long l0, long l1) {
  long l = l0 + blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= l1)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double2 a = u[l], b = v[l], c = w[l];
  b0[l] = a;
  b1[l] = b;
  b2[l] = c;
  b3[l] = make_double2(-(ky * c.y - kz * b.y), ky * c.x - kz * b.x);
  b4[l] = make_double2(-(kz * a.y - kx * c.y), kz * a.x - kx * c.x);
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
                              double2 *dw, long l0, long l1) {
  long l = l0 + blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= l1)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double kk = kx * kx + ky * ky + kz * kz, m = g.mfx[i] * g.mfy[j] * g.mfz[k], visc = nu * kk;
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
/* Yao & Hussain 2020a section 2 (as cases2/yao.py): omega on the physical slab [i][kloc][j] (rows 2Q), x = i dx,
   y = -Ly/2 + j dy, z = (k0 + kloc) dz (z >= Lz/2 wrapped to z - Lz). Tube 1 x_c = (-x0 + p cos s cos al,
   p cos s sin al, s), tube 2 its mirror in x = 0 with the opposite sign; w = w0 [1 - g(r/rc)], r = distance to
   the nearest axis point (Newton on (x - x_c(s)) . t(s) = 0 from s = z), direction t / |t| there. */
__global__ static void yaow(G g, double Lx, double Ly, double Lz, double *w0_, double *w1_, double *w2_) {
  const double x0 = 0.81, pp = 0.2, rc = 2.0 / 3, w0 = 20, K0 = 0.5 * exp(2.0) * log(2.0);
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= g.X1 * g.kc * g.ny)
    return;
  long row = l / g.ny, j = l % g.ny, i = row / g.kc, k = g.k0 + row % g.kc, r = row * 2 * g.Q + j;
  double x = i * Lx / g.nx, y = -Ly / 2 + j * Ly / g.ny, z = k * Lz / g.nz;
  if (z >= Lz / 2)
    z -= Lz;
  double o[3] = {0, 0, 0};
  for (int tb = 0; tb < 2; tb++) {
    double xc0 = tb == 0 ? -x0 : x0, al = tb == 0 ? pi / 3 : pi - pi / 3, sg = tb == 0 ? 1 : -1;
    double ca = pp * cos(al), sa = pp * sin(al), s = z, cx, cy, dx, dy, dz, tx, ty;
    for (int it = 0; it < 50; it++) {
      cx = xc0 + ca * cos(s);
      cy = sa * cos(s);
      dx = x - cx;
      dy = y - cy;
      dz = z - s;
      dz -= Lz * nearbyint(dz / Lz);
      tx = -ca * sin(s);
      ty = -sa * sin(s);
      double gg = dx * tx + dy * ty + dz, dg = -(tx * tx + ty * ty + 1) + dx * (-ca * cos(s)) + dy * (-sa * cos(s));
      double h = gg / dg;
      s -= h;
      if (fabs(h) < 1e-14)
        break;
    }
    cx = xc0 + ca * cos(s);
    cy = sa * cos(s);
    dz = z - s;
    dz -= Lz * nearbyint(dz / Lz);
    double rr = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + dz * dz), q = rr / rc, f;
    tx = -ca * sin(s);
    ty = -sa * sin(s);
    if (q >= 1)
      continue;
    f = q > 0 ? w0 * (1 - exp(-K0 / q * exp(1 / (q - 1)))) : w0;
    f *= sg / sqrt(tx * tx + ty * ty + 1);
    o[0] += f * tx;
    o[1] += f * ty;
    o[2] += f;
  }
  w0_[r] = o[0];
  w1_[r] = o[1];
  w2_[r] = o[2];
}
/* u = i k x omega / |k|^2 (omega spectral in a, b, c; the Nyquist planes and the mean set to zero) */
__global__ static void biot(G g, const double2 *a, const double2 *b, const double2 *c, double2 *u,
                            double2 *v, double2 *w) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x, i, j, k;
  double kx, ky, kz;
  if (l >= g.nf)
    return;
  mode(g, l, &i, &j, &k, &kx, &ky, &kz);
  double kk = kx * kx + ky * ky + kz * kz;
  double m = kk > 0 && i != g.Nx && 2 * j != g.ny && k != g.Nz ? 1 / kk : 0;
  double2 A = a[l], B = b[l], C = c[l];
  u[l] = make_double2(-(ky * C.y - kz * B.y) * m, (ky * C.x - kz * B.x) * m);
  v[l] = make_double2(-(kz * A.y - kx * C.y) * m, (kz * A.x - kx * C.x) * m);
  w[l] = make_double2(-(kx * B.y - ky * A.y) * m, (kx * B.x - ky * A.x) * m);
}
/* Gamma_z = int_{z=0, x<0} omega_z = int v(0, y, 0) dy - int v(-pi, y, 0) dy (Stokes; v even in x): rows
   i = 0 and i = Nx of plane k = 0 (rank with k0 = 0), one block */
__global__ static void gzsum(G g, const double *v, double *out) {
  __shared__ double sh[nthread];
  double acc = 0;
  for (long j = threadIdx.x; j < g.ny; j += blockDim.x)
    acc += v[j] - v[g.Nx * g.kc * 2 * g.Q + j];
  sh[threadIdx.x] = acc;
  __syncthreads();
  for (int q = blockDim.x / 2; q > 0; q /= 2) {
    if (threadIdx.x < q)
      sh[threadIdx.x] += sh[threadIdx.x + q];
    __syncthreads();
  }
  if (threadIdx.x == 0)
    out[0] = sh[0];
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
/* lambda2 of the velocity gradient g_[i][j] = d_j u_i (as boxmp/lam2) */
__device__ static double lam2of(const double g_[3][3]) {
  double s[3][3], o[3][3], M[3][3];
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
  return fmax(fmin(e0, e1), fmin(fmax(e0, e1), e2));
}
__global__ static void llam2(long n, Grad G_, float *lam) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n)
    return;
  double g_[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      g_[i][j] = G_.g[i][j][l];
  lam[l] = (float)lam2of(g_);
}
/* -V: |w| and lambda2 / ls for all planes from six gradient components in float (g6: d_j u_0, d_j u_1) and three in
   double physical slabs (d_j u_2, rows 2Q); writes the crop halves (planes i < ni, order [k][i][j]) and full float fields */
struct Grad6 {
  const float *g[6];
  const double *h[3];
};
__global__ static void wsnap2(long X1, long kc, long ny, long Q, long ni, Grad6 G_, double ls, __half *wa, __half *wz,
                              __half *l2, float *fa, float *fl) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= X1 * kc * ny)
    return;
  long j = l % ny, k = l / ny % kc, i = l / ny / kc, r = (l / ny) * 2 * Q + j;
  double g_[3][3];
  for (int a = 0; a < 6; a++)
    g_[a / 3][a % 3] = G_.g[a][l];
  for (int b = 0; b < 3; b++)
    g_[2][b] = G_.h[b][r];
  double ox = g_[2][1] - g_[1][2], oy = g_[0][2] - g_[2][0], oz = g_[1][0] - g_[0][1];
  double wn = sqrt(ox * ox + oy * oy + oz * oz), lm = lam2of(g_) / ls;
  fa[l] = (float)wn;
  fl[l] = (float)lm;
  if (i < ni) {
    long o = (k * ni + i) * ny + j;
    wa[o] = __double2half(wn);
    wz[o] = __double2half(oz);
    l2[o] = __double2half(lm);
  }
}
/* -V D: block mean over D x D in (x, y) of a full float field [i][kloc][j]: coarse node I (x = I D dx, 0 <= I <= (X1-1)/D)
   and J (y = -pi + J D dy) average the fine points within D/2 with weights 1/2 at the two ends (mirror in x at 0 and pi:
   |w| and lambda2 are even there; periodic in y); output half [kloc][I][J] */
__global__ static void wdown(long X1, long kc, long ny, long D, const float *f, __half *o) {
  long Xc = (X1 - 1) / D + 1, Yc = ny / D, l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= kc * Xc * Yc)
    return;
  long J = l % Yc, I = l / Yc % Xc, k = l / Yc / Xc, h = D / 2, n2 = 2 * (X1 - 1);
  double s = 0;
  for (long a = -h; a <= h; a++) {
    long i = I * D + a;
    i = i < 0 ? -i : i > X1 - 1 ? n2 - i : i;
    double wa = a == -h || a == h ? 0.5 : 1;
    for (long b = -h; b <= h; b++) {
      long j = ((J * D + b) % ny + ny) % ny;
      double wb = b == -h || b == h ? 0.5 : 1;
      s += wa * wb * f[(i * kc + k) * ny + j];
    }
  }
  o[l] = __float2half((float)(s / (D * D)));
}
/* -W snapshot: |w|, w_z, lambda2 / ls as half from the gradient planes i < ni of the slab ([i][kloc][j]), written
   in the order [kloc][i][j] (z, x, y) */
__global__ static void wsnap(long ni, long kc, long ny, Grad G_, double ls, __half *wa, __half *wz, __half *l2) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= ni * kc * ny)
    return;
  long j = l % ny, k = l / ny % kc, i = l / ny / kc, o = (k * ni + i) * ny + j;
  double g_[3][3];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      g_[a][b] = G_.g[a][b][l];
  double ox = g_[2][1] - g_[1][2], oy = g_[0][2] - g_[2][0], oz = g_[1][0] - g_[0][1];
  wa[o] = __double2half(sqrt(ox * ox + oy * oy + oz * oz));
  wz[o] = __double2half(oz);
  l2[o] = __double2half(lam2of(g_) / ls);
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
  int yao = 0, sflag = 0, kdump = 0; /* -I yao, -S, -k */
  long fg[3] = {0, 0, 0}; /* -G: grid of a k. restart file */
  double Wt[64], HW = 1.5, dD = 0, dA = 0, dB = 0; /* -W snapshot times, -w half width in x, -D interval,from,to */
  long VD = 0; /* -V: with -W also the whole quarter block-averaged by VD in x, y (z by the reader) */
  int nW = 0;
  double cfl = 0, dtmax, h, Gz = 0; /* -c */
  int dtknown = 0;
  char *tlist = NULL;
  long ne, nL, nT = 0, *Ts = NULL;
  double nu, dt, T_, t, dd, F, X, Y, Z, x, N, s, cx, cy, cz, kc, sums[6], mx[2],
      *host, *dP, *P, t0, t1, tio;
  double2 *u[3], *du[3], *b[7], *dg[6];
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
  syn6 = getenv("BOXSYM_SYN6") != NULL && atoi(getenv("BOXSYM_SYN6")) != 0;
  nocomm = getenv("BOXSYM_NOCOMM") != NULL && atoi(getenv("BOXSYM_NOCOMM")) != 0;
  nofuse = getenv("BOXSYM_NOFUSE") != NULL && atoi(getenv("BOXSYM_NOFUSE")) != 0;
  selfx = !usempi && !(getenv("BOXSYM_SELFCOPY") != NULL && atoi(getenv("BOXSYM_SELFCOPY")) != 0);
  MPI_Type_contiguous(2, MPI_DOUBLE, &cplx);
  MPI_Type_commit(&cplx);
  MPI_Type_contiguous(2, MPI_FLOAT, &cplxf);
  MPI_Type_commit(&cplxf);
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
                "  [-I yao (Yao & Hussain 2020a IC on the GPU, instead of -i/-r)] [-c <Courant number> (adaptive dt <= -s;\n"
                "  -d -e -L -T are then times, q.STEP.t holds t)] [-S (all-to-all in float)]\n"
                "  [-k (dumps are k.STEP: spectral, double complex [3][nx/2+1][nz/2+1][ny/2+1], unnormalized)]\n"
                "  [-r k.STEP [-G nx,ny,nz (its grid if not this one: zero-padded, its Nyquist modes dropped)]]\n"
                "  [-W t1,t2,... (with -c: snapshots w.STEP_{wabs,wz,lam2}.npy, half, [z][x][y] of the quarter, x < -w)]\n"
                "  [-w <x half width, default 1.5>] [-D dt,t0,t1 (with -c: diagnostics also at t0 + k dt in [t0, t1])]\n"
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
    if (argv[0][1] == 'S') {
      sflag = 1;
      continue;
    }
    if (argv[0][1] == 'k') {
      kdump = 1;
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
    if (argv[0][1] == 'W') {
      char *c = *++argv, *e;
      for (;;) {
        if (nW == 64)
          fail("-W: at most 64 times%s%s", "", "");
        Wt[nW++] = strtod(c, &e);
        if (e == c || (*e != ',' && *e != '\0'))
          fail("-W '%s' is not a list of times%s", *argv, "");
        if (*e == '\0')
          break;
        c = e + 1;
      }
      continue;
    }
    if (argv[0][1] == 'D') {
      if (sscanf(*++argv, "%lf,%lf,%lf", &dD, &dA, &dB) != 3 || dD <= 0)
        fail("-D '%s': need interval,from,to%s", *argv, "");
      continue;
    }
    if (argv[0][1] == 'G') {
      if (sscanf(*++argv, "%ld,%ld,%ld", &fg[0], &fg[1], &fg[2]) != 3)
        fail("-G '%s': need nx,ny,nz%s", *argv, "");
      continue;
    }
    if (argv[0][1] == 'I') {
      if (strcmp(*++argv, "yao") != 0)
        fail("-I '%s': only 'yao' is known%s", *argv, "");
      yao = 1;
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
    case 'c':
      cfl = x;
      break;
    case 'w':
      HW = x;
      break;
    case 'V':
      VD = x;
      break;
    default:
      fail("unknown option '%s'%s", *argv, "");
    }
    argv++;
  }
  if (cfl < 0 || (g.f ? 1.5 : 1.0) * cfl > 3.3)
    fail("-c: need 0 < C <= 3.3 (2.2 with -f)%s%s", "", "");
  if (nx < 4 || ny < 4 || nz < 4 || nx % 2 || ny % 2 || nz % 2 || nu < 0 || dt <= 0 ||
      T_ <= 0 || X <= 0 || Y <= 0 || Z <= 0 || (ipath != NULL) + (rpath != NULL) + yao != 1)
    fail("need -x -y -z (all even), -n, -s, -t and one of -i, -r, -I%s%s", "", "");
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
  /* BOXSYM_RB: the FFT buffer R holds one of nrb batches of x and z lines (less memory, more launches) */
  nrb = getenv("BOXSYM_RB") != NULL ? atoi(getenv("BOXSYM_RB")) : 1;
  if (nrb < 1)
    nrb = 1;
  if (nrb > 4)
    nrb = 4;
  if (nrb > g.pc)
    nrb = g.pc > 0 ? g.pc : 1;
  if (nrb > 1 && sflag) /* with float exchange buffers the received float field shares its buffer with the double
                          x-stage output, which needs all lines read first */
    fail("BOXSYM_RB > 1 needs the double exchange (no -S)%s%s", "", "");
  {
    long mx_ = 0, mz_ = 0;
    for (int k = 0; k < nrb; k++) {
      long a_ = first(g.kc * g.Q, nrb, k + 1) - first(g.kc * g.Q, nrb, k), b_ = first(g.pc, nrb, k + 1) - first(g.pc, nrb, k);
      mx_ = a_ > mx_ ? a_ : mx_;
      mz_ = b_ > mz_ ? b_ : mz_;
    }
    nR = nx * mx_ > nz * mz_ * g.Q ? nx * mx_ : nz * mz_ * g.Q;
  }
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
  dtmax = dt;
  nd = dd > 0 && cfl == 0 ? steps("-d", dd, dt) : 0;
  ne = ee > 0 && cfl == 0 ? steps("-e", ee, dt) : 0;
  nL = LL > 0 && cfl == 0 ? steps("-L", LL, dt) : 0;
  double *Tt = NULL; /* -T times (adaptive) */
  if (tlist != NULL) {
    char *c = tlist;
    if ((Ts = (long *)malloc((strlen(tlist) + 1) * sizeof(long))) == NULL)
      fail("malloc failed%s%s", "", "");
    for (;;) {
      double v = strtod(c, &end);
      if (end == c || (*end != ',' && *end != '\0'))
        fail("-T '%s' is not a list of times%s", tlist, "");
      if ((Tt = (double *)realloc(Tt, (nT + 1) * sizeof(double))) == NULL)
        fail("malloc failed%s%s", "", "");
      Tt[nT] = v;
      Ts[nT++] = cfl == 0 ? steps("-T", v, dt) : -1;
      if (*end == '\0')
        break;
      c = end + 1;
    }
  }
  nt = cfl == 0 ? steps("-t", T_, dt) : 0;
  for (int c = 0; c < 3; c++) {
    cuda(cudaMalloc(&u[c], g.nf * sizeof(double2)), "cudaMalloc");
    cuda(cudaMalloc(&du[c], g.nf * sizeof(double2)), "cudaMalloc");
    cuda(cudaMemset(du[c], 0, g.nf * sizeof(double2)), "cudaMemset");
  }
  for (int c = 0; c < 6; c++)
    cuda(cudaMalloc(&b[c], nb * sizeof(double2)), "cudaMalloc");
  /* b[6] (N_x of nl_a/nl_b, one more field); without it (BOXSYM_NL=0 or no memory) the stage is synth5, cross,
     analysis */
  nl = !(getenv("BOXSYM_NL") != NULL && atoi(getenv("BOXSYM_NL")) == 0);
  if (nl && cudaMalloc(&b[6], nb * sizeof(double2)) != cudaSuccess) {
    cudaGetLastError();
    nl = 0;
  }
  if (!nl)
    b[6] = NULL;
  cuda(cudaMalloc(&R, nR * sizeof(double2)), "cudaMalloc");
  for (int q = 0; q < 2; q++)
    cuda(cudaMalloc(&T[q], g.nf * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&dP, 6 * nblock * sizeof(double)), "cudaMalloc");
  {
    /* the separable mask, as mask(): exp(-36 (2p/nx)^36) per direction (-f) or 3p < nx */
    long nm = g.X1 + g.Q + g.Z1;
    double *hm = (double *)malloc(nm * sizeof(double)), *dm;
    if (hm == NULL)
      fail("malloc failed%s%s", "", "");
    for (long l = 0; l < nm; l++) {
      long p = l < g.X1 ? l : l < g.X1 + g.Q ? l - g.X1 : l - g.X1 - g.Q, n = l < g.X1 ? nx : l < g.X1 + g.Q ? ny : nz;
      double x = 2.0 * p / n;
      hm[l] = g.f ? exp(-36 * pow(x, 36)) : 3 * p < n;
    }
    cuda(cudaMalloc(&dm, nm * sizeof(double)), "cudaMalloc");
    cuda(cudaMemcpy(dm, hm, nm * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy");
    free(hm);
    g.mfx = dm;
    g.mfy = dm + g.X1;
    g.mfz = dm + g.X1 + g.Q;
  }
  if ((P = (double *)malloc(6 * nblock * sizeof(double))) == NULL ||
      (host = (double *)malloc(nloc * sizeof(double))) == NULL)
    fail("malloc failed%s%s", "", "");
  ws = 0;
  if (nrb == 1) {
    plan(&px, nx, nx, nx, nx, nx, 1, g.kc * g.Q, CUFFT_Z2Z, &ws);
    plan(&pz, nz, nz, nz, nz, nz, 1, g.pc * g.Q, CUFFT_Z2Z, &ws);
    pxb[0] = px;
    pzb[0] = pz;
  } else {
    for (int k = 0; k < nrb; k++) {
      plan(&pxb[k], nx, nx, nx, nx, nx, 1, xba(k + 1) - xba(k), CUFFT_Z2Z, &ws);
      plan(&pzb[k], nz, nz, nz, nz, nz, 1, (zba(k + 1) - zba(k)) * g.Q, CUFFT_Z2Z, &ws);
    }
    px = pxb[0];
    pz = pzb[0];
  }
  /* chunked stage boundaries (BOXSYM_HC): every rank needs at least nch p rows */
  nch = getenv("BOXSYM_HC") != NULL ? atoi(getenv("BOXSYM_HC")) : 1;
  if (nch < 1 || !nl)
    nch = 1;
  if (nch > 4)
    nch = 4;
  if (nch > g.X1 / size)
    nch = g.X1 / size;
  if (nrb > 1 && nch > 1 && nch < nrb) /* a chunk must fit R */
    nch = nrb;
  for (int j = 0; nch > 1 && j < nch; j++) {
    plan(&pzc[j], nz, nz, nz, nz, nz, 1, (chb(j) - cha(j)) * g.Q, CUFFT_Z2Z, &ws);
    if ((cnt_s[j] = (int *)malloc(4 * size * sizeof(int))) == NULL)
      fail("malloc failed%s%s", "", "");
    dsp_s[j] = cnt_s[j] + size;
    cnt_r[j] = dsp_s[j] + size;
    dsp_r[j] = cnt_r[j] + size;
    for (int q = 0; q < size; q++) {
      long kq = first(g.Z1, size, q + 1) - first(g.Z1, size, q), pq = first(g.X1, size, q + 1) - first(g.X1, size, q);
      long qa = first(pq, nch, j), qb = first(pq, nch, j + 1);
      cnt_s[j][q] = (chb(j) - cha(j)) * kq * g.Q;
      dsp_s[j][q] = sdz[q] + cha(j) * kq * g.Q;
      cnt_r[j][q] = (qb - qa) * g.kc * g.Q;
      dsp_r[j][q] = rdz[q] + qa * g.kc * g.Q;
    }
  }
  plan(&py[0], ny, 2 * g.Q, 2 * g.Q, g.Q, g.Q, 1, g.X1 * g.kc, CUFFT_D2Z, &ws);
  plan(&py[1], ny, g.Q, g.Q, 2 * g.Q, 2 * g.Q, 1, g.X1 * g.kc, CUFFT_Z2D, &ws);
  cuda(cudaMalloc(&work, ws > 0 ? ws : 1), "cudaMalloc");
  for (int q = 0; q < 4; q++) {
    cufft(cufftSetWorkArea(*all[q], work), "cufftSetWorkArea");
    cufft(cufftSetStream(*all[q], st), "cufftSetStream");
  }
  for (int k = 1; k < nrb; k++) {
    cufft(cufftSetWorkArea(pxb[k], work), "cufftSetWorkArea");
    cufft(cufftSetStream(pxb[k], st), "cufftSetStream");
    cufft(cufftSetWorkArea(pzb[k], work), "cufftSetWorkArea");
    cufft(cufftSetStream(pzb[k], st), "cufftSetStream");
  }
  for (int j = 0; nch > 1 && j < nch; j++) {
    cufft(cufftSetWorkArea(pzc[j], work), "cufftSetWorkArea");
    cufft(cufftSetStream(pzc[j], st), "cufftSetStream");
  }
  cuda(cudaMemGetInfo(&fr, &tot), "cudaMemGetInfo");
  char host_name[256] = "?";
  gethostname(host_name, sizeof host_name);
  fprintf(stderr,
          "boxsymmp: rank %d/%d %s gpu %d: p %ld+%ld k %ld+%ld, fields %.3f GB, R+T %.3f GB, "
          "work %.3f GB, gpu used %.3f GB\n",
          rank, size, host_name, local % ndev, g.p0, g.pc, g.k0, g.kc,
          ((6 + nl) * nb + 6 * g.nf) * sizeof(double2) / 1e9, (nR + 2 * g.nf) * sizeof(double2) / 1e9,
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
  if (yao) {
    yaow<<<grid(nloc), nthread, 0, st>>>(g, X * pi, Y * pi, Z * pi, (double *)b[0], (double *)b[1], (double *)b[2]);
    launched("yaow");
    analysis(3, b, wx, wz);
    biot<<<grid(g.nf), nthread, 0, st>>>(g, b[0], b[1], b[2], u[0], u[1], u[2]);
    launched("biot");
    devsync();
    name = (char *)"-I yao";
  } else if (strncmp(strrchr(name, '/') ? strrchr(name, '/') + 1 : name, "k.", 2) == 0) {
    long fx = fg[0] ? fg[0] : nx, fy = fg[0] ? fg[1] : ny, fz = fg[0] ? fg[2] : nz;
    long fX1 = fx / 2 + 1, fZ1 = fz / 2 + 1, fQ = fy / 2 + 1, same = fx == nx && fy == ny && fz == nz;
    long pk = same ? fX1 : fx / 2, rk = same ? fZ1 : fz / 2, qk = same ? fQ : fy / 2;
    double sc_ = N / ((double)fx * fy * fz); /* the file holds unnormalized FFT coefficients of its grid */
    if (fx > nx || fy > ny || fz > nz || fx % 2 || fy % 2 || fz % 2)
      fail("-G: the k. file grid must be even and not larger than -x -y -z%s%s", "", "");
    double2 *hk = (double2 *)malloc(g.nf * sizeof(double2)), *pl = (double2 *)malloc(fZ1 * fQ * sizeof(double2));
    if (hk == NULL || pl == NULL)
      fail("malloc failed%s%s", "", "");
    if ((fd = open(name, O_RDONLY)) < 0)
      fail("fail to open '%s'%s", name, "");
    if (fstat(fd, &sb) != 0 || sb.st_size != (off_t)(3 * fX1 * fZ1 * fQ * sizeof(double2)))
      fail("'%s' has wrong size (not a k. file of this grid or -G?)%s", name, "");
    for (int c = 0; c < 3; c++) {
      memset(hk, 0, g.nf * sizeof(double2));
      for (long p = g.p0; p < g.p0 + g.pc && p < pk; p++) {
        if (rw(0, fd, pl, fZ1 * fQ * sizeof(double2), (off_t)(((c * fX1 + p) * fZ1) * fQ * sizeof(double2))) != 0)
          fail("fail to read '%s'%s", name, "");
        for (long r = 0; r < rk; r++)
          for (long q = 0; q < qk; q++)
            hk[((p - g.p0) * g.Z1 + r) * g.Q + q] = make_double2(pl[r * fQ + q].x * sc_, pl[r * fQ + q].y * sc_);
      }
      cuda(cudaMemcpy(u[c], hk, g.nf * sizeof(double2), cudaMemcpyHostToDevice), "cudaMemcpy");
    }
    close(fd);
    free(hk);
    free(pl);
  } else {
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
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    fprintf(stderr, "boxsymmp: read '%s' in %.1f s\n", name, MPI_Wtime() - tio);
  start = tstep;
  t = tstep * dt;
  if (rpath != NULL) {
    /* q.STEP.t (written with every dump): t dt */
    FILE *ft;
    double tf, df;
    sprintf(path, "%s.t", rpath);
    if ((ft = fopen(path, "r")) != NULL) {
      if (fscanf(ft, "%lf %lf", &tf, &df) != 2)
        fail("cannot read '%s'%s", path, "");
      fclose(ft);
      if (cfl > 0) {
        t = tf;
        dt = fmin(df, dtmax);
        dtknown = fg[0] == 0 || (fg[0] == nx && fg[1] == ny && fg[2] == nz); /* a padded restart sets dt anew */
      } else if (fabs(tf - t) > 1e-9 * fmax(1, t))
        fail("'%s': t differs from STEP * -s (an adaptive run? use -c)%s", path, "");
    } else if (cfl > 0)
      fail("-c needs '%s' (t of the checkpoint)%s", path, "");
  }
  double *dG;
  cuda(cudaMalloc(&dG, sizeof(double)), "cudaMalloc");
  if (getenv("BOXSYM_A2AONLY") != NULL) {
    /* timing: K steps worth of exchanges alone (5 stages x (5 or 6 synthesis + 3 analysis)), same buffers and
       precision as the RK stages, back to back on the exchange stream */
    int K = atoi(getenv("BOXSYM_A2AONLY")), ns = syn6 ? 6 : 5;
    lowp = sflag;
    for (int it = 0; it <= K; it++) {
      if (it == 1) {
        cuda(cudaStreamSynchronize(sc), "cudaStreamSynchronize");
        MPI_Barrier(MPI_COMM_WORLD);
        t0 = MPI_Wtime();
      }
      for (int r = 0; r < 5; r++) {
        for (int c = 0; c < ns; c++)
          a2a(tb(c), scz, sdz, b[c], rcz, rdz);
        for (int c = 0; c < 3; c++)
          a2a(b[c], rcz, rdz, tb(c), scz, sdz);
      }
    }
    cuda(cudaStreamSynchronize(sc), "cudaStreamSynchronize");
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    lowp = 0;
    double by = (double)(ns + 3) * 5 * g.X1 * g.Z1 * g.Q * (sflag ? 8 : 16);
    if (rank == 0)
      fprintf(stderr, "boxsymmp: exchanges only (%d per step, %s): %.4f s/step, %.1f GB/step in all, %.1f GB/s per rank (incl. self)\n",
              (ns + 3) * 5, sflag ? "float" : "double", (t1 - t0) / K, by / 1e9, by / size / ((t1 - t0) / K) / 1e9);
  }
  devsync();
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for (;;) {
    int last = cfl == 0 ? tstep >= nt : t >= T_ - 1e-9 * fmax(1, T_);
    int diag = tstep % 10 == 0 || tstep == start || last, dump, espec, lam;
#define MULT(t, P) (fabs((t) - (P) * nearbyint((t) / (P))) <= 1e-9 * fmax(1, (t)))
    if (cfl == 0) {
      dump = nd > 0 && tstep % nd == 0;
      for (long q = 0; q < nT; q++)
        if (Ts[q] == tstep)
          dump = 1;
      espec = ne > 0 && tstep % ne == 0;
      lam = nL > 0 && tstep % nL == 0 && t >= L0 - 1e-9;
    } else {
      dump = dd > 0 && MULT(t, dd);
      for (long q = 0; q < nT; q++)
        if (fabs(t - Tt[q]) <= 1e-9 * fmax(1, t))
          dump = 1;
      espec = ee > 0 && MULT(t, ee);
      lam = LL > 0 && MULT(t, LL) && t >= L0 - 1e-9;
    }
    dump = dump && (rpath == NULL || tstep != start);
    int snap = 0;
    if (cfl > 0) {
      for (int q = 0; q < nW; q++)
        if (fabs(t - Wt[q]) <= 1e-9 * fmax(1, t))
          snap = rpath == NULL || tstep != start;
      if (dD > 0 && t >= dA - 1e-9 * fmax(1, t) && t <= dB + 1e-9 * fmax(1, t)) {
        double r_ = (t - dA) / dD;
        if (fabs(r_ - nearbyint(r_)) * dD <= 1e-9 * fmax(1, t))
          diag = 1;
      }
    }
    /* diagnostics on u^n: from the first RK stage's synthesis (fuse), or (dump, last step) a synthesis of their own */
    int fuse = diag && !dump && !last && !nofuse;
    auto diagnose = [&]() {
      norms<<<nblock, nthread, 0, st>>>(g, u[0], u[1], u[2], dP);
      launched("norms");
      collect(6, 0, dP, P, sums);
      peaks<<<nblock, nthread, 0, st>>>(g, s, cx, cy, cz, (double *)dg[0], (double *)dg[1],
                                 (double *)dg[2], (double *)dg[3], (double *)dg[4],
                                 (double *)dg[5], dP);
      launched("peaks");
      collect(2, 1, dP, P, mx);
      cuda(cudaMemsetAsync(dG, 0, sizeof(double), st), "cudaMemsetAsync");
      if (g.k0 == 0 && g.kc > 0)
        gzsum<<<1, nthread, 0, st>>>(g, (const double *)dg[1], dG);
      launched("gzsum");
      devsync();
      cuda(cudaMemcpy(&Gz, dG, sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
      MPI_Allreduce(MPI_IN_PLACE, &Gz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      Gz *= s * Y * pi / ny;
      if (cfl > 0 && !nocomm && (tstep % 10 == 0 || !dtknown)) {
        /* the same rule in a continuous and a restarted run (dt of the checkpoint from q.STEP.t); the out
           line shows the dt of the next steps and its C */
        double dn = cfl / mx[0];
        dt = fmin(dtmax, fmin(dn, dtknown ? 1.1 * dt : dn));
        dtknown = 1;
      }
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
                    sums[5] / sums[0],
                    Gz,
                    dt};
      char sv[64];
      if (!nocomm && !(isfinite(E) && isfinite(Om) && isfinite(mx[0]) && isfinite(mx[1]))) {
        snprintf(sv, sizeof sv, "%ld", tstep);
        fail("non-finite E, Omega, C or wmax at step %s%s", sv, "");
      }
      if (!nocomm && (g.f ? 1.5 : 1.0) * dt * mx[0] > 3.3) {
        snprintf(sv, sizeof sv, "%ld (C = %g)", tstep, dt * mx[0]);
        fail("Courant number too large at step %s%s: reduce -s", sv, "");
      }
      if (rank == 0) {
        if (tstep == start)
          printf("step t E Omega C wmax keta Rlambda lambda L tailx taily tailz Gz dt\n");
        printf("% 10ld", tstep);
        for (int i = 0; i < (int)(sizeof q / sizeof *q); i++)
          printf(" % .16e", q[i]);
        printf("\n");
        if (fflush(stdout) != 0 || ferror(stdout))
          fail("fail to write stdout%s%s", "", "");
      }
    };
    if ((diag && !fuse) || dump) {
      spread<<<grid(g.nf), nthread, 0, st>>>(g, u[0], u[1], u[2], b[0], b[1], b[2], b[3], b[4], b[5]);
      launched("spread");
      synthesis(6, b, ox6, oz6);
    }
    if (diag && !fuse) {
      for (int c = 0; c < 6; c++)
        dg[c] = b[c];
      diagnose();
    }
    if (dump && kdump) {
      /* k.STEP: rank's p rows are one contiguous block per component */
      sprintf(path, "k.%08ld.tmp", tstep);
      tio = MPI_Wtime();
      bad = 0;
      double2 *hk = (double2 *)malloc(g.nf * sizeof(double2));
      if (hk == NULL)
        fail("malloc failed%s%s", "", "");
      if (rank == 0 &&
          ((fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644)) < 0 ||
           ftruncate(fd, (off_t)(3 * g.X1 * g.Z1 * g.Q * sizeof(double2))) != 0 || close(fd) != 0))
        bad = 1;
      if (any(bad))
        fail("fail to open '%s'%s", path, "");
      if ((fd = open(path, O_WRONLY)) < 0)
        bad = 1;
      for (int c = 0; c < 3 && !bad; c++) {
        cuda(cudaMemcpy(hk, u[c], g.nf * sizeof(double2), cudaMemcpyDeviceToHost), "cudaMemcpy");
        if (g.nf > 0 && rw(1, fd, hk, g.nf * sizeof(double2), (off_t)((c * g.X1 + g.p0) * g.Z1 * g.Q * sizeof(double2))) != 0)
          bad = 1;
      }
      free(hk);
      if (fd >= 0 && (fsync(fd) != 0 || close(fd) != 0))
        bad = 1;
      if (any(bad))
        fail("fail to write '%s'%s", path, "");
      char final[FILENAME_MAX];
      sprintf(final, "k.%08ld", tstep);
      if (rank == 0 && rename(path, final) != 0)
        bad = 1;
      if (any(bad))
        fail("fail to rename '%s'%s", path, "");
      if (rank == 0) {
        FILE *ft;
        sprintf(path, "k.%08ld.t", tstep);
        if ((ft = fopen(path, "w")) == NULL || fprintf(ft, "%.17g %.17g\n", t, dt) < 0 || fclose(ft) != 0)
          bad = 1;
      }
      if (any(bad))
        fail("fail to write '%s'%s", path, "");
      if (rank == 0)
        fprintf(stderr, "boxsymmp: wrote '%s' in %.1f s\n", final, MPI_Wtime() - tio);
    }
    if (dump && !kdump) {
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
      if (rank == 0) {
        FILE *ft;
        sprintf(path, "q.%08ld.t", tstep);
        if ((ft = fopen(path, "w")) == NULL || fprintf(ft, "%.17g %.17g\n", t, dt) < 0 || fclose(ft) != 0)
          bad = 1;
      }
      if (any(bad))
        fail("fail to write '%s'%s", path, "");
      if (rank == 0)
        fprintf(stderr, "boxsymmp: wrote '%s' in %.1f s\n", final, MPI_Wtime() - tio);
    }
    if (espec) {
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
    if (snap) {
      /* velocity gradient for the planes i <= imax (x <= HW) of every slab, as the lambda2 crop; then |w|, w_z,
         lambda2 / 1000 as half in w.STEP_{wabs,wz,lam2}.npy, shape (nz/2+1, imax+1, ny), each rank one block */
      double ts_ = MPI_Wtime(), dx_ = X * pi / nx, ls_ = 1000;
      long ni = (long)(HW / dx_) + 1;
      if (ni > g.X1)
        ni = g.X1;
      long nl = ni * g.kc * ny;
      float *LG[3][3];
      __half *wh;
      /* the gradient planes live in du (dead at the start of a step: the first RK stage has A = 0; zeroed again
         below) and the half output in R (free after the syntheses) when they fit, else in their own memory */
      int own = VD == 0 && (3 * (size_t)nl * sizeof(float) > (size_t)g.nf * sizeof(double2) ||
                            3 * (size_t)nl * sizeof(__half) > (size_t)nR * sizeof(double2));
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          if (own)
            cuda(cudaMalloc(&LG[i][j], (nl > 0 ? nl : 1) * sizeof(float)), "cudaMalloc");
          else
            LG[i][j] = (float *)du[i] + (size_t)j * nl; /* three gradient planes per du[i] */
        }
      if (own)
        cuda(cudaMalloc(&wh, 3 * (nl > 0 ? nl : 1) * sizeof(__half)), "cudaMalloc");
      else
        wh = (__half *)R;
      long nfull = g.X1 * g.kc * ny, Xc = VD > 0 ? (g.X1 - 1) / VD + 1 : 0, Yc = VD > 0 ? ny / VD : 0, nd = g.kc * Xc * Yc;
      __half *wd = NULL;
      if (VD > 0) {
        /* six gradient components in float, two per du[c] (2 nfull <= 4 nf floats), the last three stay double in b0..b2;
           crop halves and the coarse output in b3, full float |w| and lambda2 in b4, b5 */
        if (2 * (size_t)nfull * sizeof(float) > (size_t)g.nf * sizeof(double2) || (size_t)nfull * sizeof(float) > (size_t)nb * sizeof(double2) ||
            (3 * (size_t)nl + 2 * (size_t)nd) * sizeof(__half) > (size_t)nb * sizeof(double2) || nx % (2 * VD) || ny % VD)
          fail("-V: does not fit or does not divide the grid%s%s", "", "");
        for (int a = 0; a < 6; a++)
          LG[a / 3][a % 3] = (float *)du[a / 2] + (size_t)(a % 2) * nfull;
        wh = (__half *)b[3];
        wd = wh + 3 * nl;
      }
      for (int hh = 0; hh < 2; hh++) {
        int m = hh == 0 ? 6 : 3, ox9[6], oz9[6];
        for (int q = 0; q < m; q++) {
          int ij = 6 * hh + q, i = ij / 3, j = ij % 3;
          lderiv<<<grid(g.nf), nthread, 0, st>>>(g, j, u[i], b[q], s);
          launched("lderiv");
          ox9[q] = ux[i] ^ (j == 0);
          oz9[q] = uz[i] ^ (j == 2);
        }
        synthesis(m, b, ox9, oz9);
        for (int q = 0; q < m && nl > 0 && !(VD > 0 && hh == 1); q++) {
          int ij = 6 * hh + q;
          long np_ = VD > 0 ? nfull : nl;
          ltofloat<<<grid(np_), nthread, 0, st>>>(np_ / ny, ny, g.Q, (double *)b[q], LG[ij / 3][ij % 3]);
          launched("ltofloat");
        }
      }
      if (VD > 0) {
        Grad6 g6;
        for (int a = 0; a < 6; a++)
          g6.g[a] = LG[a / 3][a % 3];
        for (int c = 0; c < 3; c++)
          g6.h[c] = (const double *)b[c];
        wsnap2<<<grid(nfull), nthread, 0, st>>>(g.X1, g.kc, ny, g.Q, ni, g6, ls_, wh, wh + nl, wh + 2 * nl, (float *)b[4],
                                                (float *)b[5]);
        launched("wsnap2");
        wdown<<<grid(nd), nthread, 0, st>>>(g.X1, g.kc, ny, VD, (const float *)b[4], wd);
        wdown<<<grid(nd), nthread, 0, st>>>(g.X1, g.kc, ny, VD, (const float *)b[5], wd + nd);
        launched("wdown");
      }
      Grad gg;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          gg.g[i][j] = LG[i][j];
      if (nl > 0 && VD == 0)
        wsnap<<<grid(nl), nthread, 0, st>>>(ni, g.kc, ny, gg, ls_, wh, wh + nl, wh + 2 * nl);
      launched("wsnap");
      devsync();
      __half *hh_ = (__half *)malloc(3 * (nl > 0 ? nl : 1) * sizeof(__half));
      if (hh_ == NULL)
        fail("malloc failed%s%s", "", "");
      cuda(cudaMemcpy(hh_, wh, 3 * nl * sizeof(__half), cudaMemcpyDeviceToHost), "cudaMemcpy");
      const char *nm_[3] = {"wabs", "wz", "lam2"};
      char hd[128], shp[64], wp[FILENAME_MAX];
      snprintf(shp, sizeof shp, "(%ld, %ld, %ld)", g.Z1, ni, ny);
      int hn = npyhead(hd, "<f2", shp);
      bad = 0;
      for (int q = 0; q < 3; q++) {
        sprintf(wp, "w.%08ld_%s.npy", tstep, nm_[q]);
        if (rank == 0 && ((fd = open(wp, O_WRONLY | O_CREAT | O_TRUNC, 0644)) < 0 ||
                          ftruncate(fd, (off_t)(hn + 2 * g.Z1 * ni * ny)) != 0 || rw(1, fd, hd, hn, 0) != 0 ||
                          close(fd) != 0))
          bad = 1;
        if (any(bad))
          fail("fail to write '%s'%s", wp, "");
        if ((fd = open(wp, O_WRONLY)) < 0 ||
            (nl > 0 && rw(1, fd, hh_ + q * nl, 2 * nl, (off_t)(hn + 2 * g.k0 * ni * ny)) != 0) || close(fd) != 0)
          bad = 1;
        if (any(bad))
          fail("fail to write '%s'%s", wp, "");
      }
      if (VD > 0) {
        __half *hd_ = (__half *)malloc(2 * (nd > 0 ? nd : 1) * sizeof(__half));
        if (hd_ == NULL)
          fail("malloc failed%s%s", "", "");
        cuda(cudaMemcpy(hd_, wd, 2 * nd * sizeof(__half), cudaMemcpyDeviceToHost), "cudaMemcpy");
        const char *dn_[2] = {"wabs", "lam2"};
        snprintf(shp, sizeof shp, "(%ld, %ld, %ld)", g.Z1, Xc, Yc);
        int hn2 = npyhead(hd, "<f2", shp);
        for (int q = 0; q < 2; q++) {
          sprintf(wp, "d.%08ld_%s.npy", tstep, dn_[q]);
          if (rank == 0 && ((fd = open(wp, O_WRONLY | O_CREAT | O_TRUNC, 0644)) < 0 ||
                            ftruncate(fd, (off_t)(hn2 + 2 * g.Z1 * Xc * Yc)) != 0 || rw(1, fd, hd, hn2, 0) != 0 || close(fd) != 0))
            bad = 1;
          if (any(bad))
            fail("fail to write '%s'%s", wp, "");
          if ((fd = open(wp, O_WRONLY)) < 0 ||
              (nd > 0 && rw(1, fd, hd_ + q * nd, 2 * nd, (off_t)(hn2 + 2 * g.k0 * Xc * Yc)) != 0) || close(fd) != 0)
            bad = 1;
          if (any(bad))
            fail("fail to write '%s'%s", wp, "");
        }
        free(hd_);
      }
      if (rank == 0) {
        FILE *fr;
        sprintf(wp, "w.%08ld.txt", tstep);
        if ((fr = fopen(wp, "w")) != NULL) {
          fprintf(fr, "boxsymmp snapshot, step %ld, t = %.10g, t* = t Gamma0 / (2 pi b^2) = %.6f (b 1.62, Gamma0 7.680385)\n"
                      "files w.%08ld_{wabs,wz,lam2}.npy: float16, shape (%ld, %ld, %ld), C order [k][i][j]\n"
                      "  k: z = k dz, k = 0 .. nz/2 (z in [0, pi], the tube axis), dz = %.10g\n"
                      "  i: x = i dx, i = 0 .. %ld (x in [0, %.6g]: collision normal, x = 0 is the collision plane), dx = %.10g\n"
                      "  j: y = -pi + j dy, j = 0 .. ny - 1 (propagation direction), dy = %.10g\n"
                      "wabs = |omega|; wz = omega_z (signed); lam2 = lambda2 / %g (lambda2 of S^2 + W^2, negative in vortex cores)\n"
                      "Quarter of the full box: the full 2 pi box follows by the mirror symmetries\n"
                      "  x -> -x: |w|, lambda2 even; omega_z odd (the tube at x < 0 has omega_z > 0)\n"
                      "  z -> -z: |w|, lambda2, omega_z even\n"
                      "grid %ld^3 (2 pi periodic), Re = Gamma0 / nu = %.6g\n",
                  tstep, t, t * 7.680385244 / (2 * pi * 1.62 * 1.62), tstep, g.Z1, ni, ny, Z * pi / nz, ni - 1,
                  (ni - 1) * dx_, dx_, Y * pi / ny, ls_, nx, 7.680385244 / nu);
          if (VD > 0)
            fprintf(fr, "also d.%08ld_{wabs,lam2}.npy: the whole quarter, float16, shape (%ld, %ld, %ld), order [k][I][J]: k as above (all\n"
                        "  z planes), x = I %ld dx, y = -pi + J %ld dy; block mean over %ld x %ld fine points in (x, y) centred on the node,\n"
                        "  weights 1/2 at the two ends, mirrored at x = 0 and pi; average k the same way to get the 3D block mean\n",
                    tstep, g.Z1, (g.X1 - 1) / VD + 1, ny / VD, VD, VD, VD, VD);
          fclose(fr);
        }
        fprintf(stderr, "boxsymmp: snapshot w.%08ld (%ld x %ld x %ld, 3 x %.2f GB) in %.1f s\n", tstep, g.Z1, ni, ny,
                2e-9 * g.Z1 * ni * ny, MPI_Wtime() - ts_);
      }
      free(hh_);
      if (own) {
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            cuda(cudaFree(LG[i][j]), "cudaFree");
        cuda(cudaFree(wh), "cudaFree");
      }
      for (int c = 0; c < 3; c++)
        cuda(cudaMemsetAsync(du[c], 0, g.nf * sizeof(double2), st), "cudaMemsetAsync");
    }
    if (lam && (rpath == NULL || tstep != start || L0start)) {
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
    if (last)
      break;
    /* step size: dt, or (adaptive) cut to land on the next output time */
    double te = T_;
    int land = 0;
    auto steph = [&]() {
    h = dt;
    if (cfl > 0) {
      double P_[3] = {dd, ee, LL};
      for (int q = 0; q < 3; q++)
        if (P_[q] > 0)
          te = fmin(te, (floor(t / P_[q] + 1e-9) + 1) * P_[q]);
      for (long q = 0; q < nT; q++)
        if (Tt[q] > t + 1e-9 * fmax(1, t))
          te = fmin(te, Tt[q]);
      for (int q = 0; q < nW; q++)
        if (Wt[q] > t + 1e-9 * fmax(1, t))
          te = fmin(te, Wt[q]);
      if (dD > 0 && t < dB - 1e-9 * fmax(1, t))
        te = fmin(te, t < dA - 1e-9 * fmax(1, t) ? dA : dA + (floor((t - dA) / dD + 1e-9) + 1) * dD);
      if (te - t <= h * (1 + 1e-9)) {
        h = te - t;
        land = 1;
      }
    }
    g.fp = F > 0 ? h / F : 1;
    };
    if (!fuse)
      steph();
    lowp = sflag;
    for (int r = 0; r < 5; r++) {
      if (syn6) {
        spread<<<grid(g.nf), nthread, 0, st>>>(g, u[0], u[1], u[2], b[0], b[1], b[2], b[3], b[4], b[5]);
        launched("spread");
        synthesis(6, b, ox6, oz6);
      } else {
        if (nl) {
          if (!headdone)
            spreadw<<<grid(g.nf), nthread, 0, st>>>(g, u[0], u[1], u[2], B_(OX_), B_(OY_), 0, g.nf);
          launched("spreadw");
          nl_a(b, u, s);
        } else {
          spread5<<<grid(g.nf), nthread, 0, st>>>(g, u[0], u[1], u[2], b[0], b[1], b[2], b[3], b[4], 0, g.nf);
          launched("spread5");
          synth5(b);
        }
      }
      if (r == 0 && fuse) {
        for (int c = 0; c < 6; c++)
          dg[c] = nl && !syn6 ? B_(c) : b[c];
        diagnose();
        steph();
      }
      if (nl && !syn6) {
        int head = nch > 1 && r < 4; /* chunked boundary to the next stage of this step */
        headdone = 0;
        nl_b(b, s, head);
        if (!head) {
          update<<<grid(g.nf), nthread, 0, st>>>(g, r == 4, nu, h, A[r], B[r], B_(N0_), B_(U_), B_(V_), u[0], u[1],
                                                 u[2], du[0], du[1], du[2], 0, g.nf);
          launched("update");
          continue;
        }
        /* next stage's map: U' receives into the dead W, OX' = old OZ (dead), OY' = old N0 (rows already used) */
        int nm[7] = {bm[W_], bm[OX_], bm[OY_], bm[OZ_], bm[N0_], bm[U_], bm[V_]};
        for (int j = 0; j < nch; j++) {
          long l0 = cha(j) * g.Z1 * g.Q, l1 = chb(j) * g.Z1 * g.Q;
          azstagec(7, B_(V_), uz[2], j);
          update<<<grid(l1 - l0), nthread, 0, st>>>(g, 0, nu, h, A[r], B[r], B_(N0_), B_(U_), B_(V_), u[0], u[1],
                                                    u[2], du[0], du[1], du[2], l0, l1);
          launched("update");
          spreadw<<<grid(l1 - l0), nthread, 0, st>>>(g, u[0], u[1], u[2], b[nm[OX_]], b[nm[OY_]], l0, l1);
          launched("spreadw");
          zsendc(0, u[0], b[nm[U_]], oz6[0], j);
        }
        memcpy(bm, nm, sizeof bm);
        headdone = 1;
        continue;
      }
      cross<<<grid(g.nr), nthread, 0, st>>>(g, s, (double *)b[0], (double *)b[1], (double *)b[2],
                                     (double *)b[3], (double *)b[4], (double *)b[5]);
      launched("cross");
      analysis(3, b, ux, uz);
      update<<<grid(g.nf), nthread, 0, st>>>(g, r == 4, nu, h, A[r], B[r], b[0], b[1], b[2], u[0],
                                      u[1], u[2], du[0], du[1], du[2], 0, g.nf);
      launched("update");
    }
    lowp = 0;
    tstep++;
    t = cfl == 0 ? tstep * dt : land ? te : t + h;
  }
  devsync();
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if (rank == 0)
    fprintf(stderr, "boxsymmp: %ld steps on %d ranks in %.2f s: %.4f s/step\n", tstep - start, size,
            t1 - t0, tstep > start ? (t1 - t0) / (tstep - start) : 0.0);
  if (!usempi)
    ncclCommDestroy(nccl);
  MPI_Finalize();
}
