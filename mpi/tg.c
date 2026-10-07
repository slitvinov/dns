#include <fftw3.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
enum { SIN, COS };
enum { EVEN, ODD };
static long N, N1, n3, nl, s0, nloc;
static int rank, size, *cnt, *off, *sc, *sd;
static fftw_plan plans[2][2][2][2][2][2];
static double *wsyn[2][2], *wana[2][2], *xa, *ya, *buf;
static void parallel_loop(void *(*work)(char *), char *jobdata, size_t elsize,
                          int njobs, void *data) {
  (void)data;
#pragma omp parallel for
  for (int i = 0; i < njobs; ++i)
    work(jobdata + elsize * i);
}
static long xidx(long i, long j, long k) { return (i * N1 + j) * N1 + k; }
static long yidx(long i, long j, long k) { return (i * nl + j) * N1 + k; }
static int wavenumber(int c, long n) { return c == EVEN ? 2 * n : 2 * n + 1; }
static void kind(int dir, int c, int type, fftw_r2r_kind *kd, long *n,
                 long *coff, long *poff) {
  if (c == EVEN && type == COS) {
    *kd = FFTW_REDFT00;
    *n = N + 1;
    *coff = 0;
    *poff = 0;
  } else if (c == EVEN && type == SIN) {
    *kd = FFTW_RODFT00;
    *n = N - 1;
    *coff = 1;
    *poff = 1;
  } else if (c == ODD && type == COS) {
    *kd = dir == 0 ? FFTW_REDFT10 : FFTW_REDFT01;
    *n = N;
    *coff = 0;
    *poff = 0;
  } else {
    *kd = dir == 0 ? FFTW_RODFT10 : FFTW_RODFT01;
    *n = N;
    *coff = 0;
    *poff = 1;
  }
}
static fftw_plan plan(int rnk, fftw_iodim *dims, int hrnk, fftw_iodim *hdims,
                      long in, long out, fftw_r2r_kind *kd) {
  double *x = fftw_alloc_real(nloc), *y = fftw_alloc_real(nloc);
  fftw_plan p = fftw_plan_guru_r2r(rnk, dims, hrnk, hdims, x + in, y + out, kd,
                                   FFTW_MEASURE | FFTW_UNALIGNED);
  if (p == NULL) {
    fprintf(stderr, "tg: error: fftw_plan_guru_r2r failed\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  fftw_free(x);
  fftw_free(y);
  return p;
}
static void zero(double *a) {
#pragma omp parallel for
  for (long l = 0; l < nloc; l++)
    a[l] = 0;
}
static void transpose(int dir, double *a, double *b) {
  long i, j;
  int r;
  if (dir == 0)
    MPI_Alltoallv(a, sc, sd, MPI_DOUBLE, buf, sc, sd, MPI_DOUBLE,
                  MPI_COMM_WORLD);
#pragma omp parallel for private(j, r)
  for (i = 0; i < nl; i++)
    for (r = 0; r < size; r++)
      for (j = 0; j < cnt[r]; j++)
        if (dir == 0)
          memcpy(&b[xidx(i, off[r] + j, 0)], &buf[sd[r] + (i * cnt[r] + j) * N1],
                 N1 * sizeof(double));
        else
          memcpy(&buf[sd[r] + (i * cnt[r] + j) * N1], &a[xidx(i, off[r] + j, 0)],
                 N1 * sizeof(double));
  if (dir == 1)
    MPI_Alltoallv(buf, sc, sd, MPI_DOUBLE, b, sc, sd, MPI_DOUBLE,
                  MPI_COMM_WORLD);
}
static void transform(int dir, int c, const int *t, double *a, double *b) {
  long n[3], coff[3], poff[3], cy, py, cx, px;
  fftw_iodim d1[2], h1[1], d2[1], h2[2];
  fftw_r2r_kind kd[3], k1[2];
  fftw_plan *p;
  for (int d = 0; d < 3; d++)
    kind(dir, c, t[d], &kd[d], &n[d], &coff[d], &poff[d]);
  cy = coff[0] * nl * N1 + coff[2];
  py = poff[0] * nl * N1 + poff[2];
  cx = coff[1] * N1;
  px = poff[1] * N1;
  p = plans[dir][c][t[0]][t[1]][t[2]];
  if (p[0] == NULL) {
    d1[0].n = n[0];
    d1[0].is = d1[0].os = nl * N1;
    d1[1].n = n[2];
    d1[1].is = d1[1].os = 1;
    h1[0].n = nl;
    h1[0].is = h1[0].os = N1;
    k1[0] = kd[0];
    k1[1] = kd[2];
    d2[0].n = n[1];
    d2[0].is = d2[0].os = N1;
    h2[0].n = nl;
    h2[0].is = h2[0].os = N1 * N1;
    h2[1].n = N1;
    h2[1].is = h2[1].os = 1;
    p[0] = plan(2, d1, 1, h1, dir == 0 ? cy : py, dir == 0 ? py : cy, k1);
    p[1] = plan(1, d2, 2, h2, dir == 0 ? cx : px, dir == 0 ? px : cx, &kd[1]);
  }
  if (dir == 0) {
    zero(ya);
    fftw_execute_r2r(p[0], a + cy, ya + py);
    transpose(0, ya, xa);
    zero(b);
    fftw_execute_r2r(p[1], xa + cx, b + px);
  } else {
    zero(xa);
    fftw_execute_r2r(p[1], a + px, xa + cx);
    transpose(1, xa, ya);
    zero(b);
    fftw_execute_r2r(p[0], ya + py, b + cy);
  }
}
static void scale(double **w, const int *t, double *a) {
  long i, j, k;
#pragma omp parallel for private(j, k)
  for (i = 0; i < N1; i++)
    for (j = 0; j < nl; j++)
      for (k = 0; k < N1; k++)
        a[yidx(i, j, k)] *= w[t[0]][i] * w[t[1]][s0 + j] * w[t[2]][k];
}
static void synthesis(int c, const int *t, const double *u, double *tmp,
                      double *f) {
#pragma omp parallel for
  for (long l = 0; l < nloc; l++)
    tmp[l] = u[l];
  scale(wsyn[c], t, tmp);
  transform(0, c, t, tmp, f);
}
static void analysis(int c, const int *t, double *f, double *u) {
  transform(1, c, t, f, u);
  scale(wana[c], t, u);
}
static void moments(double *fe, double *fo, double *S) {
  long i, j, k, l;
  double m[9] = {0}, w, x, y, z;
#pragma omp parallel for private(j, k, l, w, x, y, z) reduction(+ : m[:9])
  for (i = 0; i < nl; i++)
    for (j = 0; j < N1; j++)
      for (k = 0; k < N1; k++) {
        w = (s0 + i == 0 || s0 + i == N ? 0.5 : 1) *
            (j == 0 || j == N ? 0.5 : 1) * (k == 0 || k == N ? 0.5 : 1);
        l = xidx(i, j, k);
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
  MPI_Allreduce(MPI_IN_PLACE, m, 9, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (int n = 3; n < 9; n++)
    S[n] = (n % 2 ? -1 : 1) * (m[n] / m[0]) / pow(m[2] / m[0], n / 2.0);
}
int main(int argc, char **argv) {
  long M, i, j, k, l, tstep, ne, nd;
  int c, d, kcut, provided, r;
  double nu, dt, T, t, x, e, dd, *u[2][3], *um[2][3], *w[2][3], *du[2][3],
      *U[2][3], *W[2][3], *F[3], *tmp;
  char *end;
  static const int tu[3][3] = {
      {SIN, COS, COS}, {COS, SIN, COS}, {COS, COS, SIN}};
  static const int tw[3][3] = {
      {COS, SIN, SIN}, {SIN, COS, SIN}, {SIN, SIN, COS}};
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (provided < MPI_THREAD_FUNNELED) {
    if (rank == 0)
      fprintf(stderr, "tg: error: MPI_THREAD_FUNNELED is not supported\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  M = 0;
  nu = -1;
  dt = -1;
  T = 0;
  e = 0;
  dd = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    if (argv[0][1] == 'h') {
      if (rank == 0)
        fprintf(stderr,
                "Usage: tg -M <modes> -n <viscosity> -t <end time> -s "
                "<time step> [-e <spectrum interval>] [-d <dump interval>]\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (argv[1] == NULL) {
      if (rank == 0)
        fprintf(stderr, "tg: error: %s needs an argument\n", argv[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    x = strtod(argv[1], &end);
    if (*end != '\0') {
      if (rank == 0)
        fprintf(stderr, "tg: error: '%s' is not a number\n", argv[1]);
      MPI_Abort(MPI_COMM_WORLD, 1);
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
      T = x;
      break;
    case 'e':
      e = x;
      break;
    case 'd':
      dd = x;
      break;
    default:
      if (rank == 0)
        fprintf(stderr, "tg: error: unknown option '%s'\n", *argv);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    argv++;
  }
  if (M < 4 || M % 2 != 0 || nu < 0 || dt <= 0 || T <= 0) {
    if (rank == 0)
      fprintf(stderr, "tg: error: need -M (even, >= 4), -n, -s, -t\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  N = M / 2;
  N1 = N + 1;
  n3 = N1 * N1 * N1;
  if (N1 < size) {
    if (rank == 0)
      fprintf(stderr, "tg: error: more ranks than M/2 + 1\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
  fftw_threads_set_callback(parallel_loop, NULL);
  cnt = malloc(size * sizeof(int));
  off = malloc(size * sizeof(int));
  sc = malloc(size * sizeof(int));
  sd = malloc(size * sizeof(int));
  for (r = 0; r < size; r++) {
    cnt[r] = N1 / size + (r < N1 % size);
    off[r] = r == 0 ? 0 : off[r - 1] + cnt[r - 1];
  }
  nl = cnt[rank];
  s0 = off[rank];
  nloc = nl * N1 * N1;
  for (r = 0; r < size; r++) {
    sc[r] = nl * cnt[r] * N1;
    sd[r] = nl * off[r] * N1;
  }
  kcut = M / 3;
  ne = e > 0 ? lround(e / dt) : 0;
  nd = dd > 0 ? lround(dd / dt) : 0;
  for (c = 0; c < 2; c++)
    for (d = 0; d < 2; d++) {
      wsyn[c][d] = malloc(N1 * sizeof(double));
      wana[c][d] = malloc(N1 * sizeof(double));
      for (i = 0; i < N1; i++) {
        double s = c == EVEN && d == COS && (i == 0 || i == N) ? 1 : 0.5;
        wsyn[c][d][i] = s;
        wana[c][d][i] = 1 / (2 * N * s);
      }
    }
  for (c = 0; c < 2; c++)
    for (d = 0; d < 3; d++) {
      u[c][d] = fftw_alloc_real(nloc);
      um[c][d] = fftw_alloc_real(nloc);
      w[c][d] = fftw_alloc_real(nloc);
      du[c][d] = fftw_alloc_real(nloc);
      U[c][d] = fftw_alloc_real(nloc);
      W[c][d] = fftw_alloc_real(nloc);
      memset(u[c][d], 0, nloc * sizeof(double));
    }
  for (d = 0; d < 3; d++)
    F[d] = fftw_alloc_real(nloc);
  tmp = fftw_alloc_real(nloc);
  xa = fftw_alloc_real(nloc);
  ya = fftw_alloc_real(nloc);
  buf = fftw_alloc_real(nloc);
  if (rank == 0) {
    u[ODD][0][yidx(0, 0, 0)] = 1;
    u[ODD][1][yidx(0, 0, 0)] = -1;
  }
  t = 0;
  tstep = 0;
  for (;;) {
    for (c = 0; c < 2; c++)
#pragma omp parallel for private(j, k, l)
      for (i = 0; i < N1; i++)
        for (j = 0; j < nl; j++)
          for (k = 0; k < N1; k++) {
            double m = wavenumber(c, i), n = wavenumber(c, s0 + j),
                   p = wavenumber(c, k);
            l = yidx(i, j, k);
            w[c][0][l] = -n * u[c][2][l] + p * u[c][1][l];
            w[c][1][l] = -p * u[c][0][l] + m * u[c][2][l];
            w[c][2][l] = -m * u[c][1][l] + n * u[c][0][l];
          }
    if (tstep % 10 == 0) {
      double energy = 0, Omega = 0, Pal = 0, sum[3];
      for (c = 0; c < 2; c++)
#pragma omp parallel for private(j, k, l, d) reduction(+ : energy, Omega, Pal)
        for (i = 0; i < N1; i++)
          for (j = 0; j < nl; j++)
            for (k = 0; k < N1; k++) {
              long g[3] = {wavenumber(c, i), wavenumber(c, s0 + j),
                           wavenumber(c, k)};
              l = yidx(i, j, k);
              for (d = 0; d < 3; d++) {
                double su = 1, sw = 1;
                for (int e = 0; e < 3; e++) {
                  su *= tu[d][e] == COS && g[e] == 0 ? 1 : 0.5;
                  sw *= tw[d][e] == COS && g[e] == 0 ? 1 : 0.5;
                }
                energy += su * u[c][d][l] * u[c][d][l];
                Omega += sw * w[c][d][l] * w[c][d][l];
                Pal += sw * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) *
                       w[c][d][l] * w[c][d][l];
              }
            }
      sum[0] = energy;
      sum[1] = Omega;
      sum[2] = Pal;
      MPI_Allreduce(MPI_IN_PLACE, sum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      double S[9], Sb[9];
      for (int q = 0; q < 2; q++) {
        for (c = 0; c < 2; c++) {
          static const int tc[3] = {COS, COS, COS}, ts[3] = {SIN, COS, COS};
#pragma omp parallel for private(j, k, l)
          for (i = 0; i < N1; i++)
            for (j = 0; j < nl; j++)
              for (k = 0; k < N1; k++) {
                double m = wavenumber(c, i);
                l = yidx(i, j, k);
                du[c][0][l] = q == 0 ? m * u[c][0][l] : -m * m * u[c][0][l];
              }
          synthesis(c, q == 0 ? tc : ts, du[c][0], tmp, F[c]);
        }
        moments(F[EVEN], F[ODD], q == 0 ? S : Sb);
      }
      if (rank == 0) {
        printf("% 10ld % .16e % .16e % .16e", tstep, t, sum[0] / 2, sum[1] / 2);
        for (int n = 3; n < 9; n++)
          printf(" % .6e", S[n]);
        printf(" % .6e % .6e % .16e", Sb[4], Sb[6], sum[2] / 2);
        printf("\n");
        fflush(stdout);
      }
    }
    if (ne > 0 && tstep % ne == 0) {
      long nb = 2 * (long)(sqrt(3.0) * M) + 2;
      double *E = calloc(nb, sizeof(double));
      char path[FILENAME_MAX];
      FILE *file;
      for (c = 0; c < 2; c++)
        for (i = 0; i < N1; i++)
          for (j = 0; j < nl; j++)
            for (k = 0; k < N1; k++) {
              long g[3] = {wavenumber(c, i), wavenumber(c, s0 + j),
                           wavenumber(c, k)};
              l = yidx(i, j, k);
              long b = (long)(2 * sqrt((double)(g[0] * g[0] + g[1] * g[1] +
                                                g[2] * g[2])));
              for (d = 0; d < 3; d++) {
                double su = 1;
                for (int q = 0; q < 3; q++)
                  su *= tu[d][q] == COS && g[q] == 0 ? 1 : 0.5;
                E[b] += su * u[c][d][l] * u[c][d][l] / 2;
              }
            }
      MPI_Reduce(rank == 0 ? MPI_IN_PLACE : E, E, nb, MPI_DOUBLE, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      if (rank == 0) {
        sprintf(path, "e.%08ld", tstep);
        if ((file = fopen(path, "w")) == NULL) {
          fprintf(stderr, "tg: error: fail to open '%s'\n", path);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fprintf(file, "# t = %.16e\n", t);
        for (long b = 0; b < nb; b++)
          fprintf(file, "%.1f %.16e\n", b / 2.0, E[b]);
        if (fclose(file) != 0) {
          fprintf(stderr, "tg: error: fail to close '%s'\n", path);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
      free(E);
    }
    if (nd > 0 && tstep % nd == 0) {
      char path[FILENAME_MAX];
      FILE *file = NULL;
      double *g = NULL, *h = NULL;
      int *gc = NULL, *gd = NULL;
      if (rank == 0) {
        sprintf(path, "u.%08ld", tstep);
        if ((file = fopen(path, "w")) == NULL) {
          fprintf(stderr, "tg: error: fail to open '%s'\n", path);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        g = malloc(n3 * sizeof(double));
        h = malloc(n3 * sizeof(double));
        gc = malloc(size * sizeof(int));
        gd = malloc(size * sizeof(int));
        for (r = 0; r < size; r++) {
          gc[r] = cnt[r] * N1 * N1;
          gd[r] = off[r] * N1 * N1;
        }
      }
      for (c = 0; c < 2; c++)
        for (d = 0; d < 3; d++) {
          MPI_Gatherv(u[c][d], nloc, MPI_DOUBLE, g, gc, gd, MPI_DOUBLE, 0,
                      MPI_COMM_WORLD);
          if (rank == 0) {
            for (r = 0; r < size; r++)
              for (i = 0; i < N1; i++)
                for (j = 0; j < cnt[r]; j++)
                  memcpy(&h[(i * N1 + off[r] + j) * N1],
                         &g[gd[r] + (i * cnt[r] + j) * N1],
                         N1 * sizeof(double));
            if (fwrite(h, sizeof(double), n3, file) != (size_t)n3) {
              fprintf(stderr, "tg: error: fail to write '%s'\n", path);
              MPI_Abort(MPI_COMM_WORLD, 1);
            }
          }
        }
      if (rank == 0) {
        if (fclose(file) != 0) {
          fprintf(stderr, "tg: error: fail to close '%s'\n", path);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        free(g);
        free(h);
        free(gc);
        free(gd);
      }
    }
    if (t > T)
      break;
    for (c = 0; c < 2; c++)
      for (d = 0; d < 3; d++) {
        synthesis(c, tu[d], u[c][d], tmp, U[c][d]);
        synthesis(c, tw[d], w[c][d], tmp, W[c][d]);
      }
    for (c = 0; c < 2; c++) {
      double **ue = U[EVEN], **uo = U[ODD], **wa = W[c], **wb = W[1 - c];
#pragma omp parallel for
      for (l = 0; l < nloc; l++) {
        F[0][l] = ue[1][l] * wa[2][l] - ue[2][l] * wa[1][l] +
                  uo[1][l] * wb[2][l] - uo[2][l] * wb[1][l];
        F[1][l] = ue[2][l] * wa[0][l] - ue[0][l] * wa[2][l] +
                  uo[2][l] * wb[0][l] - uo[0][l] * wb[2][l];
        F[2][l] = ue[0][l] * wa[1][l] - ue[1][l] * wa[0][l] +
                  uo[0][l] * wb[1][l] - uo[1][l] * wb[0][l];
      }
      for (d = 0; d < 3; d++)
        analysis(c, tu[d], F[d], du[c][d]);
    }
    for (c = 0; c < 2; c++)
#pragma omp parallel for private(j, k, l, d)
      for (i = 0; i < N1; i++)
        for (j = 0; j < nl; j++)
          for (k = 0; k < N1; k++) {
            double m = wavenumber(c, i), n = wavenumber(c, s0 + j),
                   p = wavenumber(c, k), kk = m * m + n * n + p * p, P, cv,
                   v, f[3];
            l = yidx(i, j, k);
            for (d = 0; d < 3; d++)
              f[d] = du[c][d][l];
            if (i > kcut || s0 + j > kcut || k > kcut)
              f[0] = f[1] = f[2] = 0;
            P = kk > 0 ? (m * f[0] + n * f[1] + p * f[2]) / kk : 0;
            f[0] -= m * P;
            f[1] -= n * P;
            f[2] -= p * P;
            for (d = 0; d < 3; d++) {
              v = u[c][d][l];
              if (tstep == 0) {
                cv = nu * dt * kk / 2;
                u[c][d][l] = ((1 - cv) * v + dt * f[d]) / (1 + cv);
              } else {
                cv = nu * dt * kk;
                u[c][d][l] = ((1 - cv) * um[c][d][l] + 2 * dt * f[d]) / (1 + cv);
              }
              um[c][d][l] = v;
            }
          }
    t += dt;
    tstep++;
  }
  MPI_Finalize();
}
