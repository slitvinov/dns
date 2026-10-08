#define _GNU_SOURCE
#include <complex.h>
#include <fenv.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { nvars = 4 };
static const double pi = 3.141592653589793238;
static const double a[] = {1 / 6.0, 1 / 3.0, 1 / 3.0, 1 / 6.0};
static const double b[] = {0.5, 0.5, 1.0};
static int rank;
static void parallel_loop(void *(*work)(char *), char *jobdata, size_t elsize,
                          int njobs, void *data) {
  (void)data;
#pragma omp parallel for
  for (int i = 0; i < njobs; ++i)
    work(jobdata + elsize * i);
}
static void fail(const char *msg, const char *path) {
  fprintf(stderr, "fourier: error: %s '%s'\n", msg, path);
  MPI_Abort(MPI_COMM_WORLD, 1);
}
static void c2r(fftw_plan plan, long m, fftw_complex *hat, double *real,
                fftw_complex *work) {
#pragma omp parallel for
  for (long i = 0; i < m; i++)
    work[i] = hat[i];
  fftw_mpi_execute_dft_c2r(plan, work, real);
}
static double cabs2(fftw_complex z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}
static long steps(const char *opt, double x, double dt) {
  long n = lround(x / dt);
  if (fabs(n * dt - x) > 1e-9 * x) {
    fprintf(stderr, "fourier: error: %s %g is not a multiple of -s %g\n", opt, x,
            dt);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return n;
}
int main(int argc, char **argv) {
  fftw_plan fplan, bplan;
  MPI_File fh;
  MPI_Offset size;
  FILE *file;
  char path[FILENAME_MAX], *input_path, *end;
  double energy, Omega, sum[2], dx, L, invn3, kmax, nu, dt, T, t;
  fftw_complex *curlX, *curlY, *curlZ, *dU, *dV, *dW, *P_hat, *U_hat, *U_hat0,
      *U_hat1, *V_hat, *V_hat0, *V_hat1, *W_hat, *W_hat0, *W_hat1, *dump_hat;
  int *dealias, rk, Verbose, Dump, provided, nproc;
  long idump, tstep, nt;
  ptrdiff_t alloc, n0, s0, n1, s1;
  size_t ivar;
  double *CU, *CV, *CW, *kk, *kx, *kz, *U, *U_tmp, *V, *V_tmp, *W, *W_tmp,
      *dump, *row;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  input_path = NULL;
  dt = -1;
  T = 0;
  nu = -1;
  Verbose = 0;
  Dump = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    switch (argv[0][1]) {
    case 'h':
      if (rank == 0)
        fprintf(stderr,
                "Usage: fourier [-v] [-d] -i <input.raw> -n <viscosity> -t "
                "<end time> -s <time step>\n");
      MPI_Finalize();
      exit(1);
    case 'v':
      Verbose = 1;
      break;
    case 'd':
      Dump = 1;
      break;
    case 'i':
    case 'n':
    case 's':
    case 't':
      if (argv[1] == NULL)
        fail("needs an argument", argv[0]);
      if (argv[0][1] == 'i') {
        input_path = argv[1];
      } else {
        double x = strtod(argv[1], &end);
        if (*end != '\0')
          fail("not a number", argv[1]);
        if (argv[0][1] == 'n')
          nu = x;
        else if (argv[0][1] == 's')
          dt = x;
        else
          T = x;
      }
      argv++;
      break;
    default:
      fail("unknown option", *argv);
    }
  }
  if (T == 0 || nu == -1 || dt == -1 || input_path == NULL)
    fail("need -i, -n, -s, -t", "");
  fftw_init_threads();
  fftw_mpi_init();
  fftw_plan_with_nthreads(omp_get_max_threads());
  fftw_threads_set_callback(parallel_loop, NULL);
  if (MPI_File_open(MPI_COMM_WORLD, input_path, MPI_MODE_RDONLY, MPI_INFO_NULL,
                    &fh) != MPI_SUCCESS)
    fail("fail to open", input_path);
  MPI_File_get_size(fh, &size);
  long n = size / sizeof(double) / nvars;
  n = lround(cbrt((double)n));
  if (n * n * n * nvars * (MPI_Offset)sizeof(double) != size)
    fail("wrong file", input_path);
  long nf = n / 2 + 1;
  long n3 = n * n * n;
  alloc = fftw_mpi_local_size_3d_transposed(n, n, nf, MPI_COMM_WORLD, &n0, &s0,
                                            &n1, &s1);
  if (Verbose && rank == 0)
    fprintf(stderr, "fourier: n = %ld, ranks %d, threads %d\n", n, nproc,
            omp_get_max_threads());
  long m = n1 * n * nf;
  U = fftw_alloc_real(2 * alloc);
  V = fftw_alloc_real(2 * alloc);
  W = fftw_alloc_real(2 * alloc);
  row = fftw_alloc_real(n0 * n * n + 1);
  double *var[3] = {U, V, W};
  for (int d = 0; d < 3; d++) {
    if (MPI_File_read_at_all(fh, (d * n3 + s0 * n * n) * sizeof(double), row,
                             n0 * n * n, MPI_DOUBLE,
                             MPI_STATUS_IGNORE) != MPI_SUCCESS)
      fail("fail to read", input_path);
    for (long i = 0; i < n0; i++)
      for (long j = 0; j < n; j++)
        memcpy(&var[d][(i * n + j) * 2 * nf], &row[(i * n + j) * n],
               n * sizeof(double));
  }
  if (MPI_File_close(&fh) != MPI_SUCCESS)
    fail("fail to close", input_path);
  L = 2 * pi;
  dx = L / n;
  invn3 = 1.0 / n3;
  dump = fftw_alloc_real(2 * alloc);
  U_tmp = fftw_alloc_real(2 * alloc);
  V_tmp = fftw_alloc_real(2 * alloc);
  W_tmp = fftw_alloc_real(2 * alloc);
  CU = fftw_alloc_real(2 * alloc);
  CV = fftw_alloc_real(2 * alloc);
  CW = fftw_alloc_real(2 * alloc);
  kx = malloc(n * sizeof(double));
  kz = malloc(nf * sizeof(double));
  kk = malloc(m * sizeof(double));
  dealias = malloc(m * sizeof(int));
  U_hat = fftw_alloc_complex(alloc);
  V_hat = fftw_alloc_complex(alloc);
  W_hat = fftw_alloc_complex(alloc);
  P_hat = fftw_alloc_complex(alloc);
  U_hat0 = fftw_alloc_complex(alloc);
  V_hat0 = fftw_alloc_complex(alloc);
  W_hat0 = fftw_alloc_complex(alloc);
  U_hat1 = fftw_alloc_complex(alloc);
  V_hat1 = fftw_alloc_complex(alloc);
  W_hat1 = fftw_alloc_complex(alloc);
  dU = fftw_alloc_complex(alloc);
  dV = fftw_alloc_complex(alloc);
  dW = fftw_alloc_complex(alloc);
  curlX = fftw_alloc_complex(alloc);
  curlY = fftw_alloc_complex(alloc);
  curlZ = fftw_alloc_complex(alloc);
  dump_hat = fftw_alloc_complex(alloc);
  struct {
    fftw_complex *var;
    const char *name;
  } list[nvars] = {{U_hat, "U"}, {V_hat, "V"}, {W_hat, "W"}, {P_hat, "P"}};
  memset(P_hat, 0, alloc * sizeof(fftw_complex));
  fplan = fftw_mpi_plan_dft_r2c_3d(n, n, n, CU, dU, MPI_COMM_WORLD,
                                   FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
  bplan = fftw_mpi_plan_dft_c2r_3d(n, n, n, dU, CU, MPI_COMM_WORLD,
                                   FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);
  if (fplan == NULL || bplan == NULL)
    fail("fftw_mpi_plan failed", "");
  for (long i = 0; i < n / 2; i++) {
    kx[i] = i;
    kz[i] = i;
  }
  kz[n / 2] = n / 2;
  for (long i = -n / 2; i < 0; i++)
    kx[i + n] = i;
  kmax = 2. / 3. * (n / 2 + 1);
#pragma omp parallel for collapse(3)
  for (long j = 0; j < n1; j++)
    for (long i = 0; i < n; i++)
      for (long k = 0; k < nf; k++) {
        long l = (j * n + i) * nf + k;
        double y = kx[j + s1], x = kx[i];
        dealias[l] = fabs(x) < kmax && fabs(y) < kmax && fabs(kz[k]) < kmax;
        kk[l] = x * x + y * y + kz[k] * kz[k];
      }
  fftw_mpi_execute_dft_r2c(fplan, U, U_hat);
  fftw_mpi_execute_dft_r2c(fplan, V, V_hat);
  fftw_mpi_execute_dft_r2c(fplan, W, W_hat);
  idump = 0;
  t = 0.0;
  nt = steps("-t", T, dt);
  tstep = 0;
  for (;;) {
    if (tstep % 10 == 0) {
      energy = 0.0;
      Omega = 0.0;
#pragma omp parallel for reduction(+ : energy, Omega)
      for (long l = 0; l < m; l++) {
        double h = l % nf == 0 || l % nf == n / 2 ? 0.5 : 1;
        double e = cabs2(U_hat[l]) + cabs2(V_hat[l]) + cabs2(W_hat[l]);
        energy += h * e;
        Omega += h * kk[l] * e;
      }
      sum[0] = energy * invn3 * invn3;
      sum[1] = Omega * invn3 * invn3;
      MPI_Allreduce(MPI_IN_PLACE, sum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (rank == 0) {
        printf("% 10ld % .16e % .16e % .16e\n", tstep, t, sum[0], sum[1]);
        if (fflush(stdout) != 0 || ferror(stdout))
          fail("fail to write", "stdout");
      }
      if (Dump) {
        sprintf(path, "%08ld.raw", tstep);
        if (MPI_File_open(MPI_COMM_WORLD, path,
                          MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
                          &fh) != MPI_SUCCESS)
          fail("fail to open", path);
        for (ivar = 0; ivar < sizeof list / sizeof *list; ivar++) {
          c2r(bplan, m, list[ivar].var, dump, dump_hat);
          for (long i = 0; i < n0; i++)
            for (long j = 0; j < n; j++)
              for (long k = 0; k < n; k++)
                row[(i * n + j) * n + k] = dump[(i * n + j) * 2 * nf + k] * invn3;
          if (MPI_File_write_at_all(fh, (ivar * n3 + s0 * n * n) * sizeof(double),
                                    row, n0 * n * n, MPI_DOUBLE,
                                    MPI_STATUS_IGNORE) != MPI_SUCCESS)
            fail("fail to write", path);
        }
        if (MPI_File_close(&fh) != MPI_SUCCESS)
          fail("fail to close", path);
        if (rank == 0) {
          sprintf(path, "a.%08ld.xdmf2", idump);
          if ((file = fopen(path, "w")) == NULL)
            fail("fail to open", path);
          fprintf(file,
                  "<Xdmf\n"
                  "    Version=\"2\">\n"
                  "  <Domain>\n"
                  "    <Grid>\n"
                  "      <Time\n"
                  "          Value=\"%+.16e\"/>\n"
                  "      <Topology\n"
                  "          TopologyType=\"3DCoRectMesh\"\n"
                  "          Dimensions=\"%ld %ld %ld\"/>\n"
                  "      <Geometry\n"
                  "          GeometryType=\"ORIGIN_DXDYDZ\">\n"
                  "        <DataItem\n"
                  "            Dimensions=\"3\">\n"
                  "          0\n"
                  "          0\n"
                  "          0\n"
                  "        </DataItem>\n"
                  "        <DataItem\n"
                  "            Dimensions=\"3\">\n"
                  "          %.16e\n"
                  "          %.16e\n"
                  "          %.16e\n"
                  "        </DataItem>\n"
                  "      </Geometry>\n",
                  t, n, n, n, dx, dx, dx);
          for (ivar = 0; ivar < sizeof list / sizeof *list; ivar++)
            fprintf(file,
                    "      <Attribute\n"
                    "          name=\"%s\">\n"
                    "        <DataItem\n"
                    "            Format=\"Binary\"\n"
                    "            Seek=\"%ld\"\n"
                    "            Precision=\"8\"\n"
                    "            Dimensions=\"%ld %ld %ld\">\n"
                    "          %08ld.raw\n"
                    "        </DataItem>\n"
                    "      </Attribute>\n",
                    list[ivar].name, (long)(ivar * n3 * sizeof(double)), n, n,
                    n, tstep);
          fprintf(file, "    </Grid>\n"
                        "  </Domain>\n"
                        "</Xdmf>\n");
          if (ferror(file) || fclose(file) != 0)
            fail("fail to write", path);
        }
        idump++;
      }
    }
    if (tstep >= nt)
      break;
#pragma omp parallel for
    for (long l = 0; l < m; l++) {
      U_hat0[l] = U_hat1[l] = U_hat[l];
      V_hat0[l] = V_hat1[l] = V_hat[l];
      W_hat0[l] = W_hat1[l] = W_hat[l];
    }
    for (rk = 0; rk < 4; rk++) {
      c2r(bplan, m, U_hat, U, curlX);
      c2r(bplan, m, V_hat, V, curlX);
      c2r(bplan, m, W_hat, W, curlX);
#pragma omp parallel for collapse(3)
      for (long j = 0; j < n1; j++)
        for (long i = 0; i < n; i++)
          for (long k = 0; k < nf; k++) {
            long l = (j * n + i) * nf + k;
            double x = kx[i], y = kx[j + s1], z = kz[k];
            curlZ[l] = I * (x * V_hat[l] - y * U_hat[l]);
            curlY[l] = I * (z * U_hat[l] - x * W_hat[l]);
            curlX[l] = I * (y * W_hat[l] - z * V_hat[l]);
          }
      fftw_mpi_execute_dft_c2r(bplan, curlX, CU);
      fftw_mpi_execute_dft_c2r(bplan, curlY, CV);
      fftw_mpi_execute_dft_c2r(bplan, curlZ, CW);
#pragma omp parallel for collapse(2)
      for (long i = 0; i < n0; i++)
        for (long j = 0; j < n; j++)
          for (long k = 0; k < n; k++) {
            long l = (i * n + j) * 2 * nf + k;
            double u = U[l] * invn3, v = V[l] * invn3, w = W[l] * invn3,
                   cu = CU[l] * invn3, cv = CV[l] * invn3, cw = CW[l] * invn3;
            U_tmp[l] = v * cw - w * cv;
            V_tmp[l] = w * cu - u * cw;
            W_tmp[l] = u * cv - v * cu;
          }
      fftw_mpi_execute_dft_r2c(fplan, U_tmp, dU);
      fftw_mpi_execute_dft_r2c(fplan, V_tmp, dV);
      fftw_mpi_execute_dft_r2c(fplan, W_tmp, dW);
#pragma omp parallel for collapse(3)
      for (long j = 0; j < n1; j++)
        for (long i = 0; i < n; i++)
          for (long k = 0; k < nf; k++) {
            long l = (j * n + i) * nf + k;
            double x = kx[i], y = kx[j + s1], z = kz[k];
            dU[l] *= dealias[l] * dt;
            dV[l] *= dealias[l] * dt;
            dW[l] *= dealias[l] * dt;
            P_hat[l] =
                kk[l] > 0 ? (dU[l] * x + dV[l] * y + dW[l] * z) / kk[l] : 0.0;
            dU[l] -= P_hat[l] * x + nu * dt * kk[l] * U_hat[l];
            dV[l] -= P_hat[l] * y + nu * dt * kk[l] * V_hat[l];
            dW[l] -= P_hat[l] * z + nu * dt * kk[l] * W_hat[l];
          }
      if (rk < 3) {
#pragma omp parallel for
        for (long l = 0; l < m; l++) {
          U_hat[l] = U_hat0[l] + b[rk] * dU[l];
          V_hat[l] = V_hat0[l] + b[rk] * dV[l];
          W_hat[l] = W_hat0[l] + b[rk] * dW[l];
        }
      }
#pragma omp parallel for
      for (long l = 0; l < m; ++l) {
        U_hat1[l] += a[rk] * dU[l];
        V_hat1[l] += a[rk] * dV[l];
        W_hat1[l] += a[rk] * dW[l];
      }
    }
#pragma omp parallel for
    for (long l = 0; l < m; l++) {
      U_hat[l] = U_hat1[l];
      V_hat[l] = V_hat1[l];
      W_hat[l] = W_hat1[l];
    }
    tstep++;
    t = tstep * dt;
  }
  fftw_destroy_plan(fplan);
  fftw_destroy_plan(bplan);
  fftw_mpi_cleanup();
  fftw_cleanup_threads();
  MPI_Finalize();
}
