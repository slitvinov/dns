#define _GNU_SOURCE
#include <assert.h>
#include <complex.h>
#include <fenv.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { nvars = 4 };
static const double pi = 3.141592653589793238;
static const double a[] = {1 / 6.0, 1 / 3.0, 1 / 3.0, 1 / 6.0};
static const double b[] = {0.5, 0.5, 1.0};
static void c2r(fftw_plan fplan, long n3f, fftw_complex *hat, double *real,
                fftw_complex *work) {
  memcpy(work, hat, n3f * sizeof(fftw_complex));
  fftw_execute_dft_c2r(fplan, work, real);
}

int main(int argc, char **argv) {
  (void)argc;
  fftw_plan fplan, bplan;
  FILE *file;
  char path[FILENAME_MAX], *input_path, *end;
  long double s;
  double dx, L, invn3, k2, kmax, nu, dt, T, t;
  fftw_complex *curlX, *curlY, *curlZ, *dU, *dV, *dW, *P_hat, *U_hat, *U_hat0,
      *U_hat1, *V_hat, *V_hat0, *V_hat1, *W_hat, *W_hat0, *W_hat1, *dump_hat;
  int *dealias, rk, tstep, Verbose;
  long i, j, k, l, idump;
  size_t offset;
  size_t ivar;
  double *CU, *CV, *CW, *kk, *kx, *kz, *U, *U_tmp, *V, *V_tmp, *W, *W_tmp,
      *dump;
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  input_path = NULL;
  T = 0;
  nu = -1;
  Verbose = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr,
              "Usage: dns [-v] -i input.raw -n viscosity -t <end time>\n\n"
              "Example:\n"
              "  dns -i input.raw -m 0.1 -t 1.0\n\n");
      break;
    case 'v':
      Verbose = 1;
      break;
    case 'i':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "dns: error: -i needs an argument\n");
        exit(1);
      }
      input_path = *argv;
      break;
    case 'n':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "dns: error: -n needs an argument\n");
        exit(1);
      }
      nu = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(stderr, "dns: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      break;
    case 't':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "dns: error: -i needs an argument\n");
        exit(1);
      }
      T = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(stderr, "dns: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      break;
    default:
      fprintf(stderr, "dns: error: unknown option '%s'\n", *argv);
      exit(1);
    }
  }
  if (T == 0) {
    fprintf(stderr, "dns: error: -t is not set or invalid\n");
    exit(1);
  }
  if (nu == -1) {
    fprintf(stderr, "dns: error: -n is not set or invalid\n");
    exit(1);
  }
  if (input_path == NULL) {
    fprintf(stderr, "dns: error: -i is not set\n");
    exit(1);
  }
  if ((file = fopen(input_path, "r")) == NULL) {
    fprintf(stderr, "dns: error: fail to open '%s'\n", input_path);
    exit(1);
  }
  fseek(file, 0, SEEK_END);
  offset = ftell(file);
  rewind(file);
  long n = offset / sizeof(double) / nvars;
  n = round(powf(n, 1.0 / 3));
  if (n * n * n * nvars * sizeof(double) != offset) {
    fprintf(stderr, "dns: error: wrong file '%s'\n", input_path);
    exit(1);
  }
  if (Verbose)
    fprintf(stderr, "dns: n = %ld\n", n);
  long nf = n / 2 + 1;
  long n3 = n * n * n;
  long n3f = n * n * nf;
  U = fftw_alloc_real(n3);
  V = fftw_alloc_real(n3);
  W = fftw_alloc_real(n3);
  if (fread(U, sizeof(double), n3, file) != (size_t)n3 ||
      fread(V, sizeof(double), n3, file) != (size_t)n3 ||
      fread(W, sizeof(double), n3, file) != (size_t)n3 || fclose(file) != 0) {
    fprintf(stderr, "dns: error: fail to read '%s'\n", input_path);
    exit(1);
  }
  dt = 0.01;
  L = 2 * pi;
  dx = L / n;
  invn3 = 1.0 / n3;
  dump = fftw_alloc_real(n3);
  U_tmp = fftw_alloc_real(n3);
  V_tmp = fftw_alloc_real(n3);
  W_tmp = fftw_alloc_real(n3);
  CU = fftw_alloc_real(n3);
  CV = fftw_alloc_real(n3);
  CW = fftw_alloc_real(n3);
  kx = malloc(n * sizeof(double));
  kz = malloc(nf * sizeof(double));
  kk = malloc(n3f * sizeof(double));
  dealias = malloc(n3f * sizeof(int));
  U_hat = fftw_alloc_complex(n3f);
  V_hat = fftw_alloc_complex(n3f);
  W_hat = fftw_alloc_complex(n3f);
  P_hat = fftw_alloc_complex(n3f);
  U_hat0 = fftw_alloc_complex(n3f);
  V_hat0 = fftw_alloc_complex(n3f);
  W_hat0 = fftw_alloc_complex(n3f);
  U_hat1 = fftw_alloc_complex(n3f);
  V_hat1 = fftw_alloc_complex(n3f);
  W_hat1 = fftw_alloc_complex(n3f);
  dU = fftw_alloc_complex(n3f);
  dV = fftw_alloc_complex(n3f);
  dW = fftw_alloc_complex(n3f);
  curlX = fftw_alloc_complex(n3f);
  curlY = fftw_alloc_complex(n3f);
  curlZ = fftw_alloc_complex(n3f);
  dump_hat = fftw_alloc_complex(n3f);
  struct {
    fftw_complex *var;
    const char *name;
  } list[nvars] = {{U_hat, "U"}, {V_hat, "V"}, {W_hat, "W"}, {P_hat, "P"}};
  fplan = fftw_plan_dft_r2c_3d(n, n, n, U, U_hat,
                               FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  bplan = fftw_plan_dft_c2r_3d(n, n, n, U_hat, U, FFTW_ESTIMATE);
  for (i = 0; i < n / 2; i++) {
    kx[i] = i;
    kz[i] = i;
  }
  kz[n / 2] = n / 2;
  for (i = -n / 2; i < 0; i++)
    kx[i + n] = i;
  kmax = 2. / 3. * (n / 2 + 1);
  for (i = l = 00; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < nf; k++, l++)
        dealias[l] = (fabs(kx[i]) < kmax) && (fabs(kx[j]) < kmax) &&
                     (fabs(kz[k]) < kmax);

  for (i = l = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < nf; k++, l++) {
        k2 = kx[i] * kx[i] + kx[j] * kx[j] + kz[k] * kz[k];
        kk[l] = k2 > 0 ? k2 : 1;
      }

  fftw_execute_dft_r2c(fplan, U, U_hat);
  fftw_execute_dft_r2c(fplan, V, V_hat);
  fftw_execute_dft_r2c(fplan, W, W_hat);

  idump = 0;
  t = 0.0;
  tstep = 0;
  for (;;) {
    if (tstep % 10 == 0) {
      s = 0.0;
      for (k = 0; k < n3; k++)
        s += U[k] * U[k] + V[k] * V[k] + W[k] * W[k];
      s *= 0.5 * dx * dx * dx / L / L / L;
      fprintf(stderr, "dns: % 8ld % .4e % .4Le\n", idump, t, s);
      sprintf(path, "%08ld.raw", idump);
      file = fopen(path, "w");
      for (ivar = 0; ivar < sizeof list / sizeof *list; ivar++) {
        memcpy(dump_hat, list[ivar].var, n3f * sizeof(fftw_complex));
        fftw_execute_dft_c2r(bplan, dump_hat, dump);
        for (i = 0; i < n3; i++)
          dump[i] *= invn3;
        fwrite(dump, n3, sizeof(double), file);
      }
      fclose(file);
      sprintf(path, "a.%08ld.xdmf2", idump);
      file = fopen(path, "w");
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
              "          GeometryType=\"ORIGIn_DXDYDZ\">\n"
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
      offset = 0;
      for (ivar = 0; ivar < sizeof list / sizeof *list; ivar++) {
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
                list[ivar].name, offset, n, n, n, idump);
        offset += n3 * sizeof(double);
      }
      fprintf(file, "    </Grid>\n"
                    "  </Domain>\n"
                    "</Xdmf>\n");
      fclose(file);
      idump++;
    }
    if (t > T)
      break;
    memcpy(U_hat0, U_hat, sizeof(fftw_complex) * n3f);
    memcpy(V_hat0, V_hat, sizeof(fftw_complex) * n3f);
    memcpy(W_hat0, W_hat, sizeof(fftw_complex) * n3f);
    memcpy(U_hat1, U_hat, sizeof(fftw_complex) * n3f);
    memcpy(V_hat1, V_hat, sizeof(fftw_complex) * n3f);
    memcpy(W_hat1, W_hat, sizeof(fftw_complex) * n3f);
    for (rk = 0; rk < 4; rk++) {
      if (rk > 0) {
        c2r(bplan, n3f, U_hat, U, curlX);
        c2r(bplan, n3f, V_hat, V, curlX);
        c2r(bplan, n3f, W_hat, W, curlX);
        for (k = 0; k < n3; k++) {
          U[k] *= invn3;
          V[k] *= invn3;
          W[k] *= invn3;
        }
      }
      for (i = l = 0; i < n; i++)
        for (j = 0; j < n; j++)
          for (k = 0; k < nf; k++, l++) {
            curlZ[l] = I * (kx[i] * V_hat[l] - kx[j] * U_hat[l]);
            curlY[l] = I * (kz[k] * U_hat[l] - kx[i] * W_hat[l]);
            curlX[l] = I * (kx[j] * W_hat[l] - kz[k] * V_hat[l]);
          }
      fftw_execute_dft_c2r(bplan, curlX, CU);
      fftw_execute_dft_c2r(bplan, curlY, CV);
      fftw_execute_dft_c2r(bplan, curlZ, CW);
      for (k = 0; k < n3; k++) {
        CU[k] *= invn3;
        CV[k] *= invn3;
        CW[k] *= invn3;
      }
      for (k = 0; k < n3; k++) {
        U_tmp[k] = V[k] * CW[k] - W[k] * CV[k];
        V_tmp[k] = W[k] * CU[k] - U[k] * CW[k];
        W_tmp[k] = U[k] * CV[k] - V[k] * CU[k];
      }
      fftw_execute_dft_r2c(fplan, U_tmp, dU);
      fftw_execute_dft_r2c(fplan, V_tmp, dV);
      fftw_execute_dft_r2c(fplan, W_tmp, dW);

      for (k = 0; k < n3f; k++) {
        dU[k] *= dealias[k] * dt;
        dV[k] *= dealias[k] * dt;
        dW[k] *= dealias[k] * dt;
      }
      for (i = l = 0; i < n; i++)
        for (j = 0; j < n; j++)
          for (k = 0; k < nf; k++, l++) {
            P_hat[l] = (dU[l] * kx[i] + dV[l] * kx[j] + dW[l] * kz[k]) / kk[l];
            dU[l] -= P_hat[l] * kx[i] + nu * dt * kk[l] * U_hat[l];
            dV[l] -= P_hat[l] * kx[j] + nu * dt * kk[l] * V_hat[l];
            dW[l] -= P_hat[l] * kz[k] + nu * dt * kk[l] * W_hat[l];
          }

      if (rk < 3) {
        for (k = 0; k < n3f; k++) {
          U_hat[k] = U_hat0[k] + b[rk] * dU[k];
          V_hat[k] = V_hat0[k] + b[rk] * dV[k];
          W_hat[k] = W_hat0[k] + b[rk] * dW[k];
        }
      }
      for (k = 0; k < n3f; ++k) {
        U_hat1[k] += a[rk] * dU[k];
        V_hat1[k] += a[rk] * dV[k];
        W_hat1[k] += a[rk] * dW[k];
      }
    }
    memcpy(U_hat, U_hat1, sizeof(fftw_complex) * n3f);
    memcpy(V_hat, V_hat1, sizeof(fftw_complex) * n3f);
    memcpy(W_hat, W_hat1, sizeof(fftw_complex) * n3f);
    t += dt;
    tstep++;
  }

  fftw_destroy_plan(fplan);
  fftw_destroy_plan(bplan);
  fftw_free(CU);
  fftw_free(curlX);
  fftw_free(curlY);
  fftw_free(curlZ);
  fftw_free(CV);
  fftw_free(CW);
  fftw_free(dU);
  fftw_free(dV);
  fftw_free(dW);
  fftw_free(P_hat);
  fftw_free(U);
  fftw_free(U_hat);
  fftw_free(U_tmp);
  fftw_free(V);
  fftw_free(V_hat);
  fftw_free(V_hat1);
  fftw_free(U_hat1);
  fftw_free(W_hat1);
  fftw_free(V_hat0);
  fftw_free(U_hat0);
  fftw_free(W_hat0);
  fftw_free(V_tmp);
  fftw_free(W);
  fftw_free(W_hat);
  fftw_free(W_tmp);
  fftw_free(dump_hat);

  free(dump);
  free(dealias);
  free(kx);
  free(kz);
  free(kk);
}
