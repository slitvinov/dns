#define _GNU_SOURCE
#include <assert.h>
#include <complex.h>
#include <fenv.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { N = 1 << 5, Nf = N / 2 + 1 };
static void backward(fftw_complex *U_hat, double *U, fftw_complex *work) {
  fftw_plan plan;
  memcpy(work, U_hat, N * N * Nf * sizeof(fftw_complex));
  plan = fftw_plan_dft_c2r_3d(N, N, N, work, U, FFTW_ESTIMATE);
  fftw_execute_dft_c2r(plan, work, U);
  fftw_destroy_plan(plan);
}

int main(void) {
  fftw_plan fplan;
  long double s;
  double dx, L, invN3;
  fftw_complex *curlX, *curlY, *curlZ, *dU, *dV, *dW, *P_hat, *U_hat, *U_hat0,
      *U_hat1, *V_hat, *V_hat0, *V_hat1, *W_hat, *W_hat0, *W_hat1, *work;
  int *dealias;
  long i, j, k, l, m, offset;
  size_t idump;
  int rk;
  int tstep;
  double a[] = {1 / 6.0, 1 / 3.0, 1 / 3.0, 1 / 6.0};
  double b[] = {0.5, 0.5, 1.0};
  double *CU, *CV, *CW, *kk, *kx, *kz, kmax, nu, dt, T, t, *U, *U_tmp, *V,
      *V_tmp, *W, *W_tmp;
  double pi = 3.141592653589793238;
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  nu = 0.000625;
  T = 0.1;
  dt = 0.01;
  L = 2 * pi;
  dx = L / N;
  invN3 = 1.0 / (N * N * N);
  U = fftw_alloc_real(N * N * N);
  V = fftw_alloc_real(N * N * N);
  W = fftw_alloc_real(N * N * N);
  U_tmp = fftw_alloc_real(N * N * N);
  V_tmp = fftw_alloc_real(N * N * N);
  W_tmp = fftw_alloc_real(N * N * N);
  CU = fftw_alloc_real(N * N * N);
  CV = fftw_alloc_real(N * N * N);
  CW = fftw_alloc_real(N * N * N);
  kx = malloc(N * sizeof(double));
  kz = malloc(Nf * sizeof(double));
  kk = malloc(N * N * Nf * sizeof(double));
  dealias = malloc(N * N * Nf * sizeof(int));
  U_hat = fftw_alloc_complex(N * N * Nf);
  V_hat = fftw_alloc_complex(N * N * Nf);
  W_hat = fftw_alloc_complex(N * N * Nf);
  P_hat = fftw_alloc_complex(N * N * Nf);
  U_hat0 = fftw_alloc_complex(N * N * Nf);
  V_hat0 = fftw_alloc_complex(N * N * Nf);
  W_hat0 = fftw_alloc_complex(N * N * Nf);
  U_hat1 = fftw_alloc_complex(N * N * Nf);
  V_hat1 = fftw_alloc_complex(N * N * Nf);
  W_hat1 = fftw_alloc_complex(N * N * Nf);
  dU = fftw_alloc_complex(N * N * Nf);
  dV = fftw_alloc_complex(N * N * Nf);
  dW = fftw_alloc_complex(N * N * Nf);
  curlX = fftw_alloc_complex(N * N * Nf);
  curlY = fftw_alloc_complex(N * N * Nf);
  curlZ = fftw_alloc_complex(N * N * Nf);
  work = fftw_alloc_complex(N * N * Nf);
  struct {
    double *var;
    const char *name;
  } list[] = {{U, "U"}, {V, "V"}, {W, "W"}};

  fplan = fftw_plan_dft_r2c_3d(N, N, N, U, U_hat,
                               FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

  for (i = 0; i < N / 2; i++) {
    kx[i] = i;
    kz[i] = i;
  }
  kz[N / 2] = N / 2;
  for (i = -N / 2; i < 0; i++)
    kx[i + N] = i;
  kmax = 2. / 3. * (N / 2 + 1);
  for (i = l = 00; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < Nf; k++, l++)
        dealias[l] = (fabs(kx[i]) < kmax) && (fabs(kx[j]) < kmax) &&
                     (fabs(kz[k]) < kmax);

  for (i = l = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < Nf; k++, l++) {
        m = kx[i] * kx[i] + kx[j] * kx[j] + kz[k] * kz[k];
        kk[l] = m > 0 ? m : 1;
      }

  for (i = l = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++, l++) {
        U[l] = sin(dx * i) * cos(dx * j) * cos(dx * k);
        V[l] = -cos(dx * i) * sin(dx * j) * cos(dx * k);
        W[l] = 0.0;
      }

  fftw_execute_dft_r2c(fplan, U, U_hat);
  fftw_execute_dft_r2c(fplan, V, V_hat);
  fftw_execute_dft_r2c(fplan, W, W_hat);

  t = 0.0;
  tstep = 0;
  while (t <= T) {
    backward(U_hat, U, work);
    backward(V_hat, V, work);
    backward(W_hat, W, work);

    for (k = 0; k < N * N * N; k++) {
      U[k] *= invN3;
      V[k] *= invN3;
      W[k] *= invN3;
    }
    if (tstep % 2 == 0) {
      s = 0.0;
      for (k = 0; k < N * N * N; k++)
        s += U[k] * U[k] + V[k] * V[k] + W[k] * W[k];
      s *= 0.5 * dx * dx * dx / L / L / L;
      fprintf(stderr, "eng = %.16Le\n", s);
      FILE *file;
      char path[FILENAME_MAX];
      sprintf(path, "%08d.raw", tstep);
      file = fopen(path, "w");
      for (idump = 0; idump < sizeof list / sizeof *list; idump++)
        fwrite(list[idump].var, N * N * N, sizeof(double), file);
      fclose(file);
      sprintf(path, "a.%08d.xdmf2", tstep);
      file = fopen(path, "w");
      fprintf(file,
              "<Xdmf\n"
              "    Version=\"2\">\n"
              "  <Domain>\n"
              "    <Grid>\n"
              "      <Topology\n"
              "	  TopologyType=\"3DCoRectMesh\"\n"
              "	  Dimensions=\"%d %d %d\"/>\n"
              "      <Geometry\n"
              "	  GeometryType=\"ORIGIN_DXDYDZ\">\n"
              "	<DataItem\n"
              "	    Dimensions=\"3\">\n"
              "	  0\n"
              "	  0\n"
              "	  0\n"
              "	</DataItem>\n"
              "	<DataItem\n"
              "	    Dimensions=\"3\">\n"
              "	  %.16e\n"
              "	  %.16e\n"
              "	  %.16e\n"
              "	</DataItem>\n"
              "      </Geometry>\n",
              N, N, N, dx, dx, dx);
      offset = 0;
      for (idump = 0; idump < sizeof list / sizeof *list; idump++) {
        fprintf(file,
                "      <Attribute\n"
                "          Name=\"%s\">\n"
                "        <DataItem\n"
                "            Format=\"Binary\"\n"
                "            Seek=\"%ld\"\n"
                "            Precision=\"8\"\n"
                "            Dimensions=\"%d %d %d\">\n"
                "          %08d.raw\n"
                "        </DataItem>\n"
                "      </Attribute>\n",
                list[idump].name, offset, N, N, N, tstep);
        offset += N * N * N * sizeof(double);
      }
      fprintf(file, "    </Grid>\n"
                    "  </Domain>\n"
                    "</Xdmf>\n");
      fclose(file);
    }
    memcpy(U_hat0, U_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(V_hat0, V_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(W_hat0, W_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(U_hat1, U_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(V_hat1, V_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(W_hat1, W_hat, sizeof(fftw_complex) * N * N * Nf);
    for (rk = 0; rk < 4; rk++) {
      if (rk > 0) {
        backward(U_hat, U, work);
        backward(V_hat, V, work);
        backward(W_hat, W, work);
        for (k = 0; k < N * N * N; k++) {
          U[k] *= invN3;
          V[k] *= invN3;
          W[k] *= invN3;
        }
      }
      for (i = l = 0; i < N; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < Nf; k++, l++) {
            curlZ[l] = I * (kx[i] * V_hat[l] - kx[j] * U_hat[l]);
            curlY[l] = I * (kz[k] * U_hat[l] - kx[i] * W_hat[l]);
            curlX[l] = I * (kx[j] * W_hat[l] - kz[k] * V_hat[l]);
          }
      backward(curlX, CU, work);
      backward(curlY, CV, work);
      backward(curlZ, CW, work);
      for (k = 0; k < N * N * N; k++) {
        CU[k] *= invN3;
        CV[k] *= invN3;
        CW[k] *= invN3;
      }
      for (k = 0; k < N * N * N; k++) {
        U_tmp[k] = V[k] * CW[k] - W[k] * CV[k];
        V_tmp[k] = W[k] * CU[k] - U[k] * CW[k];
        W_tmp[k] = U[k] * CV[k] - V[k] * CU[k];
      }
      fftw_execute_dft_r2c(fplan, U_tmp, dU);
      fftw_execute_dft_r2c(fplan, V_tmp, dV);
      fftw_execute_dft_r2c(fplan, W_tmp, dW);

      for (k = 0; k < N * N * Nf; k++) {
        dU[k] *= dealias[k] * dt;
        dV[k] *= dealias[k] * dt;
        dW[k] *= dealias[k] * dt;
      }
      for (i = l = 0; i < N; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < Nf; k++, l++) {
            P_hat[l] = (dU[l] * kx[i] + dV[l] * kx[j] + dW[l] * kz[k]) / kk[l];
            dU[l] -= P_hat[l] * kx[i] + nu * dt * kk[l] * U_hat[l];
            dV[l] -= P_hat[l] * kx[j] + nu * dt * kk[l] * V_hat[l];
            dW[l] -= P_hat[l] * kz[k] + nu * dt * kk[l] * W_hat[l];
          }

      if (rk < 3) {
        for (k = 0; k < N * N * Nf; k++) {
          U_hat[k] = U_hat0[k] + b[rk] * dU[k];
          V_hat[k] = V_hat0[k] + b[rk] * dV[k];
          W_hat[k] = W_hat0[k] + b[rk] * dW[k];
        }
      }
      for (k = 0; k < N * N * Nf; ++k) {
        U_hat1[k] += a[rk] * dU[k];
        V_hat1[k] += a[rk] * dV[k];
        W_hat1[k] += a[rk] * dW[k];
      }
    }
    memcpy(U_hat, U_hat1, sizeof(fftw_complex) * N * N * Nf);
    memcpy(V_hat, V_hat1, sizeof(fftw_complex) * N * N * Nf);
    memcpy(W_hat, W_hat1, sizeof(fftw_complex) * N * N * Nf);
    t += dt;
    tstep++;
  }

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
  fftw_free(V_tmp);
  fftw_free(W);
  fftw_free(W_hat);
  fftw_free(W_tmp);
  fftw_free(work);

  free(dealias);
  free(kx);
  free(kz);
  free(kk);
}
