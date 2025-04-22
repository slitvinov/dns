#define _GNU_SOURCE
#include <assert.h>
#include <complex.h>
#include <fenv.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { N = 1 << 5, Nf = N / 2 + 1, tot = N * N * N };
#define MALLOC(var, nelem)                                                     \
  if ((var = fftw_malloc((nelem) * sizeof *var)) == NULL) {                    \
    fprintf(stderr, "%s:%d: fftwf_malloc failed\n", __FILE__, __LINE__);       \
    exit(1);                                                                   \
  }

static void forward(double *U, fftw_complex *U_hat) {
  fftw_plan plan;
  plan = fftw_plan_dft_r2c_3d(N, N, N, U, U_hat, FFTW_ESTIMATE);
  fftw_execute_dft_r2c(plan, U, U_hat);
  fftw_destroy_plan(plan);
}
static void backward(fftw_complex *U_hat, double *U) {
  fftw_plan plan;
  plan = fftw_plan_dft_c2r_3d(N, N, N, U_hat, U, FFTW_ESTIMATE);
  fftw_execute_dft_c2r(plan, U_hat, U);
  fftw_destroy_plan(plan);
}

int main(void) {
  double dx, L, s;
  fftw_complex *curlX, *curlY, *curlZ, *dU, *dV, *dW, *P_hat, *U_hat, *U_hat0,
      *U_hat1, *V_hat, *V_hat0, *V_hat1, *W_hat, *W_hat0, *W_hat1;
  int *dealias;
  long i, j, k, m;
  int rk;
  int tstep;
  int z;
  double a[] = {1 / 6.0, 1 / 3.0, 1 / 3.0, 1 / 6.0};
  double b[] = {0.5, 0.5, 1.0};
  double *CU;
  double *CV;
  double *CW;
  double *kk, *kx, *kz, kmax, nu, dt, T, t;
  double *U;
  double *U_tmp;
  double *V;
  double *V_tmp;
  double *W;
  double *W_tmp;
  double pi = 3.141592653589793238;

  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  nu = 0.000625;
  T = 0.1;
  dt = 0.01;
  L = 2 * pi;
  dx = L / N;
  U = fftw_alloc_real(2 * N * N * Nf);
  V = fftw_alloc_real(2 * N * N * Nf);
  W = fftw_alloc_real(2 * N * N * Nf);
  U_tmp = fftw_alloc_real(2 * N * N * Nf);
  V_tmp = fftw_alloc_real(2 * N * N * Nf);
  W_tmp = fftw_alloc_real(2 * N * N * Nf);
  CU = fftw_alloc_real(2 * N * N * Nf);
  CV = fftw_alloc_real(2 * N * N * Nf);
  CW = fftw_alloc_real(2 * N * N * Nf);

  kx = malloc(N * sizeof(double));
  kz = malloc(N * sizeof(double));
  kk = malloc(N * N * Nf * sizeof(double));

  dealias = malloc(N * N * Nf * sizeof(int));
  MALLOC(U_hat, N * N * Nf);
  MALLOC(V_hat, N * N * Nf);
  MALLOC(W_hat, N * N * Nf);
  MALLOC(P_hat, N * N * Nf);
  MALLOC(U_hat0, N * N * Nf);
  MALLOC(V_hat0, N * N * Nf);
  MALLOC(W_hat0, N * N * Nf);
  MALLOC(U_hat1, N * N * Nf);
  MALLOC(V_hat1, N * N * Nf);
  MALLOC(W_hat1, N * N * Nf);
  MALLOC(dU, N * N * Nf);
  MALLOC(dV, N * N * Nf);
  MALLOC(dW, N * N * Nf);
  MALLOC(curlX, N * N * Nf);
  MALLOC(curlY, N * N * Nf);
  MALLOC(curlZ, N * N * Nf);
  for (i = 0; i < N / 2; i++) {
    kx[i] = i;
    kz[i] = i;
  }
  kz[N / 2] = N / 2;
  for (i = -N / 2; i < 0; i++)
    kx[i + N] = i;
  memset(U, 0, 2 * N * N * Nf * sizeof(double));
  memset(V, 0, 2 * N * N * Nf * sizeof(double));
  memset(W, 0, 2 * N * N * Nf * sizeof(double));
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++) {
        z = (i * N + j) * N + k;
        U[z] = sin(dx * i) * cos(dx * j) * cos(dx * k);
        V[z] = -cos(dx * i) * sin(dx * j) * cos(dx * k);
        W[z] = 0.0;
      }

  forward(U, U_hat);
  forward(V, V_hat);
  forward(W, W_hat);

  kmax = 2. / 3. * (N / 2 + 1);
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < Nf; k++) {
        z = (i * N + j) * Nf + k;
        dealias[z] = (fabs(kx[i]) < kmax) && (fabs(kx[j]) < kmax) &&
                     (fabs(kz[k]) < kmax);
      }

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < Nf; k++) {
        z = (i * N + j) * Nf + k;
        m = kx[i] * kx[i] + kx[j] * kx[j] + kz[k] * kz[k];
        kk[z] = m > 0 ? m : 1;
      }
  t = 0.0;
  tstep = 0;
  while (t <= T) {
    if (tstep % 2 == 0) {
      s = 0.0;
      for (k = 0; k < N * N * N; k++)
        s += U[k] * U[k] + V[k] * V[k] + W[k] * W[k];
      s *= 0.5 * dx * dx * dx / L / L / L;
      fprintf(stderr, "k = %.16e\n", s);
    }
    t += dt;
    tstep++;
    memcpy(U_hat0, U_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(V_hat0, V_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(W_hat0, W_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(U_hat1, U_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(V_hat1, V_hat, sizeof(fftw_complex) * N * N * Nf);
    memcpy(W_hat1, W_hat, sizeof(fftw_complex) * N * N * Nf);
    for (rk = 0; rk < 4; rk++) {
      if (rk > 0) {
        backward(U_hat, U);
        backward(V_hat, V);
        backward(W_hat, W);
        for (k = 0; k < N * N * N; k++) {
          U[k] /= tot;
          V[k] /= tot;
          W[k] /= tot;
        }
      }
      for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < Nf; k++) {
            z = (i * N + j) * Nf + k;
            curlZ[z] = I * (kx[i] * V_hat[z] - kx[j] * U_hat[z]);
            curlY[z] = I * (kz[k] * U_hat[z] - kx[i] * W_hat[z]);
            curlX[z] = I * (kx[j] * W_hat[z] - kz[k] * V_hat[z]);
          }
      backward(curlX, CU);
      backward(curlY, CV);
      backward(curlZ, CW);
      for (k = 0; k < N * N * N; k++) {
        CU[k] /= tot;
        CV[k] /= tot;
        CW[k] /= tot;
      }
      for (k = 0; k < N * N * N; k++) {
        U_tmp[k] = V[k] * CW[k] - W[k] * CV[k];
        V_tmp[k] = W[k] * CU[k] - U[k] * CW[k];
        W_tmp[k] = U[k] * CV[k] - V[k] * CU[k];
      }
      forward(U_tmp, dU);
      forward(V_tmp, dV);
      forward(W_tmp, dW);

      for (k = 0; k < N * N * Nf; k++) {
        dU[k] *= dealias[k] * dt;
        dV[k] *= dealias[k] * dt;
        dW[k] *= dealias[k] * dt;
      }
      for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < Nf; k++) {
            z = (i * N + j) * Nf + k;
            P_hat[z] = (dU[z] * kx[i] + dV[z] * kx[j] + dW[z] * kz[k]) / kk[z];
            dU[z] -= P_hat[z] * kx[i] + nu * dt * kk[z] * U_hat[z];
            dV[z] -= P_hat[z] * kx[j] + nu * dt * kk[z] * V_hat[z];
            dW[z] -= P_hat[z] * kz[k] + nu * dt * kk[z] * W_hat[z];
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
  fftw_free(U_hat0);
  fftw_free(U_hat1);
  fftw_free(U_tmp);
  fftw_free(V);
  fftw_free(V_hat);
  fftw_free(V_hat0);
  fftw_free(V_hat1);
  fftw_free(V_tmp);
  fftw_free(W);
  fftw_free(W_hat);
  fftw_free(W_hat0);
  fftw_free(W_hat1);
  fftw_free(W_tmp);

  free(dealias);
  free(kx);
  free(kz);
  free(kk);
}
