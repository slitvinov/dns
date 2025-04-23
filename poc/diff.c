#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { n = 64, nf = n / 2 + 1, n3 = n * n * n, n3f = n * n * nf };
int main(void) {
  long i, j, k, l;
  double *u, *kx, *kz, dx, L;
  double pi = 3.141592653589793238;
  fftw_complex *u_hat;
  fftw_plan fplan, bplan;
  u = fftw_alloc_real(n3);
  u_hat = fftw_alloc_complex(n3f);
  kx = malloc(n * sizeof(double));
  kz = malloc(nf * sizeof(double));

  if ((fplan = fftw_plan_dft_r2c_3d(
           n, n, n, u, u_hat, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT)) == NULL) {
    fprintf(stderr, "%s:%d: error: fftw_plan_dft_r2c_3d failed\n", __FILE__,
            __LINE__);
    exit(1);
  }
  if ((bplan = fftw_plan_dft_c2r_3d(n, n, n, u_hat, u, FFTW_ESTIMATE)) ==
      NULL) {
    fprintf(stderr, "%s:%d: error: fftw_plan_dft_c2r_3d failed\n", __FILE__,
            __LINE__);
    exit(1);
  }

  for (i = 0; i < n / 2; i++) {
    kx[i] = i;
    kz[i] = i;
  }
  kz[n / 2] = n / 2;
  for (i = -n / 2; i < 0; i++)
    kx[i + n] = i;

  L = 2 * pi;
  dx = L / n;
  for (i = l = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++, l++)
        u[l] = exp(sin(dx * k));

  fftw_execute_dft_r2c(fplan, u, u_hat);
  for (i = l = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < nf; k++, l++)
        u_hat[l] *= I * kz[k];
  fftw_execute_dft_c2r(bplan, u_hat, u);

  for (i = l = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++, l++)
        printf("%.16e %.16e %.16e %.16e\n", i * dx, j * dx, k * dx, u[l] / n3);

  fftw_free(u);
  fftw_free(u_hat);
  fftw_destroy_plan(fplan);
  fftw_destroy_plan(bplan);
  free(kx);
  free(kz);
}
