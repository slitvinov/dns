#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { n0 = 2, n1 = 2, n2 = 2, m = n0 * n1 * (n2 / 2 + 1) };
#define pi 3.141592653589793238

int main(int argc, char **argv) {
  double im;
  double om;
  double re;
  double *real;
  fftw_complex *compl;
  fftw_complex *refer;
  fftw_plan plan_c2r;
  fftw_plan plan_r2c;
  int i;
  int k;
  int l;
  int i0;
  int i1;
  int i2;
  int j0;
  int j1;
  int j2;
  unsigned flags;

  real = fftw_alloc_real(n0 * n1 * n2);
  compl = fftw_alloc_complex(m);
  refer = fftw_alloc_complex(m);

  flags = FFTW_ESTIMATE;
  if ((plan_r2c = fftw_plan_dft_r2c_3d(n0, n1, n2, real, compl, flags)) ==
      NULL) {
    fprintf(stderr, "%s:%d: fftw_plan_dft_r2c_3d\n", __FILE__, __LINE__);
    exit(1);
  }
  if ((plan_c2r = fftw_plan_dft_c2r_3d(n0, n1, n2, compl, real, flags)) ==
      NULL) {
    fprintf(stderr, "%s:%d: fftw_plan_dft_c2r_3d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < n0 * n1 * n2; i++)
    real[i] = i + 1;
  k = 0;
  for (i0 = 0; i0 < n0; i0++)
    for (i1 = 0; i1 < n1; i1++)
      for (i2 = 0; i2 < n2 / 2 + 1; i2++) {
        im = 0;
        re = 0;
        l = 0;
        for (j0 = 0; j0 < n0; j0++)
          for (j1 = 0; j1 < n1; j1++)
            for (j2 = 0; j2 < n2; j2++) {
              om = -pi * (2.0 * j0 * i0 / n0 + 2.0 * j1 * i1 / n1 +
                          2.0 * j2 * i2 / n2);
              re += real[l] * cos(om);
              im += real[l] * sin(om);
              l++;
            }
        refer[k] = CMPLX(re, im);
        k++;
      }
  fftw_execute_dft_r2c(plan_r2c, real, compl);
  for (i = 0; i < m; i++)
    printf("%d %+.16e %+.16e\n", i, cimag(compl[i]) - cimag(refer[i]),
           creal(compl[i]) - creal(refer[i]));

  memset(real, 0, n0 * n1 * n2 * sizeof(double));
  fftw_execute_dft_r2c(plan_c2r, real, compl);
  for (i = 0; i < n0 * n1 * n2; i++)
    printf("%g\n", real[i] / n0 / n1 / n2);

  fftw_free(real);
  fftw_free(compl);
  fftw_free(refer);
  fftw_destroy_plan(plan_c2r);
  fftw_destroy_plan(plan_r2c);
  fftw_cleanup();
}
