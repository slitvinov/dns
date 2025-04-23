#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdlib.h>

enum { n = 10, m = n / 2 + 1 };
#define pi 3.141592653589793238

int main(int argc, char **argv) {
  int i;
  int j;
  double om;
  double *real;
  fftw_complex *compl;
  fftw_complex *refer;
  unsigned flags;
  fftw_plan plan_r2c;
  fftw_plan plan_c2r;
  double im;
  double re;

  real = fftw_alloc_real(n);
  compl = fftw_alloc_complex(m);
  refer = fftw_alloc_complex(m);

  flags = FFTW_DESTROY_INPUT;
  if ((plan_r2c = fftw_plan_dft_r2c_1d(n, real, compl, flags)) == NULL) {
    fprintf(stderr, "%s:%d: fftw_plan_dft_r2c_1d failed\n", __FILE__, __LINE__);
    exit(1);
  }
  if ((plan_c2r = fftw_plan_dft_c2r_1d(n, compl, real, flags)) == NULL) {
    fprintf(stderr, "%s:%d: fftw_plan_dft_c2r_1d failed\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < n; i++)
    real[i] = i;
  fftw_execute(plan_r2c);
  for (i = 0; i < m; i++) {
    im = 0;
    re = 0;
    for (j = 0; j < n; j++) {
      om = -2 * pi * j * i / n;
      re += real[j] * cos(om);
      im += real[j] * sin(om);
    }
    refer[i] = CMPLX(re, im);
  }

  for (i = 0; i < m; i++)
    printf("%d %+.16e %+.16e\n", i, cimag(compl[i]) - cimag(refer[i]),
           creal(compl[i]) - creal(refer[i]));

  fftw_execute(plan_c2r);
  for (i = 0; i < n; i++)
    printf("%g\n", real[i] / n);

  fftw_free(real);
  fftw_free(compl);
  fftw_free(refer);
  fftw_destroy_plan(plan_c2r);
  fftw_destroy_plan(plan_r2c);
  fftw_cleanup();
}
