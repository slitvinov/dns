#include <complex.h>
#include "fftw3.h"
#include <math.h>
#include <stdlib.h>

enum { n = 10, m = 2 * (n / 2 + 1) };

int main(int argc, char **argv) {
  int i;
  double *real;
  fftw_complex * compl ;
  unsigned flags;
  fftw_plan plan_r2c;
  fftw_plan plan_c2r;

  real = fftw_alloc_real(n);
  compl = fftw_alloc_complex(m);
  for (i = 0; i < n; i++)
    real[i] = i;
  flags = 0;
  plan_r2c = fftw_plan_dft_r2c_1d(n, real, compl, flags);
  plan_c2r = fftw_plan_dft_c2r_1d(n, compl, real, flags);
  fftw_execute(plan_r2c);

  for (i = 0; i < n / 2; i++)
    printf("%d %+.16e %+.16e\n", i, 2 * compl [i][0] / n, 2 * compl [i][1] / n);

  fftw_execute(plan_c2r);
  for (i = 0; i < n; i++)
    printf("%g\n", real[i] / n);

  fftw_free(real);
  fftw_free(compl );
  fftw_destroy_plan(plan_c2r);
  fftw_destroy_plan(plan_r2c);
  fftw_cleanup();
}
