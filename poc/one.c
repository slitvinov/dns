#include "fftw3.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>

enum { n = 10, m = 2 * (n / 2 + 1) };

int main(int argc, char **argv) {
  int i;
  double real[m];
  fftw_complex compl [m];
  unsigned flags;
  fftw_plan plan_r2c;
  fftw_plan plan_c2r;

  for (i = 0; i < m; i++)
    real[i] = i + 1;

  flags = 0;
  plan_r2c = fftw_plan_dft_r2c_1d(n, real, compl, flags);
  plan_c2r = fftw_plan_dft_c2r_1d(n, compl, real, flags);

  fftw_execute(plan_r2c);
  fftw_execute(plan_c2r);

  for (i = 0; i < n; i++)
    printf("%g\n", real[i] / n);

  fftw_cleanup();
}
