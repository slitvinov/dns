#include <complex.h>
#include "fftw3.h"
#include <math.h>
#include <stdlib.h>

enum { n0 = 5, n1 = 6, n2 = 7, m = n0 * n1 * (n2 / 2 + 1)};
#define pi 3.141592653589793238

int main(int argc, char **argv) {
  double *real;
  fftw_complex * compl ;
  unsigned flags;
  fftw_plan plan_r2c;
  fftw_plan plan_c2r;
  real = fftw_alloc_real(n0 * n1 * n2);
  compl = fftw_alloc_complex(m);

  flags = FFTW_DESTROY_INPUT;
  fftw_free(real);
  fftw_free(compl );
  fftw_destroy_plan(plan_c2r);
  fftw_destroy_plan(plan_r2c);
  fftw_cleanup();
}
