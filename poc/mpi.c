#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <stdlib.h>

enum { n0 = 3, n1 = 4, n2 = 5, m = n0 * n1 * (n2 / 2 + 1) };
#define pi 3.141592653589793238
#define MPI_CALL(command)                                                      \
  if ((ierr = (command)) != MPI_SUCCESS) {                                     \
    MPI_Error_string(ierr, err_buffer, &resultlen);                            \
    fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, err_buffer);            \
    exit(1);                                                                   \
  }

int main(int argc, char **argv) {
  char err_buffer[MPI_MAX_ERROR_STRING];
  int resultlen;
  int rank;
  double im;
  double om;
  double re;
  double *real;
  fftw_complex * compl ;
  fftw_complex *refer;
  fftw_plan plan_c2r;
  fftw_plan plan_r2c;
  int ierr;
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

  MPI_CALL(MPI_Init(&argc, &argv));
  MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  fftw_mpi_init();

  fprintf(stderr, "randk = %d\n", rank);
  /*  n = fftw_mpi_local_size_3d_transposed(N, N, Nf, MPI_COMM_WORLD, &n0, &s0, &n1,
      &s1); */
  fftw_mpi_cleanup();
  MPI_CALL(MPI_Finalize());
}
