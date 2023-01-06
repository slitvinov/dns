#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <stdlib.h>

enum { n0 = 6, n1 = 4, n2 = 3, m = n0 * n1 * (n2 / 2 + 1) };
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
  ptrdiff_t alloc_local;
  ptrdiff_t local_n0;
  ptrdiff_t local_0_start;

  MPI_CALL(MPI_Init(&argc, &argv));
  MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  fftw_mpi_init();
  /* fprintf(stderr, "rank = %d\n", rank); */
  alloc_local = fftw_mpi_local_size_3d(n0, n1, n2 / 2 + 1, MPI_COMM_WORLD,
                                       &local_n0, &local_0_start);
  fprintf(stderr,
          "%d: alloc_local, local_n0, local_0_start: %ld %ld "
          "%ld\n",
          rank, alloc_local, local_n0, local_0_start);

  if ((real = fftw_alloc_real(2 * alloc_local)) == NULL) {
    fprintf(stderr, "%s:%d: fftw_alloc_real failed\n", __FILE__, __LINE__);
    exit(1);
  }
  if ((compl = fftw_alloc_complex(alloc_local)) == NULL) {
    fprintf(stderr, "%s:%d: fftw_alloc_complex failed\n", __FILE__, __LINE__);
    exit(1);
  }

  flags = FFTW_MPI_TRANSPOSED_OUT;
  if ((plan_r2c = fftw_mpi_plan_dft_r2c_3d(n0, n1, n2, real, compl,
                                           MPI_COMM_WORLD, flags)) == NULL) {
    fprintf(stderr, "%s:%d: fftw_mpi_plan_dft_r2c_3d failed\n", __FILE__,
            __LINE__);
    exit(1);
  }
  fftw_destroy_plan(plan_r2c);
  fftw_free(real);
  fftw_free(compl );
  fftw_mpi_cleanup();
  MPI_CALL(MPI_Finalize());
}
