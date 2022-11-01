#include "fftw3-mpi.h"
#include <complex>
#include <math.h>
#include <vector>

using namespace std;
typedef double precision;

enum { M = 7, N = 1 << M, Nf = N / 2 + 1, tot = N * N * N};
int main(int argc, char *argv[]) {
  int rank;
  double L, dx;
  precision nu, dt, T;
  precision pi = 3.141592653589793238;
  ptrdiff_t alloc_local, local_n0, local_0_start, local_n1, local_1_start;
  int i;
  int j;
  int k;
  int m;
  int rk;
  int tstep;
  int z;
  precision s_in[1];
  precision s_out[1];
  precision a[4];
  precision b[3];
  precision kx[N];
  precision kz[Nf];

  MPI::Init(argc, argv);
  fftw_mpi_init();
  rank = MPI::COMM_WORLD.Get_rank();
  nu = 0.000625;
  T = 0.1;
  dt = 0.01;
  L = 2 * pi;
  dx = L / N;
  a[0] = 1. / 6.;
  a[1] = 1. / 3.;
  a[2] = 1. / 3.;
  a[3] = 1. / 6.;
  b[0] = 0.5;
  b[1] = 0.5;
  b[2] = 1.0;
  alloc_local = fftw_mpi_local_size_3d_transposed(N, N, Nf, MPI::COMM_WORLD,
                                                  &local_n0, &local_0_start,
                                                  &local_n1, &local_1_start);

  vector<precision> U(2 * alloc_local);
  vector<precision> V(2 * alloc_local);
  vector<precision> W(2 * alloc_local);
  vector<precision> U_tmp(2 * alloc_local);
  vector<precision> V_tmp(2 * alloc_local);
  vector<precision> W_tmp(2 * alloc_local);
  vector<precision> CU(2 * alloc_local);
  vector<precision> CV(2 * alloc_local);
  vector<precision> CW(2 * alloc_local);
  vector<int> dealias(2 * alloc_local);
  vector<precision> kk(2 * alloc_local);
  vector<complex<precision>> U_hat(alloc_local);
  vector<complex<precision>> V_hat(alloc_local);
  vector<complex<precision>> W_hat(alloc_local);
  vector<complex<precision>> P_hat(alloc_local);
  vector<complex<precision>> U_hat0(alloc_local);
  vector<complex<precision>> V_hat0(alloc_local);
  vector<complex<precision>> W_hat0(alloc_local);
  vector<complex<precision>> U_hat1(alloc_local);
  vector<complex<precision>> V_hat1(alloc_local);
  vector<complex<precision>> W_hat1(alloc_local);
  vector<complex<precision>> dU(alloc_local);
  vector<complex<precision>> dV(alloc_local);
  vector<complex<precision>> dW(alloc_local);
  vector<complex<precision>> curlX(alloc_local);
  vector<complex<precision>> curlY(alloc_local);
  vector<complex<precision>> curlZ(alloc_local);

  // Starting time
  MPI::COMM_WORLD.Barrier();

  for (i = 0; i < N / 2; i++) {
    kx[i] = i;
    kz[i] = i;
  }
  kz[N / 2] = N / 2;
  for (i = -N / 2; i < 0; i++)
    kx[i + N] = i;

  // fftw_plan plan_backward;
  fftw_plan rfftn, irfftn;
  rfftn = fftw_mpi_plan_dft_r2c_3d(
      N, N, N, U.data(), reinterpret_cast<fftw_complex *>(U_hat.data()),
      MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_OUT);
  irfftn = fftw_mpi_plan_dft_c2r_3d(
      N, N, N, reinterpret_cast<fftw_complex *>(U_hat.data()), U.data(),
      MPI::COMM_WORLD, FFTW_MPI_TRANSPOSED_IN);

  for (i = 0; i < local_n0; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++) {
        z = (i * N + j) * 2 * Nf + k;
        U[z] = sin(dx * (i + local_0_start)) * cos(dx * j) * cos(dx * k);
        V[z] = -cos(dx * (i + local_0_start)) * sin(dx * j) * cos(dx * k);
        W[z] = 0.0;
      }

  fftw_mpi_execute_dft_r2c(rfftn, U.data(),
                           reinterpret_cast<fftw_complex *>(U_hat.data()));
  fftw_mpi_execute_dft_r2c(rfftn, V.data(),
                           reinterpret_cast<fftw_complex *>(V_hat.data()));
  fftw_mpi_execute_dft_r2c(rfftn, W.data(),
                           reinterpret_cast<fftw_complex *>(W_hat.data()));

  precision kmax = 2. / 3. * (N / 2 + 1);
  for (i = 0; i < local_n1; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < Nf; k++) {
        z = (i * N + j) * Nf + k;
        dealias[z] = (abs(kx[i + local_1_start]) < kmax) * (abs(kx[j]) < kmax) *
                     (abs(kx[k]) < kmax);
      }

  for (i = 0; i < local_n1; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < Nf; k++) {
        z = (i * N + j) * Nf + k;
        m = kx[i + local_1_start] * kx[i + local_1_start] + kx[j] * kx[j] +
                kx[k] * kx[k];
        kk[z] = m > 0 ? m : 1;
      }

  complex<precision> one(0, 1);
  double t = 0.0;
  tstep = 0;
  while (t < T - 1e-8) {
    t += dt;
    tstep++;
    for (i = 0; i < local_n1; i++)
      for (j = 0; j < N; j++)
        for (k = 0; k < Nf; k++) {
          z = (i * N + j) * Nf + k;
          U_hat0[z] = U_hat[z];
          V_hat0[z] = V_hat[z];
          W_hat0[z] = W_hat[z];
          U_hat1[z] = U_hat[z];
          V_hat1[z] = V_hat[z];
          W_hat1[z] = W_hat[z];
        }
    for (rk = 0; rk < 4; rk++) {
      if (rk > 0) {
        fftw_mpi_execute_dft_c2r(
            irfftn, reinterpret_cast<fftw_complex *>(U_hat.data()), U.data());
        fftw_mpi_execute_dft_c2r(
            irfftn, reinterpret_cast<fftw_complex *>(V_hat.data()), V.data());
        fftw_mpi_execute_dft_c2r(
            irfftn, reinterpret_cast<fftw_complex *>(W_hat.data()), W.data());
        for (k = 0; k < (int)U.size(); k++) {
          U[k] /= tot;
          V[k] /= tot;
          W[k] /= tot;
        }
      }
      // Compute curl
      for (i = 0; i < local_n1; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < Nf; k++) {
            z = (i * N + j) * Nf + k;
            curlZ[z] =
                one * (kx[i + local_1_start] * V_hat[z] - kx[j] * U_hat[z]);
            curlY[z] =
                one * (kz[k] * U_hat[z] - kx[i + local_1_start] * W_hat[z]);
            curlX[z] = one * (kx[j] * W_hat[z] - kz[k] * V_hat[z]);
          }
      fftw_mpi_execute_dft_c2r(
          irfftn, reinterpret_cast<fftw_complex *>(curlX.data()), CU.data());
      fftw_mpi_execute_dft_c2r(
          irfftn, reinterpret_cast<fftw_complex *>(curlY.data()), CV.data());
      fftw_mpi_execute_dft_c2r(
          irfftn, reinterpret_cast<fftw_complex *>(curlZ.data()), CW.data());
      for (k = 0; k < (int)CU.size(); k++) {
        CU[k] /= tot;
        CV[k] /= tot;
        CW[k] /= tot;
      }

      // Cross
      for (i = 0; i < local_n0; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < N; k++) {
            z = (i * N + j) * 2 * Nf + k;
            U_tmp[z] = V[z] * CW[z] - W[z] * CV[z];
            V_tmp[z] = W[z] * CU[z] - U[z] * CW[z];
            W_tmp[z] = U[z] * CV[z] - V[z] * CU[z];
          }

      fftw_mpi_execute_dft_r2c(rfftn, U_tmp.data(),
                               reinterpret_cast<fftw_complex *>(dU.data()));
      fftw_mpi_execute_dft_r2c(rfftn, V_tmp.data(),
                               reinterpret_cast<fftw_complex *>(dV.data()));
      fftw_mpi_execute_dft_r2c(rfftn, W_tmp.data(),
                               reinterpret_cast<fftw_complex *>(dW.data()));

      for (i = 0; i < local_n1; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < Nf; k++) {
            z = (i * N + j) * Nf + k;
            dU[z] *= (dealias[z] * dt);
            dV[z] *= (dealias[z] * dt);
            dW[z] *= (dealias[z] * dt);
          }
      //
      for (i = 0; i < local_n1; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < Nf; k++) {
            z = (i * N + j) * Nf + k;
            P_hat[z] = (dU[z] * kx[i + local_1_start] + dV[z] * kx[j] +
                        dW[z] * kz[k]) /
                       kk[z];
            dU[z] -=
                (P_hat[z] * kx[i + local_1_start] + nu * dt * kk[z] * U_hat[z]);
            dV[z] -= (P_hat[z] * kx[j] + nu * dt * kk[z] * V_hat[z]);
            dW[z] -= (P_hat[z] * kz[k] + nu * dt * kk[z] * W_hat[z]);
          }

      if (rk < 3) {
        for (i = 0; i < local_n1; i++)
          for (j = 0; j < N; j++)
            for (k = 0; k < Nf; k++) {
              z = (i * N + j) * Nf + k;
              U_hat[z] = U_hat0[z] + b[rk] * dU[z];
              V_hat[z] = V_hat0[z] + b[rk] * dV[z];
              W_hat[z] = W_hat0[z] + b[rk] * dW[z];
            }
      }
      for (i = 0; i < local_n1; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < Nf; k++) {
            z = (i * N + j) * Nf + k;
            U_hat1[z] += a[rk] * dU[z];
            V_hat1[z] += a[rk] * dV[z];
            W_hat1[z] += a[rk] * dW[z];
          }
    }
    for (i = 0; i < local_n1; i++)
      for (j = 0; j < N; j++)
        for (k = 0; k < Nf; k++) {
          z = (i * N + j) * Nf + k;
          U_hat[z] = U_hat1[z];
          V_hat[z] = V_hat1[z];
          W_hat[z] = W_hat1[z];
        }

    if (tstep % 2 == 0) {
      s_in[0] = 0.0;
      for (i = 0; i < local_n0; i++)
        for (j = 0; j < N; j++)
          for (k = 0; k < N; k++) {
            z = (i * N + j) * 2 * Nf + k;
            s_in[0] += (U[z] * U[z] + V[z] * V[z] + W[z] * W[z]);
          }
      s_in[0] *= (0.5 * dx * dx * dx / L / L / L);

      MPI::COMM_WORLD.Reduce(s_in, s_out, 1, MPI::DOUBLE,
                             MPI::SUM, 0);
      if (rank == 0)
	fprintf(stderr, "k = %.16e\n", s_out[0]);
    }
  }
  MPI::Finalize();
}
