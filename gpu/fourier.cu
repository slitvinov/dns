#include <cufft.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { nvars = 4, nthread = 256, nblock = 1024 };
static const double pi = 3.141592653589793238;
static const double a[] = {1 / 6.0, 1 / 3.0, 1 / 3.0, 1 / 6.0};
static const double b[] = {0.5, 0.5, 1.0};
static void cuda(cudaError_t e, const char *what) {
  if (e != cudaSuccess) {
    fprintf(stderr, "fourier: error: %s: %s\n", what, cudaGetErrorString(e));
    exit(1);
  }
}
static void cufft(cufftResult r, const char *what) {
  if (r != CUFFT_SUCCESS) {
    fprintf(stderr, "fourier: error: %s: cufft error %d\n", what, (int)r);
    exit(1);
  }
}
static void launched(const char *what) { cuda(cudaGetLastError(), what); }
__device__ static double wave(long i, long n) { return i < n / 2 ? i : i - n; }
__global__ static void curl(long n, long nf, long n3f, const double2 *U,
                            const double2 *V, const double2 *W, double2 *cX,
                            double2 *cY, double2 *cZ, double2 *sU, double2 *sV,
                            double2 *sW) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n3f)
    return;
  long i = l / (n * nf), j = l / nf % n, k = l % nf;
  double ki = wave(i, n), kj = wave(j, n), kz = k;
  double2 u = U[l], v = V[l], w = W[l];
  cZ[l] = make_double2(-(ki * v.y - kj * u.y), ki * v.x - kj * u.x);
  cY[l] = make_double2(-(kz * u.y - ki * w.y), kz * u.x - ki * w.x);
  cX[l] = make_double2(-(kj * w.y - kz * v.y), kj * w.x - kz * v.x);
  sU[l] = u;
  sV[l] = v;
  sW[l] = w;
}
__global__ static void cross(long n3, double s, double *U, double *V, double *W,
                             const double *CU, const double *CV,
                             const double *CW) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n3)
    return;
  double u = U[l] * s, v = V[l] * s, w = W[l] * s;
  double cu = CU[l] * s, cv = CV[l] * s, cw = CW[l] * s;
  U[l] = v * cw - w * cv;
  V[l] = w * cu - u * cw;
  W[l] = u * cv - v * cu;
}
__global__ static void project(long n, long nf, long n3f, double kmax,
                               double nu, double dt, double ak, double bk,
                               int stage, double2 *dU, double2 *dV, double2 *dW,
                               double2 *P, double2 *U, double2 *V, double2 *W,
                               const double2 *U0, const double2 *V0,
                               const double2 *W0, double2 *U1, double2 *V1,
                               double2 *W1) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l >= n3f)
    return;
  long i = l / (n * nf), j = l / nf % n, k = l % nf;
  double ki = wave(i, n), kj = wave(j, n), kz = k;
  double kk = ki * ki + kj * kj + kz * kz;
  double f = (int)(fabs(ki) < kmax && fabs(kj) < kmax && fabs(kz) < kmax) * dt;
  double2 du = dU[l], dv = dV[l], dw = dW[l], p;
  du.x *= f;
  du.y *= f;
  dv.x *= f;
  dv.y *= f;
  dw.x *= f;
  dw.y *= f;
  if (kk > 0) {
    p.x = (du.x * ki + dv.x * kj + dw.x * kz) / kk;
    p.y = (du.y * ki + dv.y * kj + dw.y * kz) / kk;
  } else
    p = make_double2(0, 0);
  double2 u = U[l], v = V[l], w = W[l];
  double g = nu * dt * kk;
  du.x -= p.x * ki + g * u.x;
  du.y -= p.y * ki + g * u.y;
  dv.x -= p.x * kj + g * v.x;
  dv.y -= p.y * kj + g * v.y;
  dw.x -= p.x * kz + g * w.x;
  dw.y -= p.y * kz + g * w.y;
  P[l] = p;
  if (stage) {
    U[l] = make_double2(U0[l].x + bk * du.x, U0[l].y + bk * du.y);
    V[l] = make_double2(V0[l].x + bk * dv.x, V0[l].y + bk * dv.y);
    W[l] = make_double2(W0[l].x + bk * dw.x, W0[l].y + bk * dw.y);
  }
  U1[l].x += ak * du.x;
  U1[l].y += ak * du.y;
  V1[l].x += ak * dv.x;
  V1[l].y += ak * dv.y;
  W1[l].x += ak * dw.x;
  W1[l].y += ak * dw.y;
}
__global__ static void scale(long n3, double s, double *x) {
  long l = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (l < n3)
    x[l] *= s;
}
__global__ static void diag(long n, long nf, long n3f, const double2 *U,
                            const double2 *V, const double2 *W, double *part) {
  __shared__ double se[nthread], so[nthread];
  double e = 0, o = 0;
  for (long l = blockIdx.x * (long)blockDim.x + threadIdx.x; l < n3f;
       l += (long)gridDim.x * blockDim.x) {
    long i = l / (n * nf), j = l / nf % n, k = l % nf;
    double ki = wave(i, n), kj = wave(j, n), kz = k;
    double kk = ki * ki + kj * kj + kz * kz;
    double h = k == 0 || k == n / 2 ? 0.5 : 1;
    double2 u = U[l], v = V[l], w = W[l];
    double s = u.x * u.x + u.y * u.y + (v.x * v.x + v.y * v.y) +
               (w.x * w.x + w.y * w.y);
    e += h * s;
    o += h * kk * s;
  }
  se[threadIdx.x] = e;
  so[threadIdx.x] = o;
  __syncthreads();
  for (int m = blockDim.x / 2; m > 0; m /= 2) {
    if (threadIdx.x < m) {
      se[threadIdx.x] += se[threadIdx.x + m];
      so[threadIdx.x] += so[threadIdx.x + m];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    part[blockIdx.x] = se[0];
    part[gridDim.x + blockIdx.x] = so[0];
  }
}
__global__ static void total(int m, double *part) {
  __shared__ double se[nblock], so[nblock];
  se[threadIdx.x] = threadIdx.x < m ? part[threadIdx.x] : 0;
  so[threadIdx.x] = threadIdx.x < m ? part[m + threadIdx.x] : 0;
  __syncthreads();
  for (int k = blockDim.x / 2; k > 0; k /= 2) {
    if (threadIdx.x < k) {
      se[threadIdx.x] += se[threadIdx.x + k];
      so[threadIdx.x] += so[threadIdx.x + k];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    part[2 * m] = se[0];
    part[2 * m + 1] = so[0];
  }
}
static long steps(const char *opt, double x, double dt) {
  long n = lround(x / dt);
  if (fabs(n * dt - x) > 1e-9 * x) {
    fprintf(stderr, "fourier: error: %s %g is not a multiple of -s %g\n", opt, x,
            dt);
    exit(1);
  }
  return n;
}
int main(int argc, char **argv) {
  (void)argc;
  cufftHandle fplan, bplan;
  FILE *file;
  char path[FILENAME_MAX], *input_path, *end;
  long double energy, Omega;
  double dx, L, invn3, kmax, nu, dt, T, t, sum[2];
  double2 *curlX, *curlY, *curlZ, *dU, *dV, *dW, *P_hat, *U_hat, *U_hat0,
      *U_hat1, *V_hat, *V_hat0, *V_hat1, *W_hat, *W_hat0, *W_hat1, *swap;
  int rk, Verbose, Dump;
  long idump, tstep, nt;
  size_t offset, wf, wb;
  size_t ivar;
  double *CU, *CV, *CW, *U, *V, *W, *dump, *part;
  void *work;
  cudaDeviceProp prop;

  input_path = NULL;
  dt = -1;
  T = 0;
  nu = -1;
  Verbose = 0;
  Dump = 0;
  while (*++argv != NULL && argv[0][0] == '-') {
    switch (argv[0][1]) {
    case 'h':
      fprintf(stderr, "Usage: fourier [-v] [-d] -i <input.raw> -n <viscosity> -t "
                      "<end time> -s <time step>\n"
                      "\n"
                      "Options:\n"
                      "  -i <input.raw>    Input file\n"
                      "  -n <viscosity>    Viscosity\n"
                      "  -t <end time>     End time\n"
                      "  -s <time step>    Time step\n"
                      "  -v                Verbose output\n"
                      "  -d                Dump snapshots\n"
                      "  -h                Show this help message\n"
                      "\n"
                      "Example:\n"
                      "  fourier -i tgv.raw -n 0.01 -t 1.0 -s 0.001 -v\n");
      exit(1);
    case 'v':
      Verbose = 1;
      break;
    case 'd':
      Dump = 1;
      break;
    case 'i':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "fourier: error: -i needs an argument\n");
        exit(1);
      }
      input_path = *argv;
      break;
    case 'n':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "fourier: error: -n needs an argument\n");
        exit(1);
      }
      nu = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(stderr, "fourier: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      break;
    case 's':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "fourier: error: -s needs an argument\n");
        exit(1);
      }
      dt = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(stderr, "fourier: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      break;
    case 't':
      argv++;
      if (*argv == NULL) {
        fprintf(stderr, "fourier: error: -t needs an argument\n");
        exit(1);
      }
      T = strtod(*argv, &end);
      if (*end != '\0') {
        fprintf(stderr, "fourier: error: '%s' is not a number\n", *argv);
        exit(1);
      }
      break;
    default:
      fprintf(stderr, "fourier: error: unknown option '%s'\n", *argv);
      exit(1);
    }
  }
  if (T == 0) {
    fprintf(stderr, "fourier: error: -t is not set or invalid\n");
    exit(1);
  }
  if (nu == -1) {
    fprintf(stderr, "fourier: error: -n is not set or invalid\n");
    exit(1);
  }
  if (dt == -1) {
    fprintf(stderr, "fourier: error: -s is not set or invalid\n");
    exit(1);
  }
  if (input_path == NULL) {
    fprintf(stderr, "fourier: error: -i is not set\n");
    exit(1);
  }
  if ((file = fopen(input_path, "r")) == NULL) {
    fprintf(stderr, "fourier: error: fail to open '%s'\n", input_path);
    exit(1);
  }
  if (Verbose) {
    cuda(cudaGetDeviceProperties(&prop, 0), "cudaGetDeviceProperties");
    fprintf(stderr, "fourier: device: %s\n", prop.name);
  }
  if (fseek(file, 0, SEEK_END) != 0 || (offset = ftell(file)) == (size_t)-1) {
    fprintf(stderr, "fourier: error: fail to seek '%s'\n", input_path);
    exit(1);
  }
  rewind(file);
  long n = offset / sizeof(double) / nvars;
  n = round(powf(n, 1.0 / 3));
  if (n * n * n * nvars * sizeof(double) != offset) {
    fprintf(stderr, "fourier: error: wrong file '%s'\n", input_path);
    exit(1);
  }
  if (Verbose)
    fprintf(stderr, "fourier: n = %ld\n", n);
  long nf = n / 2 + 1;
  long n3 = n * n * n;
  long n3f = n * n * nf;
  int gr = (n3 + nthread - 1) / nthread;
  int gc = (n3f + nthread - 1) / nthread;
  if ((dump = (double *)malloc(n3 * sizeof(double))) == NULL) {
    fprintf(stderr, "fourier: error: fail to allocate host memory\n");
    exit(1);
  }
  cuda(cudaMalloc(&U, n3 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&V, n3 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&W, n3 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&CU, n3 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&CV, n3 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&CW, n3 * sizeof(double)), "cudaMalloc");
  cuda(cudaMalloc(&U_hat, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&V_hat, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&W_hat, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&P_hat, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&U_hat0, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&V_hat0, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&W_hat0, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&U_hat1, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&V_hat1, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&W_hat1, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&dU, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&dV, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&dW, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&curlX, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&curlY, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&curlZ, n3f * sizeof(double2)), "cudaMalloc");
  cuda(cudaMalloc(&part, (2 * nblock + 2) * sizeof(double)), "cudaMalloc");
  cuda(cudaMemset(P_hat, 0, n3f * sizeof(double2)), "cudaMemset");
  double *in[] = {U, V, W};
  for (ivar = 0; ivar < 3; ivar++) {
    if (fread(dump, sizeof(double), n3, file) != (size_t)n3) {
      fprintf(stderr, "fourier: error: fail to read '%s'\n", input_path);
      exit(1);
    }
    cuda(cudaMemcpy(in[ivar], dump, n3 * sizeof(double),
                    cudaMemcpyHostToDevice),
         "cudaMemcpy");
  }
  if (fclose(file) != 0) {
    fprintf(stderr, "fourier: error: fail to read '%s'\n", input_path);
    exit(1);
  }
  L = 2 * pi;
  dx = L / n;
  invn3 = 1.0 / n3;
  cufft(cufftCreate(&fplan), "cufftCreate");
  cufft(cufftCreate(&bplan), "cufftCreate");
  cufft(cufftSetAutoAllocation(fplan, 0), "cufftSetAutoAllocation");
  cufft(cufftSetAutoAllocation(bplan, 0), "cufftSetAutoAllocation");
  cufft(cufftMakePlan3d(fplan, n, n, n, CUFFT_D2Z, &wf), "cufftMakePlan3d");
  cufft(cufftMakePlan3d(bplan, n, n, n, CUFFT_Z2D, &wb), "cufftMakePlan3d");
  cuda(cudaMalloc(&work, wf > wb ? wf : wb), "cudaMalloc");
  cufft(cufftSetWorkArea(fplan, work), "cufftSetWorkArea");
  cufft(cufftSetWorkArea(bplan, work), "cufftSetWorkArea");
  kmax = 2. / 3. * (n / 2 + 1);

  cufft(cufftExecD2Z(fplan, U, U_hat), "cufftExecD2Z");
  cufft(cufftExecD2Z(fplan, V, V_hat), "cufftExecD2Z");
  cufft(cufftExecD2Z(fplan, W, W_hat), "cufftExecD2Z");

  idump = 0;
  t = 0.0;
  nt = steps("-t", T, dt);
  tstep = 0;
  for (;;) {
    if (tstep % 10 == 0) {
      diag<<<nblock, nthread>>>(n, nf, n3f, U_hat, V_hat, W_hat, part);
      launched("diag");
      total<<<1, nblock>>>(nblock, part);
      launched("total");
      cuda(cudaMemcpy(sum, part + 2 * nblock, 2 * sizeof(double),
                      cudaMemcpyDeviceToHost),
           "cudaMemcpy");
      energy = sum[0];
      Omega = sum[1];
      energy *= invn3 * invn3;
      Omega *= invn3 * invn3;
      printf("% 10ld % .16e % .16Le % .16Le\n", tstep, t, energy, Omega);
      if (fflush(stdout) != 0 || ferror(stdout)) {
        fprintf(stderr, "fourier: error: fail to write stdout\n");
        exit(1);
      }
      if (Dump) {
        double2 *list[nvars] = {U_hat, V_hat, W_hat, P_hat};
        const char *name[nvars] = {"U", "V", "W", "P"};
        sprintf(path, "%08ld.raw", tstep);
        if ((file = fopen(path, "w")) == NULL) {
          fprintf(stderr, "fourier: error: fail to open '%s'\n", path);
          exit(1);
        }
        for (ivar = 0; ivar < nvars; ivar++) {
          cuda(cudaMemcpy(curlX, list[ivar], n3f * sizeof(double2),
                          cudaMemcpyDeviceToDevice),
               "cudaMemcpy");
          cufft(cufftExecZ2D(bplan, curlX, CU), "cufftExecZ2D");
          scale<<<gr, nthread>>>(n3, invn3, CU);
          launched("scale");
          cuda(cudaMemcpy(dump, CU, n3 * sizeof(double),
                          cudaMemcpyDeviceToHost),
               "cudaMemcpy");
          if (fwrite(dump, sizeof(double), n3, file) != (size_t)n3) {
            fprintf(stderr, "fourier: error: fail to write '%s'\n", path);
            exit(1);
          }
        }
        if (ferror(file) || fclose(file) != 0) {
          fprintf(stderr, "fourier: error: fail to write '%s'\n", path);
          exit(1);
        }
        sprintf(path, "a.%08ld.xdmf2", idump);
        if ((file = fopen(path, "w")) == NULL) {
          fprintf(stderr, "fourier: error: fail to open '%s'\n", path);
          exit(1);
        }
        if (fprintf(file,
                    "<Xdmf\n"
                    "    Version=\"2\">\n"
                    "  <Domain>\n"
                    "    <Grid>\n"
                    "      <Time\n"
                    "          Value=\"%+.16e\"/>\n"
                    "      <Topology\n"
                    "          TopologyType=\"3DCoRectMesh\"\n"
                    "          Dimensions=\"%ld %ld %ld\"/>\n"
                    "      <Geometry\n"
                    "          GeometryType=\"ORIGIN_DXDYDZ\">\n"
                    "        <DataItem\n"
                    "            Dimensions=\"3\">\n"
                    "          0\n"
                    "          0\n"
                    "          0\n"
                    "        </DataItem>\n"
                    "        <DataItem\n"
                    "            Dimensions=\"3\">\n"
                    "          %.16e\n"
                    "          %.16e\n"
                    "          %.16e\n"
                    "        </DataItem>\n"
                    "      </Geometry>\n",
                    t, n, n, n, dx, dx, dx) < 0) {
          fprintf(stderr, "fourier: error: fail to write '%s'\n", path);
          exit(1);
        }
        offset = 0;
        for (ivar = 0; ivar < nvars; ivar++) {
          if (fprintf(file,
                      "      <Attribute\n"
                      "          name=\"%s\">\n"
                      "        <DataItem\n"
                      "            Format=\"Binary\"\n"
                      "            Seek=\"%ld\"\n"
                      "            Precision=\"8\"\n"
                      "            Dimensions=\"%ld %ld %ld\">\n"
                      "          %08ld.raw\n"
                      "        </DataItem>\n"
                      "      </Attribute>\n",
                      name[ivar], (long)offset, n, n, n, tstep) < 0) {
            fprintf(stderr, "fourier: error: fail to write '%s'\n", path);
            exit(1);
          }
          offset += n3 * sizeof(double);
        }
        if (fprintf(file, "    </Grid>\n"
                          "  </Domain>\n"
                          "</Xdmf>\n") < 0) {
          fprintf(stderr, "fourier: error: fail to write '%s'\n", path);
          exit(1);
        }
        if (ferror(file) || fclose(file) != 0) {
          fprintf(stderr, "fourier: error: fail to write '%s'\n", path);
          exit(1);
        }
        idump++;
      }
    }
    if (tstep >= nt)
      break;
    cuda(cudaMemcpy(U_hat0, U_hat, n3f * sizeof(double2),
                    cudaMemcpyDeviceToDevice),
         "cudaMemcpy");
    cuda(cudaMemcpy(V_hat0, V_hat, n3f * sizeof(double2),
                    cudaMemcpyDeviceToDevice),
         "cudaMemcpy");
    cuda(cudaMemcpy(W_hat0, W_hat, n3f * sizeof(double2),
                    cudaMemcpyDeviceToDevice),
         "cudaMemcpy");
    cuda(cudaMemcpy(U_hat1, U_hat, n3f * sizeof(double2),
                    cudaMemcpyDeviceToDevice),
         "cudaMemcpy");
    cuda(cudaMemcpy(V_hat1, V_hat, n3f * sizeof(double2),
                    cudaMemcpyDeviceToDevice),
         "cudaMemcpy");
    cuda(cudaMemcpy(W_hat1, W_hat, n3f * sizeof(double2),
                    cudaMemcpyDeviceToDevice),
         "cudaMemcpy");
    for (rk = 0; rk < 4; rk++) {
      curl<<<gc, nthread>>>(n, nf, n3f, U_hat, V_hat, W_hat, curlX, curlY,
                            curlZ, dU, dV, dW);
      launched("curl");
      cufft(cufftExecZ2D(bplan, dU, U), "cufftExecZ2D");
      cufft(cufftExecZ2D(bplan, dV, V), "cufftExecZ2D");
      cufft(cufftExecZ2D(bplan, dW, W), "cufftExecZ2D");
      cufft(cufftExecZ2D(bplan, curlX, CU), "cufftExecZ2D");
      cufft(cufftExecZ2D(bplan, curlY, CV), "cufftExecZ2D");
      cufft(cufftExecZ2D(bplan, curlZ, CW), "cufftExecZ2D");
      cross<<<gr, nthread>>>(n3, invn3, U, V, W, CU, CV, CW);
      launched("cross");
      cufft(cufftExecD2Z(fplan, U, dU), "cufftExecD2Z");
      cufft(cufftExecD2Z(fplan, V, dV), "cufftExecD2Z");
      cufft(cufftExecD2Z(fplan, W, dW), "cufftExecD2Z");
      project<<<gc, nthread>>>(n, nf, n3f, kmax, nu, dt, a[rk],
                               rk < 3 ? b[rk] : 0, rk < 3, dU, dV, dW, P_hat,
                               U_hat, V_hat, W_hat, U_hat0, V_hat0, W_hat0,
                               U_hat1, V_hat1, W_hat1);
      launched("project");
    }
    swap = U_hat;
    U_hat = U_hat1;
    U_hat1 = swap;
    swap = V_hat;
    V_hat = V_hat1;
    V_hat1 = swap;
    swap = W_hat;
    W_hat = W_hat1;
    W_hat1 = swap;
    tstep++;
    t = tstep * dt;
  }
  cuda(cudaDeviceSynchronize(), "cudaDeviceSynchronize");
  cufft(cufftDestroy(fplan), "cufftDestroy");
  cufft(cufftDestroy(bplan), "cufftDestroy");
  cuda(cudaFree(work), "cudaFree");
  cuda(cudaFree(U), "cudaFree");
  cuda(cudaFree(V), "cudaFree");
  cuda(cudaFree(W), "cudaFree");
  cuda(cudaFree(CU), "cudaFree");
  cuda(cudaFree(CV), "cudaFree");
  cuda(cudaFree(CW), "cudaFree");
  cuda(cudaFree(U_hat), "cudaFree");
  cuda(cudaFree(V_hat), "cudaFree");
  cuda(cudaFree(W_hat), "cudaFree");
  cuda(cudaFree(P_hat), "cudaFree");
  cuda(cudaFree(U_hat0), "cudaFree");
  cuda(cudaFree(V_hat0), "cudaFree");
  cuda(cudaFree(W_hat0), "cudaFree");
  cuda(cudaFree(U_hat1), "cudaFree");
  cuda(cudaFree(V_hat1), "cudaFree");
  cuda(cudaFree(W_hat1), "cudaFree");
  cuda(cudaFree(dU), "cudaFree");
  cuda(cudaFree(dV), "cudaFree");
  cuda(cudaFree(dW), "cudaFree");
  cuda(cudaFree(curlX), "cudaFree");
  cuda(cudaFree(curlY), "cudaFree");
  cuda(cudaFree(curlZ), "cudaFree");
  cuda(cudaFree(part), "cudaFree");
  free(dump);
}
