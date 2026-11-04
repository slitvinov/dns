/* standalone renderer for volumes dumped by tgr -v 1
   usage: rview vol.raw out.ppm [iso_frac=0.35] [W=1024] [ss=3] [ao=1]
                [fog=0.5] [az=35] [el=28] [dist=3] */
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "render.h"
static double wtime(void) {
  struct timespec q;
  clock_gettime(CLOCK_MONOTONIC, &q);
  return q.tv_sec + 1e-9 * q.tv_nsec;
}
int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: rview vol.raw out.ppm [iso_frac W ss ao fog az el dist]\n");
    return 1;
  }
  double frac = argc > 3 ? atof(argv[3]) : 0.35;
  int W = argc > 4 ? atoi(argv[4]) : 1024, ss = argc > 5 ? atoi(argv[5]) : 3,
      ao = argc > 6 ? atoi(argv[6]) : 1;
  double fog = argc > 7 ? atof(argv[7]) : 0.5, az = argc > 8 ? atof(argv[8]) : 35,
         el = argc > 9 ? atof(argv[9]) : 28, dist = argc > 10 ? atof(argv[10]) : 3;
  FILE *f = fopen(argv[1], "rb");
  if (!f) { perror(argv[1]); return 1; }
  fseek(f, 0, SEEK_END);
  long n = ftell(f) / sizeof(float4), L = lround(cbrt((double)n));
  rewind(f);
  if (L * L * L != n) { fprintf(stderr, "rview: not a cube: %ld float4\n", n); return 1; }
  float4 *h = (float4 *)malloc(n * sizeof(float4));
  if (fread(h, sizeof(float4), n, f) != (size_t)n) { perror("fread"); return 1; }
  fclose(f);
  float fmax = 0;
  for (long i = 0; i < n; i++) fmax = fmaxf(fmax, h[i].w);
  cudaChannelFormatDesc cd = cudaCreateChannelDesc<float4>();
  cudaArray_t arr;
  cudaMalloc3DArray(&arr, &cd, make_cudaExtent(L, L, L));
  cudaMemcpy3DParms cp;
  memset(&cp, 0, sizeof cp);
  cp.srcPtr = make_cudaPitchedPtr(h, L * sizeof(float4), L, L);
  cp.dstArray = arr;
  cp.extent = make_cudaExtent(L, L, L);
  cp.kind = cudaMemcpyHostToDevice;
  cudaMemcpy3D(&cp);
  cudaTextureObject_t tl, tp;
  maketex(arr, &tl, &tp);
  uchar4 *dimg, *himg = (uchar4 *)malloc((size_t)W * W * sizeof(uchar4));
  unsigned long long *dnf, nf;
  cudaMalloc(&dimg, (size_t)W * W * sizeof(uchar4));
  cudaMalloc(&dnf, sizeof nf);
  RParam rp = {(float)L, (float)(frac * frac * fmax), W, W, ss, ao, (float)fog};
  Cam cam = orbit(az, el, dist);
  double t0 = 0, t1 = 0;
  for (int rep = 0; rep < 2; rep++) { /* second run is timed (warm) */
    cudaMemset(dnf, 0, sizeof nf);
    cudaDeviceSynchronize();
    t0 = wtime();
    raymarch<<<dim3((W + 15) / 16, (W + 15) / 16), dim3(16, 16)>>>(tl, tp, 0, rp, cam, dimg, dnf);
    cudaDeviceSynchronize();
    t1 = wtime();
  }
  cudaError_t e = cudaGetLastError();
  if (e != cudaSuccess) { fprintf(stderr, "rview: %s\n", cudaGetErrorString(e)); return 1; }
  cudaMemcpy(himg, dimg, (size_t)W * W * sizeof(uchar4), cudaMemcpyDeviceToHost);
  cudaMemcpy(&nf, dnf, sizeof nf, cudaMemcpyDeviceToHost);
  FILE *o = fopen(argv[2], "wb");
  fprintf(o, "P6\n%d %d\n255\n", W, W);
  unsigned char *rgb = (unsigned char *)himg;
  for (long q = 0; q < (long)W * W; q++) {
    rgb[3 * q] = himg[q].x; rgb[3 * q + 1] = himg[q].y; rgb[3 * q + 2] = himg[q].z;
  }
  fwrite(rgb, 3, (size_t)W * W, o);
  fclose(o);
  fprintf(stderr, "rview: L=%ld max|w|=%.3f iso=%.2f max W=%d ss=%d ao=%d: %.1f ms, %.0f fetch/ray\n",
          L, sqrt(fmax), frac, W, ss, ao, 1e3 * (t1 - t0), (double)nf / ((double)W * W * ss * ss));
  return 0;
}
