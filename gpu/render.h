/* Isosurface ray caster for f = |omega|^2 stored as float4 (grad f, f) on
   L^3 nodes spanning the box [0,pi]^3. Shared by tgr.cu (in situ) and
   rview.cu (standalone, from a dumped volume). */
#include <math.h>
#include <stdint.h>
__host__ __device__ static float3 add3(float3 a, float3 b) {
  return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
__host__ __device__ static float3 mul3(float3 a, float s) {
  return make_float3(a.x * s, a.y * s, a.z * s);
}
__host__ __device__ static float dot3(float3 a, float3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
__host__ __device__ static float3 nrm3(float3 a) {
  return mul3(a, 1.0f / sqrtf(dot3(a, a)));
}
__host__ __device__ static float3 crs3(float3 a, float3 b) {
  return make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                     a.x * b.y - a.y * b.x);
}
__host__ __device__ static float3 mix3(float3 a, float3 b, float t) {
  return add3(mul3(a, 1 - t), mul3(b, t));
}
struct Cam {
  float3 eye, f, r, u;
  float th;
};
struct RParam {
  float L, iso;   /* nodes per axis, threshold on f */
  int W, H, ss;   /* image size, ss x ss samples per pixel */
  int ao;         /* ambient occlusion on/off */
  float fog;      /* depth fog strength */
};
static Cam orbit(double az_deg, double el_deg, double dist) {
  Cam c;
  double az = az_deg * M_PI / 180, el = el_deg * M_PI / 180,
         D = dist * M_PI;
  float3 C = make_float3(M_PI / 2, M_PI / 2, M_PI / 2);
  c.eye = add3(C, make_float3(D * cos(el) * cos(az), D * cos(el) * sin(az),
                              D * sin(el)));
  c.f = nrm3(add3(C, mul3(c.eye, -1)));
  c.r = nrm3(crs3(c.f, make_float3(0, 0, 1)));
  c.u = crs3(c.r, c.f);
  c.th = tan(15 * M_PI / 180);
  return c;
}
/* hardware trilinear (9-bit weights): cheap, for marching and AO */
__device__ static float4 fetchl(cudaTextureObject_t t, float3 p, float s) {
  return tex3D<float4>(t, p.x * s + 0.5f, p.y * s + 0.5f, p.z * s + 0.5f);
}
/* fp32 value from point samples: sum_q w_q (f_q + 1/2 g_q . (p - x_q)),
   exact for any quadratic f because sum_q w_q (p - x_q) = 0, so ridge
   peaks between nodes are kept. Used for the surface itself. */
__device__ static float fval(cudaTextureObject_t t, float3 p, float s,
                             float L) {
  float x = fminf(fmaxf(p.x * s, 0.f), L - 1.0001f),
        y = fminf(fmaxf(p.y * s, 0.f), L - 1.0001f),
        z = fminf(fmaxf(p.z * s, 0.f), L - 1.0001f);
  float ix = floorf(x), iy = floorf(y), iz = floorf(z), a = x - ix, b = y - iy,
        c = z - iz, r = 0;
  for (int q = 0; q < 8; q++) {
    int qx = q & 1, qy = q >> 1 & 1, qz = q >> 2;
    float wq = (qx ? a : 1 - a) * (qy ? b : 1 - b) * (qz ? c : 1 - c);
    float4 v = tex3D<float4>(t, ix + qx + 0.5f, iy + qy + 0.5f, iz + qz + 0.5f);
    r += wq * (v.w + 0.5f * (v.x * (a - qx) + v.y * (b - qy) + v.z * (c - qz)) / s);
  }
  return r;
}
/* C2 tricubic B-spline of grad f from 8 hardware trilinear fetches (Sigg &
   Hadwiger, GPU Gems 2 ch. 20): smooth normals without cell creases */
__device__ static float3 gradbs(cudaTextureObject_t t, float3 p, float s) {
  float u[3] = {p.x * s, p.y * s, p.z * s}, h0[3], h1[3], g0[3];
  for (int k = 0; k < 3; k++) {
    float i = floorf(u[k]), a = u[k] - i, a2 = a * a, a3 = a2 * a;
    float w0 = (1 - 3 * a + 3 * a2 - a3) / 6, w1 = (4 - 6 * a2 + 3 * a3) / 6,
          w2 = (1 + 3 * a + 3 * a2 - 3 * a3) / 6, w3 = a3 / 6;
    g0[k] = w0 + w1;
    h0[k] = i - 1 + w1 / (w0 + w1) + 0.5f;
    h1[k] = i + 1 + w3 / (w2 + w3) + 0.5f;
  }
  float3 r = make_float3(0, 0, 0);
  for (int q = 0; q < 8; q++) {
    int qx = q & 1, qy = q >> 1 & 1, qz = q >> 2;
    float w = (qx ? 1 - g0[0] : g0[0]) * (qy ? 1 - g0[1] : g0[1]) *
              (qz ? 1 - g0[2] : g0[2]);
    float4 v = tex3D<float4>(t, qx ? h1[0] : h0[0], qy ? h1[1] : h0[1],
                             qz ? h1[2] : h0[2]);
    r = add3(r, mul3(make_float3(v.x, v.y, v.z), w));
  }
  return r;
}
__device__ static unsigned hashu(unsigned x) {
  x ^= x >> 16; x *= 0x7feb352d; x ^= x >> 15; x *= 0x846ca68b; x ^= x >> 16;
  return x;
}
/* march one ray; returns 0 miss, 1 surface, 2 box cap; *tout, *nout */
__device__ static int trace(cudaTextureObject_t tl, cudaTextureObject_t tp,
                            const RParam &rp, float3 o, float3 d, float *tout,
                            float3 *nout, int *nf) {
  const float B = (float)M_PI, s = (rp.L - 1) / B, h = 0.5f / s;
  float e[3] = {o.x, o.y, o.z}, v[3] = {d.x, d.y, d.z}, t0 = -1e30f, t1 = 1e30f;
  int ea = 0;
  for (int a = 0; a < 3; a++) {
    float inv = 1 / v[a], ta = (0 - e[a]) * inv, tb = (B - e[a]) * inv;
    float lo = fminf(ta, tb), hi = fmaxf(ta, tb);
    if (lo > t0) {
      t0 = lo;
      ea = a;
    }
    t1 = fminf(t1, hi);
  }
  if (t1 <= fmaxf(t0, 0))
    return 0;
  float t = fmaxf(t0, 0);
  (*nf)++;
  if (fval(tp, add3(o, mul3(d, t)), s, rp.L) >= rp.iso) {
    float nn[3] = {0, 0, 0};
    nn[ea] = v[ea] > 0 ? -1 : 1;
    *nout = make_float3(nn[0], nn[1], nn[2]);
    *tout = t;
    return 2;
  }
  float tp_ = t;
  while (t < t1) {
    float3 pt = add3(o, mul3(d, t));
    float fl = fetchl(tl, pt, s).w;
    (*nf)++;
    /* full step far from the surface, quarter step when getting close */
    float st = fl > 0.3f * rp.iso ? 0.25f * h : h;
    if (fl >= 0.4f * rp.iso && ((*nf)++, fval(tp, pt, s, rp.L) >= rp.iso)) {
      float a = tp_, b = t;
      for (int q = 0; q < 14; q++) {
        float m = 0.5f * (a + b);
        if (fval(tp, add3(o, mul3(d, m)), s, rp.L) >= rp.iso)
          b = m;
        else
          a = m;
      }
      *tout = b;
      *nout = nrm3(mul3(gradbs(tl, add3(o, mul3(d, b)), s), -1));
      return 1;
    }
    tp_ = t;
    t += st;
  }
  return 0;
}
/* tz: optional linear texture of omega_z; surfaces are red where
   omega_z > 0 and blue where omega_z < 0 (height split if tz == 0) */
__global__ static void raymarch(cudaTextureObject_t tl, cudaTextureObject_t tp,
                                cudaTextureObject_t tz, RParam rp, Cam cam, uchar4 *img,
                                unsigned long long *nfetch) {
  int px = blockIdx.x * blockDim.x + threadIdx.x,
      py = blockIdx.y * blockDim.y + threadIdx.y;
  if (px >= rp.W || py >= rp.H)
    return;
  const float B = (float)M_PI, s = (rp.L - 1) / B;
  /* background: soft vertical gradient with vignette */
  float vy = (float)py / rp.H, vx = (float)px / rp.W;
  float vig = 1 - 0.35f * ((vx - 0.5f) * (vx - 0.5f) + (vy - 0.5f) * (vy - 0.5f));
  float3 bg = mul3(mix3(make_float3(0.026f, 0.030f, 0.042f),
                        make_float3(0.004f, 0.004f, 0.006f), vy), vig);
  float3 C = make_float3(B / 2, B / 2, B / 2), oc = add3(cam.eye, mul3(C, -1));
  float dc = sqrtf(dot3(oc, oc)), rad = 0.87f * B, tn = dc - rad, tf = dc + rad;
  float3 key = nrm3(add3(mul3(cam.f, -0.55f), add3(mul3(cam.u, 0.75f), mul3(cam.r, -0.45f))));
  float3 fill = nrm3(add3(mul3(cam.f, -0.3f), add3(mul3(cam.u, -0.5f), mul3(cam.r, 0.8f))));
  float3 acc = make_float3(0, 0, 0);
  int nf = 0;
  unsigned seed = hashu(py * 9781u + px * 6271u);
  for (int sy = 0; sy < rp.ss; sy++)
    for (int sx = 0; sx < rp.ss; sx++) {
      float ox = (sx + 0.5f) / rp.ss, oy = (sy + 0.5f) / rp.ss;
      float qx = (2 * (px + ox) / rp.W - 1) * cam.th * rp.W / rp.H,
            qy = (1 - 2 * (py + oy) / rp.H) * cam.th;
      float3 d = nrm3(add3(cam.f, add3(mul3(cam.r, qx), mul3(cam.u, qy))));
      float t;
      float3 n;
      int hit = trace(tl, tp, rp, cam.eye, d, &t, &n, &nf);
      float3 col = bg;
      if (hit) {
        float3 p = add3(cam.eye, mul3(d, t));
        if (dot3(n, d) > 0)
          n = mul3(n, -1);
        /* ambient occlusion: 12 hemisphere directions x 3 radii */
        float ao = 1;
        if (rp.ao && hit == 1) {
          float3 a1 = fabsf(n.x) < 0.9f ? make_float3(1, 0, 0) : make_float3(0, 1, 0);
          float3 tu = nrm3(crs3(n, a1)), tv = crs3(n, tu);
          float occ = 0, wsum = 0, rot = 6.2831853f * (hashu(seed++) & 0xffff) / 65536.f;
          for (int k = 0; k < 12; k++) {
            float ph = rot + 2.39996323f * k, ct = sqrtf((k + 0.5f) / 12), st_ = sqrtf(1 - ct * ct);
            float3 dir = add3(mul3(n, ct), add3(mul3(tu, st_ * cosf(ph)), mul3(tv, st_ * sinf(ph))));
            for (int j = 0; j < 3; j++) {
              float r = (j == 0 ? 2.f : j == 1 ? 5.f : 12.f) / s;
              float fo = fetchl(tl, add3(p, mul3(dir, r)), s).w;
              float w = 1.f / (j + 1);
              occ += w * fminf(fmaxf((fo - 0.5f * rp.iso) / rp.iso, 0.f), 1.f);
              wsum += w;
            }
          }
          ao = 1 - 0.85f * occ / wsum;
        }
        float wz = tz ? tex3D<float>(tz, p.x * s + 0.5f, p.y * s + 0.5f, p.z * s + 0.5f)
                      : p.z - 0.5f * B;
        float3 base = wz > 0 ? make_float3(0.72f, 0.04f, 0.03f) : make_float3(0.02f, 0.10f, 0.70f);
        if (hit == 2)
          base = mul3(base, 0.5f);
        float dk = fmaxf(dot3(n, key), 0), df = fmaxf(dot3(n, fill), 0);
        float3 hv = nrm3(add3(key, mul3(d, -1)));
        float sp = powf(fmaxf(dot3(n, hv), 0), 60) * (hit == 1 ? 0.6f : 0.1f) * ao;
        float rim = powf(1 - fabsf(dot3(n, d)), 3) * 0.05f;
        col = add3(mul3(base, ao * (0.06f + 0.18f * df) + 0.9f * dk * (0.35f + 0.65f * ao)),
                   make_float3(sp + rim, sp + rim, sp + rim));
        float fz = fminf(fmaxf((t - tn) / (tf - tn), 0.f), 1.f);
        col = mix3(col, bg, rp.fog * fz * fz);
      }
      acc = add3(acc, col);
    }
  acc = mul3(acc, 1.f / (rp.ss * rp.ss));
  /* gamma 2.2 for display */
  img[py * rp.W + px] = make_uchar4(
      (unsigned char)(255 * powf(fminf(acc.x, 1.f), 1 / 2.2f)),
      (unsigned char)(255 * powf(fminf(acc.y, 1.f), 1 / 2.2f)),
      (unsigned char)(255 * powf(fminf(acc.z, 1.f), 1 / 2.2f)), 255);
  atomicAdd(nfetch, (unsigned long long)nf);
}
static void maketex(cudaArray_t arr, cudaTextureObject_t *tl,
                    cudaTextureObject_t *tp) {
  cudaResourceDesc rd;
  cudaTextureDesc td;
  memset(&rd, 0, sizeof rd);
  rd.resType = cudaResourceTypeArray;
  rd.res.array.array = arr;
  for (int q = 0; q < 2; q++) {
    memset(&td, 0, sizeof td);
    td.addressMode[0] = td.addressMode[1] = td.addressMode[2] = cudaAddressModeClamp;
    td.filterMode = q == 0 ? cudaFilterModeLinear : cudaFilterModePoint;
    td.readMode = cudaReadModeElementType;
    td.normalizedCoords = 0;
    cudaCreateTextureObject(q == 0 ? tl : tp, &rd, &td, NULL);
  }
}
