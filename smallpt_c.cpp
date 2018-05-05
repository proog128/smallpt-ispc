#include "smallpt_c.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>

__declspec(align(16)) struct float3 { float x, y, z; };

inline float3 operator*(float3 const& l, float3 const& r) { return { l.x * r.x, l.y * r.y, l.z * r.z }; };
inline float3 operator*(float f, float3 const& v) { return { f*v.x, f*v.y,  f*v.z }; };
inline float3 operator*(float3 const& v, float f) { return f * v; };
inline float3 operator+(float3 const& l, float3 const& r) { return { l.x + r.x, l.y + r.y, l.z + r.z }; };
inline float3 operator-(float3 const& l, float3 const& r) { return { l.x - r.x, l.y - r.y, l.z - r.z }; };
inline float3 operator-(float3 const& v) { return { -v.x, -v.y, -v.z }; };

struct RNGState { unsigned int z1, z2, z3, z4; };

inline float floatbits(unsigned int x) {
  union { float f; unsigned int v; } t;
  t.v = x;
  return t.f;
}

inline void seed_rng(RNGState* state, unsigned int seed) {
  state->z1 = seed;
  state->z2 = seed ^ 0xbeeff00d;
  state->z3 = ((seed & 0xfffful) << 16) | (seed >> 16);
  state->z4 = (((seed & 0xfful) << 24) | ((seed & 0xff00ul) << 8) |
    ((seed & 0xff0000ul) >> 8) | (seed & 0xff000000ul) >> 24);
}

inline unsigned int random(RNGState* state)
{
  unsigned int b;

  b = ((state->z1 << 6) ^ state->z1) >> 13;
  state->z1 = ((state->z1 & 4294967294U) << 18) ^ b;
  b = ((state->z2 << 2) ^ state->z2) >> 27;
  state->z2 = ((state->z2 & 4294967288U) << 2) ^ b;
  b = ((state->z3 << 13) ^ state->z3) >> 21;
  state->z3 = ((state->z3 & 4294967280U) << 7) ^ b;
  b = ((state->z4 << 3) ^ state->z4) >> 12;
  state->z4 = ((state->z4 & 4294967168U) << 13) ^ b;
  return (state->z1 ^ state->z2 ^ state->z3 ^ state->z4);
}

inline float frandom(RNGState* state)
{
  unsigned int irand = random(state);
  irand &= (1ul << 23) - 1;
  return floatbits(0x3F800000 | irand) - 1.0f;
}

const float eps = 0.5;

inline float3 make_float3(float x, float y, float z) {
  float3 f = { x, y, z };
  return f;
}

inline float dot(float3 a, float3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline float3 norm(float3 v) {
  float len2 = dot(v, v);
  float invlen = 1.0f / sqrt(len2);
  v = invlen * v;
  return v;
}

inline float3 cross(float3 v0, float3 v1) {
  float3 ret;
  ret.x = v0.y * v1.z - v0.z * v1.y;
  ret.y = v0.z * v1.x - v0.x * v1.z;
  ret.z = v0.x * v1.y - v0.y * v1.x;
  return ret;
}

inline float max(float3 v) {
  return fmax(v.x, fmax(v.y, v.z));
}

struct Ray {
  float3 org;
  float3 dir;
  float3 weight;
};

enum Refl_t { DIFF, SPEC, REFR };

struct Sphere {
  float radius;
  float3 position;
  float3 emission;
  float3 color;
  Refl_t refl;
};

Sphere spheres[] = {
  {1e5, {1e5+1,40.8,81.6},      {0.0, 0.0, 0.0},    { .75,.25,.25 },       DIFF },//Left 
  {1e5, {-1e5+99,40.8,81.6 },   {0.0, 0.0, 0.0},    { .25,.25,.75 },       DIFF },//Rght 
  {1e5, {50,40.8, 1e5},         {0.0, 0.0, 0.0},    { .75,.75,.75 },       DIFF },//Back 
  {1e5, {50,40.8,-1e5 + 170 },  {0.0, 0.0, 0.0},    {0.0, 0.0, 0.0},       DIFF },//Frnt 
  {1e5, {50, 1e5, 81.6 },       {0.0, 0.0, 0.0},    { .75,.75,.75 },       DIFF },//Botm 
  {1e5, {50,-1e5 + 81.6,81.6 }, {0.0, 0.0, 0.0},    { .75,.75,.75 },       DIFF },//Top 
  {16.5,{27,16.5,47 },          {0.0, 0.0, 0.0},    {0.999, 0.999, 0.999}, SPEC },//Mirr 
  {16.5,{73,16.5,78 },          {0.0, 0.0, 0.0},    {0.999, 0.999, 0.999}, REFR },//Glas 
  {13., {50,81.6,81.6 },        {12.0, 12.0, 12.0}, {0.0, 0.0, 0.0},       DIFF } //Light
};

float intersect(const Ray& r, const Sphere& s)
{
  float3 op = s.position - r.org;
  float t;
  float b = dot(op, r.dir);
  float det = b * b - dot(op, op) + s.radius * s.radius;
  if (det < 0) {
    return 0.0f;
  } else {
    det = sqrt(det);
  }
  return (t = b - det)>eps ? t : ((t = b + det)>eps ? t : 0);
};

bool intersect(const Ray& r, float& t, int& id)
{
  int n = sizeof(spheres) / sizeof(Sphere);
  float d;
  t = INFINITY;
  for (int i = n; i--;) {
    if ((d = intersect(r, spheres[i])) && d < t) {
      t = d;
      id = i;
    }
  }
  return t<INFINITY;
}

float3 shade(Ray& ray, float t, int id, int depth, RNGState* Xi)
{
  const float3 x = ray.org + t * ray.dir;
  const float3 n = norm(x - spheres[id].position);
  const float3 nl = dot(n, ray.dir) < 0.0f ? n : -n;

  ray.org.x = INFINITY; // terminate path

  const float3 R = spheres[id].emission * ray.weight;

  float3 f = spheres[id].color;
  float max_refl = max(f);
  if(max_refl > 0.0f) {
    if (depth > 5) {
      if (frandom(Xi) < max_refl) {
        f = (1 / max_refl) * f;
      } else {
        return R;
      }
    }

    Refl_t refl = spheres[id].refl;
    if (refl == DIFF) {
      const float r1 = 2.0f * M_PI * frandom(Xi);
      const float r2 = frandom(Xi);
      const float r2s = sqrt(r2);

      const float3 w = nl;
      const float3 u = norm(cross(abs(w.x) > .1 ? make_float3(0, 1, 0) : make_float3(1, 0, 0), w));
      const float3 v = cross(w, u);

      const float3 d = norm(u*cos(r1)*r2s + v * sin(r1)*r2s + w * sqrt(1 - r2));

      ray.org = x;
      ray.dir = d;
      ray.weight = ray.weight * f;
    } else if (refl == SPEC) {
      const float3 d = ray.dir - n * 2 * dot(n, ray.dir);
      ray.org = x;
      ray.dir = d;
      ray.weight = ray.weight * f;
    } else {
      float3 reflDir = ray.dir - n * 2 * dot(n, ray.dir);
      const bool into = dot(n, nl) > 0;
      const float nc = 1, nt = 1.5;
      const float nnt = into ? nc / nt : nt / nc;
      const float ddn = dot(ray.dir, nl);
      const float cos2t = 1 - nnt * nnt*(1 - ddn * ddn);
      if (cos2t < 0) {  // TIR
        ray.org = x;
        ray.dir = reflDir;
        ray.weight = ray.weight * f;
      } else {
        const float3 tdir = norm(ray.dir*nnt - n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t))));
        const float a = nt - nc, b = nt + nc;
        const float R0 = a * a / (b*b);
        const float c = 1 - (into ? -ddn : dot(tdir, n));
        const float Re = R0 + (1 - R0)*c*c*c*c*c;
        const float Tr = 1 - Re;
        const float P = .25f + .5f*Re;
        if (frandom(Xi) < P) {
          const float RP = Re / P;
          ray.org = x;
          ray.dir = reflDir;
          ray.weight = ray.weight * f * RP;
        } else {
          const float TP = Tr / (1 - P);
          ray.org = x;
          ray.dir = tdir;
          ray.weight = ray.weight * f * TP;
        }
      }
    }
  }

  return R;
}

void render(float image[], int x0, int x1, int y0, int y1, int width, int height, int samples)
{
  float3 camDir = { 0.0f, -0.042612f, -1.0f };
  Ray cam = { { 50.0f, 52.0f, 295.6f }, norm(camDir) };

  float3 cx = { width * 0.5135f / height, 0.0f, 0.0f };
  float3 cy = norm(cross(cx, cam.dir));
  cy = cy * 0.5135f;

  for(int y=y0; y<y1; ++y) {
    for (int x = x0; x < x1; ++x) {
      RNGState Xi;
      seed_rng(&Xi, y * width + x);

      for (int sy = 0; sy < 2; sy++) {      // 2x2 subpixel rows
        for (int sx = 0; sx < 2; sx++) {    // 2x2 subpixel cols
          for (int s = 0; s < samples; ++s) {
            const float r1 = 2 * frandom(&Xi);
            const float dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
            const float r2 = 2 * frandom(&Xi);
            const float dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

            float3 d = norm(cx * ((float)((sx + .5 + dx) / 2 + x) / width - 0.5f) + cy * ((float)((sy + .5 + dy) / 2 + y) / height - 0.5f) + cam.dir);
            Ray ray = { cam.org + d * 140.0f, d, make_float3(1.0, 1.0, 1.0) };

            int depth = 0;

            while (ray.org.x != INFINITY) {
              float t;
              int id;
              if (!intersect(ray, t, id)) break;

              float3 color = shade(ray, t, id, depth, &Xi);

              const float w = 0.25f * (1.0f / samples);
              int offset = 3 * (y * width + x);
              image[offset] += color.x * w;
              image[offset + 1] += color.y * w;
              image[offset + 2] += color.z * w;

              ++depth;
            }
          }
        }
      }
    }
  }
}

void run_c(float* image, int32_t width, int32_t height, int32_t samples, bool mt)
{
  if (mt) {
    printf("Not supported.\n");
  } else {
    render(image, 0, width, 0, height, width, height, samples);
  }
}
