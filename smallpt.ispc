typedef float<3> float3;
typedef float<4> float4;

const int inf = 0x7f800000;
const float pi = 0x1.921fb54442d18p+1;
const float eps = 0.5;

inline float3 make_float3(float x, float y, float z) {
  float3 f = { x, y, z };
  return f;
}

inline float dot(float3 a, float3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline uniform float dot(uniform float3 a, uniform float3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline float3 norm(float3 v) {
  float len2 = dot(v, v);
  float invlen = rsqrt(len2);
  v *= invlen;
  return v;
}

inline uniform float3 norm(uniform float3 v) {
  uniform float len2 = dot(v, v);
  uniform float invlen = rsqrt(len2);
  v *= invlen;
  return v;
}

inline uniform float3 cross(uniform float3 v0, uniform float3 v1) {
  uniform float3 ret;
  ret.x = v0.y * v1.z - v0.z * v1.y;
  ret.y = v0.z * v1.x - v0.x * v1.z;
  ret.z = v0.x * v1.y - v0.y * v1.x;
  return ret;
}

inline float3 cross(float3 v0, float3 v1) {
  float3 ret;
  ret.x = v0.y * v1.z - v0.z * v1.y;
  ret.y = v0.z * v1.x - v0.x * v1.z;
  ret.z = v0.x * v1.y - v0.y * v1.x;
  return ret;
}

inline float max(float3 v) {
  return max(v.x, max(v.y, v.z));
}

struct Ray {
  float3 org;
  float3 dir;
  float3 weight;
};

enum Refl_t { DIFF, SPEC, REFR };

struct Sphere {
  uniform float radius;
  uniform float3 position;
  uniform float3 emission;
  uniform float3 color;
  uniform Refl_t refl;
};

uniform Sphere spheres[] = {
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

float intersect(const Ray& r, uniform const Sphere& s)
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
  uniform int n = sizeof(spheres) / sizeof(Sphere);
  float d;
  t = floatbits(inf);
  for (uniform int i = n; i--;) {
    if ((d = intersect(r, spheres[i])) && d < t) {
      t = d;
      id = i;
    }
  }
  return t<inf;
}

float3 shade(Ray& ray, float t, int id, int depth, varying RNGState * uniform Xi)
{
  const float3 x = ray.org + t * ray.dir;
  const float3 n = norm(x - spheres[id].position);
  const float3 nl = dot(n, ray.dir) < 0.0f ? n : -n;

  ray.org.x = floatbits(inf); // terminate path

  const float3 R = spheres[id].emission * ray.weight;

  float3 f = spheres[id].color;
  float max_refl = max(f);
  cif (max_refl > 0.0f) {
    if (depth > 5) {
      if (frandom(Xi) < max_refl) {
        f = f * (1 / max_refl);
      } else {
        return R;
      }
    }

    Refl_t refl = spheres[id].refl;
    if (refl == DIFF) {
      const float r1 = 2.0 * pi * frandom(Xi);
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
      uniform const float nc = 1, nt = 1.5;
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
        const float P = .25 + .5*Re;
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

void render(uniform float image[], uniform int x0, uniform int x1, uniform int y0, uniform int y1, uniform int width, uniform int height, uniform int samples)
{
  uniform float3 camDir = { 0.0f, -0.042612f, -1.0f };
  uniform Ray cam = { { 50.0f, 52.0f, 295.6f }, norm(camDir) };

  uniform float3 cx = { width * 0.5135 / height, 0.0f, 0.0f };
  uniform float3 cy = norm(cross(cx, cam.dir));
  cy *= 0.5135;

  foreach_tiled(y = y0 ... y1, x = x0 ... x1) {
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

          while (ray.org.x != floatbits(inf)) {
            float t;
            int id;
            if (!intersect(ray, t, id)) break;

            float3 color = shade(ray, t, id, depth, &Xi);

            const float w = 0.25 * (1.0 / samples);
            int offset = 3 * (y * width + x);
            image[offset] += color.r * w;
            image[offset + 1] += color.g * w;
            image[offset + 2] += color.b * w;

            ++depth;
          }
        }
      }
    }
  }
}

void task run_task(uniform float image[], uniform int width, uniform int height, uniform int samples)
{
  const uniform int x0 = min(taskIndex0 * 32, width);
  const uniform int x1 = min((taskIndex0 + 1) * 32, width);
  const uniform int y0 = min(taskIndex1 * 32, height);
  const uniform int y1 = min((taskIndex1 + 1) * 32, height);

  render(image, x0, x1, y0, y1, width, height, samples);
}

export void run(uniform float image[], uniform int width, uniform int height, uniform int samples)
{
  launch[(width-1)/32+1, (height-1)/32+1] run_task(image, width, height, samples);
  //render(image, 0, width, 0, height, width, height, samples);
}
