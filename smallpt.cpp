#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>

#include "smallpt_ispc.h"
#include "smallpt_c.h"

using namespace ispc;

inline float clamp(float x)
{
  return std::fmin(1.0f, std::fmax(0.0f, x));
}

inline int toInt(double x)
{
  return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}

int main(int argc, char* argv[])
{
  int w = argc >= 2 ? atoi(argv[1]) : 320;
  int h = argc >= 3 ? atoi(argv[2]) : 240;
  int samps4 = argc >= 4 ? atoi(argv[3]) : 8;
  int mode = argc >= 5 ? atoi(argv[4]) : 0;
  int mt = argc >= 6 ? atoi(argv[5]) : 0;

  printf("Rendering %d spp at %dx%d in %s (%s) mode...\n", 4 * samps4, w, h, mode == 0 ? "ISPC" : "MSVC", mt == 0 ? "single-threaded" : "multi-threaded");

  std::vector<float> image(w * h * 3, 0.0f);

  clock_t begin = clock();
  
  if (mode == 0) {
    run(image.data(), w, h, samps4, mt);
  } else {
    run_c(image.data(), w, h, samps4, mt);
  }
  
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  printf("Total time: %.1f s (%.2f ms/sample)\n", elapsed_secs, elapsed_secs * 1000. / (4 * samps4));

  FILE *f = fopen("image.pfm", "wb");
  fprintf(f, "PF\n%d %d\n%d\n", w, h, -1);
  fwrite(image.data(), sizeof(float), w*h * 3, f);
  fclose(f);
}
