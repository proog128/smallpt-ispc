#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>

#include "smallpt_ispc.h"
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
  int w = 320, h = 240;
  int samps = argc == 2 ? atoi(argv[1]) / 4 : 32;

  std::vector<float> image(w * h * 3, 0.0f);

  clock_t begin = clock();
  
  run(image.data(), w, h, samps);
  
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  printf("%.4f ms", elapsed_secs * 1000. / (4 * samps));

  FILE *f = fopen("image.pfm", "wb");
  fprintf(f, "PF\n%d %d\n%d\n", w, h, -1);
  fwrite(image.data(), sizeof(float), w*h * 3, f);
  fclose(f);

  getchar();
}
