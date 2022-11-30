#include <assert.h>
#include "WolfMCExecutor.h"
#include <iostream>
#include <fstream>

int main(int argc, char **argv) {
  LocalDOF<2> s1({1.0 / 2.0, sin(M_PI / 3)});
  assert(abs(s1.Cos6Theta() - 1.0) < 1e-15);
  LocalDOF<2> s2({cos(M_PI / 6), 1.0 / 2.0});
  assert(abs(s2.Cos6Theta() + 1.0) < 1e-15);
  LocalDOF<2> s3({cos(M_PI / 12), sin(M_PI / 12)});
  assert(abs(s3.Cos6Theta()) < 1e-15);

  return 0;
}