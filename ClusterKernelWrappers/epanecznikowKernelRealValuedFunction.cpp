#include <cmath>
#include "epanecznikowKernelRealValuedFunction.h"

Point EpanecznikowKernelRealValuedFunction::GetValue(const Point &x) {
  // This is one-dimensional Epanecznikow kernel.
  if(x.size() != 1){
    return {};
  }

  if(fabs(x[0]) <= 1){
    return {(3.0 / 4.0) * (1.0 - pow(x[0], 2))};
  }

  return {0};
}