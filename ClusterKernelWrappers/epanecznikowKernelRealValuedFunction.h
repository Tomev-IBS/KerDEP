#ifndef DEDSTA_EPANECZNIKOWKERNELREALVALUEDFUNCTION_H
#define DEDSTA_EPANECZNIKOWKERNELREALVALUEDFUNCTION_H

#include "ClusterKernelsKDE/include/ClusterKernelsKDE/RealValuedFunction.h"

class EpanecznikowKernelRealValuedFunction : public RealValuedFunction {
  Point GetValue(const Point &x);
};

#endif //DEDSTA_EPANECZNIKOWKERNELREALVALUEDFUNCTION_H
