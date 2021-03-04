#ifndef KERDEP_EPANECZNIKOWKERNELREALVALUEDFUNCTION_H
#define KERDEP_EPANECZNIKOWKERNELREALVALUEDFUNCTION_H

#include "ClusterKernelsKDE/include/ClusterKernelsKDE/RealValuedFunction.h"

class EpanecznikowKernelRealValuedFunction : public RealValuedFunction {
  Point GetValue(const Point &x);
};

#endif //KERDEP_EPANECZNIKOWKERNELREALVALUEDFUNCTION_H
