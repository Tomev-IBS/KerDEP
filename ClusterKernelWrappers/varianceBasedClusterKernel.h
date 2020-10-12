#ifndef KERDEP_VARIANCEBASEDCLUSTERKERNEL_H
#define KERDEP_VARIANCEBASEDCLUSTERKERNEL_H

#include "ClusterKernelsKDE/include/ClusterKernelsKDE/ClusterKernel.h"
#include "ClusterKernelWrappers/epanecznikowKernelRealValuedFunction.h"

class VarianceBasedClusterKernel : public ClusterKernel {
  public:
    VarianceBasedClusterKernel(ClusterKernelStreamElement *stream_element, RealValuedFunction *kernel = new EpanecznikowKernelRealValuedFunction());
    VarianceBasedClusterKernel(unsigned int cardinality, Point elements_sum, Point elements_squared_sum, RealValuedFunction *kernel = new EpanecznikowKernelRealValuedFunction());
    Point GetMean() override;
    void Update(ClusterKernelStreamElement *stream_element) override;
    ClusterKernel* Merge(ClusterKernel *other_cluster_kernel) override;
    Point GetValue(const Point &pt) override;
    double GetWeight() override;
    unsigned int GetCardinality() override;
    Point GetElementsSum();
    Point GetElementsSquaredSum();
    void SetBandwidth(const Point &bandwidth) override;
  private:
    unsigned int cardinality_ = 0;
    Point elements_sum_ = {};
    Point elements_squared_sum_ = {};
    RealValuedFunctionPtr kernel_;
    Point bandwidth_ = {};
};

#endif //KERDEP_VARIANCEBASEDCLUSTERKERNEL_H
