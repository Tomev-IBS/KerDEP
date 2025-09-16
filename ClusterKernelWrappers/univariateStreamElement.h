#ifndef DEDSTA_UNIVARIATESTREAMELEMENT_H
#define DEDSTA_UNIVARIATESTREAMELEMENT_H

#include "ClusterKernelsKDE/include/ClusterKernelsKDE/ClusterKernelStreamElement.h"

class UnivariateStreamElement : public ClusterKernelStreamElement {
  public:
    UnivariateStreamElement(const Point &pt);
    Point GetMean() override;
  private:
    Point coordinates_;
};

#endif //DEDSTA_UNIVARIATESTREAMELEMENT_H
