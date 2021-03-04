#ifndef KERDEP_UNIVARIATESTREAMELEMENT_H
#define KERDEP_UNIVARIATESTREAMELEMENT_H

#include "ClusterKernelsKDE/include/ClusterKernelsKDE/ClusterKernelStreamElement.h"

class UnivariateStreamElement : public ClusterKernelStreamElement {
  public:
    UnivariateStreamElement(const Point &pt);
    Point GetMean() override;
  private:
    Point coordinates_;
};

#endif //KERDEP_UNIVARIATESTREAMELEMENT_H
