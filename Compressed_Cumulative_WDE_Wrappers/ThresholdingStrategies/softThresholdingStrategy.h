//
// Created by Tomev on 08/02/2021.
//

#ifndef DEDSTA_SOFTTHRESHOLDINGSTRATEGY_H
#define DEDSTA_SOFTTHRESHOLDINGSTRATEGY_H

#include "ThresholdingStrategies/hardThresholdingStrategy.h"

class SoftThresholdingStrategy : public HardThresholdingStrategy {
  public:
    SoftThresholdingStrategy() = default;
    double ComputeThresholdedCoefficient(const EmpiricalCoefficientData &wavelet_coefficient) override;
  protected:
    double threshold_coefficient_ = 0.4; // Hardle's book, page 146, chap 10.3
};

#endif //DEDSTA_SOFTTHRESHOLDINGSTRATEGY_H
