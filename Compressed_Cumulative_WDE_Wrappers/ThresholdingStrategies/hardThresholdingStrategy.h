//
// Created by Tomev on 07/02/2021.
//

#ifndef KERDEP_HARDTHRESHOLDINGSTRATEGY_H
#define KERDEP_HARDTHRESHOLDINGSTRATEGY_H

#include "ThresholdingStrategyInterface.h"

#include <unordered_map>

class HardThresholdingStrategy : public ThresholdingStrategyInterface {
  public:
    HardThresholdingStrategy() = default;
    void UpdateThresholds(const vector<EmpiricalCoefficientData> &wavelet_coefficients) override;
    double ComputeThresholdedCoefficient(const EmpiricalCoefficientData &wavelet_coefficient) override;
  protected:
    std::unordered_map<int, double> thresholds_;
    double threshold_coefficient_ = 0.6; // Hardle's book, page 146, chap 10.3
};

#endif //KERDEP_HARDTHRESHOLDINGSTRATEGY_H
