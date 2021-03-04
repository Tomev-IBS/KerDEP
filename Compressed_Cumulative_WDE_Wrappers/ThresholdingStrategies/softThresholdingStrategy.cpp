//
// Created by Tomev on 08/02/2021.
//

#include "softThresholdingStrategy.h"

double SoftThresholdingStrategy::ComputeThresholdedCoefficient(const EmpiricalCoefficientData &wavelet_coefficient) {

  if(wavelet_coefficient.coefficient_ > thresholds_[wavelet_coefficient.j_]){
    return wavelet_coefficient.coefficient_ - thresholds_[wavelet_coefficient.j_];
  }

  if(wavelet_coefficient.coefficient_ < - thresholds_[wavelet_coefficient.j_]){
    return wavelet_coefficient.coefficient_ + thresholds_[wavelet_coefficient.j_];
  }

  return 0;
}
