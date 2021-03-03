//
// Created by Tomev on 07/02/2021.
//

#include "hardThresholdingStrategy.h"

void HardThresholdingStrategy::UpdateThresholds(const vector<EmpiricalCoefficientData> &wavelet_coefficients) {

  thresholds_.clear();

  // I assume that wavelet coefficients are sorted according to j.
  int j0 = wavelet_coefficients[0].j_;
  int j1 = wavelet_coefficients[wavelet_coefficients.size() - 1].j_;

  int i = 0;
  for(int j = j0; j <= j1; ++j){

    double max_abs_value_of_coefficient = fabs(wavelet_coefficients[i].coefficient_);

    while(wavelet_coefficients[i].j_ == j){
      if(fabs(wavelet_coefficients[i].coefficient_) > max_abs_value_of_coefficient){
        max_abs_value_of_coefficient = fabs(wavelet_coefficients[i].coefficient_);
      }
      ++i;
    }

    thresholds_[j] = threshold_coefficient_ * max_abs_value_of_coefficient;
  }
}

double HardThresholdingStrategy::ComputeThresholdedCoefficient(const EmpiricalCoefficientData &wavelet_coefficient) {
  if(fabs(wavelet_coefficient.coefficient_) > thresholds_[wavelet_coefficient.j_]){
    return wavelet_coefficient.coefficient_;
  }
  return 0;
}
