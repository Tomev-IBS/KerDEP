//
// Created by Tomev on 23/01/2021.
//

#include "kerDepCcWde.h"

#include <cmath>
#include <algorithm>


#include "math_helpers.h"

KerDEP_CC_WDE::KerDEP_CC_WDE(const unsigned int &maximal_number_of_empirical_coefficients,
                             const double &weights_modifier_, WaveletDensityEstimator *(*wde_factory_method)(
                             const vector<double> &), const unsigned int &block_size)
    : CompressedCumulativeWaveletDensityEstimator(maximal_number_of_empirical_coefficients, weights_modifier_,
                                                  wde_factory_method, block_size) { }

void KerDEP_CC_WDE::PerformStep(point *pt) {
  if(block_size_ == block.size()){
    block.clear();
  }

  block.push_back((*pt)[0]);

  if(block_size_ == block.size()){
    UpdateEstimator(block);
  }
}

vector<point> KerDEP_CC_WDE::GetErrorDomain() const {
  vector<point> error_domain = {};

  auto standard_deviation = stDev(block);
  double minimal_value = *std::min_element(block.begin(), block.end()) - 5 * standard_deviation;
  double maximal_value = *std::max_element(block.begin(), block.end()) + 5 * standard_deviation;

  double step_size = (maximal_value - minimal_value) / 1000;

  for(auto current_value = minimal_value; current_value < maximal_value; current_value += step_size){
    error_domain.push_back({current_value});
  }

  return error_domain;
}

std::vector<double> KerDEP_CC_WDE::GetEstimatorValuesOnDomain(std::vector<point> domain) const {
  vector<double> estimator_values_on_domain = {};

  for(const point& pt : domain){
    double value_on_pt = GetValue(pt);
    estimator_values_on_domain.push_back(value_on_pt);
  }

  return estimator_values_on_domain;
}

unsigned int KerDEP_CC_WDE::GetCurrentCoefficientsNumber() const {
  unsigned int current_coefficients_number = 0;

  if(estimators_.empty()){
    return current_coefficients_number;
  }

  for(auto estimator : estimators_){
    current_coefficients_number += estimator->GetEmpiricalCoefficientsNumber();
  }

  return current_coefficients_number;
}

