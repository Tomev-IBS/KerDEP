//
// Created by Tomev on 23/01/2021.
//

#include "DEDSTACcWde.h"
#include "TranslatedDilatedScalingFunction.h"

#include <algorithm>

#include "math_helpers.h"

DEDSTA_CC_WDE::DEDSTA_CC_WDE(const unsigned int &maximal_number_of_empirical_coefficients,
                             const double &weights_modifier_, WaveletDensityEstimator *(*wde_factory_method)(
                             const vector<double> &), const unsigned int &block_size)
    : CompressedCumulativeWaveletDensityEstimator(maximal_number_of_empirical_coefficients, weights_modifier_,
                                                  wde_factory_method, block_size) { }

void DEDSTA_CC_WDE::PerformStep(point *pt) {
  if(block_size_ == block.size()){
    block.clear();
  }

  block.push_back((*pt)[0]);

  if(block_size_ == block.size()){
    UpdateEstimator(block);
  }
}

vector<point> DEDSTA_CC_WDE::GetErrorDomain() const {

  if(estimators_.empty()){
    return GetErrorDomainFromBlock();
  }

  return GetErrorDomainFromCoefficients();
}

std::vector<double> DEDSTA_CC_WDE::GetEstimatorValuesOnDomain(std::vector<point> domain) const {
  vector<double> estimator_values_on_domain = {};

  for(const point& pt : domain){
    double value_on_pt = GetValue(pt);
    estimator_values_on_domain.push_back(value_on_pt);
  }

  return estimator_values_on_domain;
}

unsigned int DEDSTA_CC_WDE::GetCurrentCoefficientsNumber() const {
  unsigned int current_coefficients_number = 0;

  if(estimators_.empty()){
    return current_coefficients_number;
  }

  for(auto estimator : estimators_){
    current_coefficients_number += estimator->GetEmpiricalCoefficientsNumber();
  }

  return current_coefficients_number;
}

std::pair<double, double> DEDSTA_CC_WDE::GetEstimatorSupport() const {

  if(estimators_.empty()){
    return {0, 0};
  }

  auto coefficients = estimators_[0]->GetEmpiricalScalingCoefficients();

  auto phi_jk = TranslatedDilatedScalingFunction(coefficients[0].j_, coefficients[0].k_);
  auto support = phi_jk.GetOriginalScalingFunctionSupport();

  phi_jk.UpdateIndices(coefficients[coefficients.size() - 1].j_, coefficients[coefficients.size() - 1].k_);
  support.second = phi_jk.GetTranslatedDilatedScalingFunctionSupport().second;

  for(unsigned int i = 1; i < estimators_.size(); ++i){

    coefficients = estimators_[i]->GetEmpiricalScalingCoefficients();

    auto coefficient = coefficients[0];
    phi_jk.UpdateIndices(coefficient.j_, coefficient.k_);
    auto current_support = phi_jk.GetTranslatedDilatedScalingFunctionSupport();

    if(support.first > current_support.first){
      support.first = current_support.first;
    }

    coefficient = coefficients[coefficients.size() - 1];
    phi_jk.UpdateIndices(coefficient.j_, coefficient.k_);
    current_support = phi_jk.GetTranslatedDilatedScalingFunctionSupport();

    if(support.second < current_support.second){
      support.second = current_support.second;
    }

  }

  return support;
}

vector<point> DEDSTA_CC_WDE::GetErrorDomainFromBlock() const {

  vector<point> error_domain = {};

  double standard_deviation = StDev(block);
  double minimal_value = *std::min_element(block.begin(), block.end()) - 5 * standard_deviation;
  double maximal_value = *std::max_element(block.begin(), block.end()) + 5 * standard_deviation;

  auto sections_number = 1000;

  double step_size = (maximal_value - minimal_value) / sections_number; // We want hard codded 1000 sections

  for(auto current_value = minimal_value; current_value < maximal_value; current_value += step_size){
    error_domain.push_back({current_value});
  }

  return error_domain;
}

vector<point> DEDSTA_CC_WDE::GetErrorDomainFromCoefficients() const {

  vector<point> error_domain = {};

  auto support = GetEstimatorSupport();

  double minimal_value = support.first;
  double maximal_value = support.second;

  auto sections_number = 1000;

  double step_size = (maximal_value - minimal_value) / sections_number; // We want hard codded 1000 sections

  for(auto current_value = minimal_value; current_value < maximal_value; current_value += step_size){
    error_domain.push_back({current_value});
  }

  return error_domain;
}

