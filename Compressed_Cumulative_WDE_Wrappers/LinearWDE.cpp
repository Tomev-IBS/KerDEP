//
// Created by Tomev on 19/01/2021.
//
// Implementation of Linear Wavelet Density Estimator strategy.
//

#include "Comressed_Cumulative_WDE_Wrappers/LinearWDE.h"

#include <cmath>
#include <numeric>
#include <iostream> // For debug

using std::inner_product;
using std::accumulate;

using std::cout, std::endl; // For debug

// TR TODO: Make separate file if it's used elsewhere.
/** Computes standard deviation estimator from given vector.
 * @brief Computes standard deviation estimator from given vector.
 * @param values - Vector of values.
 * @return Standard deviation estimator of values.
 */
double stDev(const vector<double> &values){

  if(values.size() < 1){
    return 0;
  }

  auto mean = accumulate(values.begin(), values.end(), 0.0) / values.size();
  auto squares_sum = inner_product(values.begin(), values.end(), values.begin(), 0.0);

  return sqrt(squares_sum / values.size() - pow(mean, 2));
}

LinearWDE::LinearWDE(const double &threshold)
    :  coefficient_threshold_(threshold) { }

/** Updates linear WDE data.
 * @brief Updates linear WDE data.
 * @param values_block - Block of new values which will be used during update process.
 */
void LinearWDE::UpdateWDEData(const vector<double> &values_block) {
  ComputeOptimalResolutionIndex(values_block);
  ComputeTranslations(values_block);
  ComputeEmpiricalScalingCoefficients(values_block);
}

/** Computes optimal resolition index.
 *
 * @brief Computes optimal resolition index.
 * @param values_block - Values from which resolution (dilation) index is computed.
 * @return Optimal resolution index.
 */
void LinearWDE::ComputeOptimalResolutionIndex(const vector<double> &values_block) {
    resolution_index_ = log2(values_block.size()) / 3.0 - 2.0 - log2(stDev(values_block));
}

/** Computes k_min and k_max (translation indices).
 * @brief Computes k_min and k_max (translation indices).
 * @param values_block - SORTED vector of values.
 */
void LinearWDE::ComputeTranslations(const vector<double> &values_block) {

  auto translated_dilated_scaling_function = TranslatedDilatedScalingFunction(0, 0);
  auto support = translated_dilated_scaling_function.GetOriginalScalingFunctionSupport();
  int support_min = support.first;
  int support_max = support.second;

  k_min_ = ceil(pow(2, resolution_index_) * values_block[0] - support_max);
  k_max_ = floor(pow(2, resolution_index_) * values_block[values_block.size() - 1] - support_min);
}

/** Computes most important (according to threshold) empirical scaling function coefficients.
 * @brief Computes most important (according to threshold) empirical scaling function coefficients.
 * @param values - Vector of values in the block.
 */
void LinearWDE::ComputeEmpiricalScalingCoefficients(const vector<double> &values) {

  empirical_scaling_coefficients_ = {};

  if(values.empty()){
    return;
  }

  auto scaling_function = TranslatedDilatedScalingFunction(resolution_index_, 0);

  for(int k = k_min_; k <= k_max_; ++k){

    cout << "k = " << k << " from [" << k_min_ << ", " << k_max_ << "]\n";

    double coefficient = 0;

    scaling_function.UpdateIndices(resolution_index_, k);

    for(auto val : values){
      coefficient += scaling_function.GetValue(val);
    }

    coefficient /= values.size();

    if(fabs(coefficient) > coefficient_threshold_){
      EmpiricalCoefficientData data;
      data.coefficient_ = coefficient;
      data.j_ = resolution_index_;
      data.k_ = k;
      empirical_scaling_coefficients_.push_back(data);
    }
  }

  cout << "Got " << empirical_scaling_coefficients_.size() << " coefficients.\n";
}

/** Computes value of decomposed function in 1D point.
 * @brief Computes value of decomposed function in 1D point.
 * @param x - 1D point in which the value of function should be computed.
 * @return Value of DWT in x.
 */
double LinearWDE::GetValue(const double &x) const {
  double result = 0;

  auto scaling_function = TranslatedDilatedScalingFunction(resolution_index_, 0);

  for(auto data : empirical_scaling_coefficients_){
    scaling_function.UpdateIndices(data.j_, data.k_);
    result += data.coefficient_ * scaling_function.GetValue(x);
  }

  return weight_ * result;
}

/** Computes coefficients for the lower resolution.
 * @brief Computes coefficients for the lower resolution.
 */
void LinearWDE::LowerCoefficientsResolution() {
  --resolution_index_;

  auto scaling_function = TranslatedDilatedScalingFunction(0, 0);
  auto filter_coefficients = scaling_function.GetFilterCoefficients();

  vector<EmpiricalCoefficientData> LowerResolutionEmpiricalCoefficients = {};

  ComputeLowerResolutionTranslations(filter_coefficients.size() - 1);

  for(int k = k_min_; k <= k_max_; ++k){

    double lower_res_coefficient = 0;

    for(auto val : empirical_scaling_coefficients_){
      int l = val.k_;
      int filter_coefficient_index = l - 2 * k;

      if(filter_coefficient_index < 0 || filter_coefficient_index >= filter_coefficients.size()){
        continue;
      }

      lower_res_coefficient += val.coefficient_ * filter_coefficients[filter_coefficient_index];
    }

    cout << "Coefficient for k = " << k << " is " << lower_res_coefficient << ".\n";

    if(fabs(lower_res_coefficient) > coefficient_threshold_){
      auto data = EmpiricalCoefficientData();
      data.coefficient_ = lower_res_coefficient;
      data.j_ = resolution_index_;
      data.k_ = k;
      LowerResolutionEmpiricalCoefficients.push_back(data);
    }

  }

  empirical_scaling_coefficients_ = LowerResolutionEmpiricalCoefficients;

  cout << "Got " << empirical_scaling_coefficients_.size() << " coefficients.\n";
}

void LinearWDE::ComputeLowerResolutionTranslations(const int &number_of_filter_coefficients) {
  int min_coefficient_number = 0;
  int max_coefficient_number = number_of_filter_coefficients - 1;

  int k_min = ceil((empirical_scaling_coefficients_[0].k_ - max_coefficient_number) / 2);
  int k_max = empirical_scaling_coefficients_[empirical_scaling_coefficients_.size() - 1].k_;
  k_max = floor((k_max - min_coefficient_number) / 2);

  k_min_ = k_min;
  k_max_ = k_max;
}

int LinearWDE::GetResolutionIndex() const {
  return resolution_index_;
}

void LinearWDE::MultiplyWeight(const double &multiplier) {
  weight_ *= multiplier;
}

double LinearWDE::GetWeight() const {
  return weight_;
}

vector<EmpiricalCoefficientData> LinearWDE::GetEmpiricalCoefficients() const {
  return empirical_scaling_coefficients_;
}
