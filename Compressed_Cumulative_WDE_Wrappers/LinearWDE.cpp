//
// Created by Tomev on 19/01/2021.
//
// Implementation of Linear Wavelet Density Estimator strategy.
//

#include "Compressed_Cumulative_WDE_Wrappers/LinearWDE.h"

#include <cmath>
#include <algorithm>

#include "math_helpers.h"

using std::cout, std::endl;

TranslatedDilatedScalingFunction LinearWDE::translated_dilated_scaling_function_ = TranslatedDilatedScalingFunction(0, 0);

LinearWDE::LinearWDE(const double &threshold)
    :  coefficient_threshold_(threshold) { }

LinearWDE::LinearWDE(vector<EmpiricalCoefficientData> empirical_scaling_coefficients, const double &threshold)
    :  coefficient_threshold_(threshold) , empirical_scaling_coefficients_(empirical_scaling_coefficients) {
  if(!empirical_scaling_coefficients.empty()){
    scaling_function_resolution_index_ = empirical_scaling_coefficients[0].j_;
    k_min_ = empirical_scaling_coefficients[0].k_;
    k_max_ = empirical_scaling_coefficients[empirical_scaling_coefficients.size() - 1].k_;
  }
}

/** Updates linear WDE data.
 * @brief Updates linear WDE data.
 * @param values_block - Block of new values which will be used during update process.
 */
void LinearWDE::UpdateWDEData(const vector<double> &values_block) {
  auto prepared_block = PrepareBlockData(values_block);
  ComputeOptimalResolutionIndex(prepared_block);
  ComputeTranslations(prepared_block);
  ComputeEmpiricalScalingCoefficients(prepared_block);
}

vector<StreamElementData> LinearWDE::PrepareBlockData(const vector<double> &values_block) {
  auto prepared_values_block = vector<StreamElementData>();

  for(auto val : values_block){
    prepared_values_block.emplace_back(StreamElementData());
    prepared_values_block.back().value_ = val;
  }

  std::sort(prepared_values_block.begin(), prepared_values_block.end());

  return prepared_values_block;
}

/** Computes optimal resolition index.
 *
 * @brief Computes optimal resolition index.
 * @param values_block - Values from which resolution (dilation) index is computed.
 * @return Optimal resolution index.
 */
void LinearWDE::ComputeOptimalResolutionIndex(const vector<StreamElementData> &values_block) {
  std::vector<double> unweighted_values = {};

  for(auto val : values_block){
    unweighted_values.push_back(val.value_);
  }

  scaling_function_resolution_index_ = log2(values_block.size()) / 3.0 - 2.0 - log2(StDev(unweighted_values));
}

/** Computes k_min and k_max (translation indices).
 * @brief Computes k_min and k_max (translation indices).
 * @param values_block - SORTED vector of values.
 */
void LinearWDE::ComputeTranslations(const vector<StreamElementData> &values_block) {

  auto support = translated_dilated_scaling_function_.GetOriginalScalingFunctionSupport();
  int support_min = support.first;
  int support_max = support.second;

  k_min_ = ceil(pow(2, scaling_function_resolution_index_) * values_block[0].value_ - support_max);
  k_max_ = floor(pow(2, scaling_function_resolution_index_) * values_block[values_block.size() - 1].value_ - support_min);
}

/** Computes most important (according to threshold) empirical scaling function coefficients.
 * @brief Computes most important (according to threshold) empirical scaling function coefficients.
 * @param values - Vector of values in the block.
 */
void LinearWDE::ComputeEmpiricalScalingCoefficients(const vector<StreamElementData> &values) {

  empirical_scaling_coefficients_ = {};

  if(values.empty()){
    return;
  }

  for(int k = k_min_; k <= k_max_; ++k){

    double coefficient = 0;

    translated_dilated_scaling_function_.UpdateIndices(scaling_function_resolution_index_, k);

    for(auto val : values){
      coefficient += translated_dilated_scaling_function_.GetValue(val.value_);
    }

    coefficient /= values.size();

    EmpiricalCoefficientData data;
    data.coefficient_ = coefficient;
    data.j_ = scaling_function_resolution_index_;
    data.k_ = k;
    empirical_scaling_coefficients_.push_back(data);
  }

}

/** Computes value of decomposed function in 1D point.
 * @brief Computes value of decomposed function in 1D point.
 * @param x - 1D point in which the value of function should be computed.
 * @return Value of DWT in x.
 */
double LinearWDE::GetValue(const double &x) const {
  double result = 0;

  for(auto data : empirical_scaling_coefficients_){
    translated_dilated_scaling_function_.UpdateIndices(data.j_, data.k_);
    result += data.coefficient_ * translated_dilated_scaling_function_.GetValue(x);
  }

  return weight_ * result;
}

/** Computes coefficients for the lower resolution.
 * @brief Computes coefficients for the lower resolution.
 */
void LinearWDE::LowerCoefficientsResolution() {
  auto filter_coefficients = translated_dilated_scaling_function_.GetFilterCoefficients();

  auto lower_resolution_translations = ComputeLowerResolutionTranslations(filter_coefficients.size() - 1,
                                                                          empirical_scaling_coefficients_);

  k_min_ = lower_resolution_translations[0];
  k_max_ = lower_resolution_translations[1];

  empirical_scaling_coefficients_ = ComputeLowerResolutionScalingCoefficients(empirical_scaling_coefficients_);

  --scaling_function_resolution_index_;
}

vector<int> LinearWDE::ComputeLowerResolutionTranslations(const int &number_of_filter_coefficients,
                                                          const vector<EmpiricalCoefficientData> &scaling_coefficients)
                                                          const {
  int min_coefficient_number = 0;
  int max_coefficient_number = number_of_filter_coefficients - 1;

  int k_min = ceil((scaling_coefficients[0].k_ - max_coefficient_number) / 2);
  int k_max = scaling_coefficients[scaling_coefficients.size() - 1].k_;
  k_max = floor((k_max - min_coefficient_number) / 2);

  return {k_min, k_max};
}

int LinearWDE::GetResolutionIndex() const {
  return scaling_function_resolution_index_;
}

vector<EmpiricalCoefficientData> LinearWDE::GetEmpiricalScalingCoefficients() const {
  return empirical_scaling_coefficients_;
}

double LinearWDE::GetWeight() const {
  return weight_;
}

void LinearWDE::SetWeight(const double &new_weight){
  weight_ = new_weight;
}

void LinearWDE::MultiplyWeight(const double &multiplier) {
  weight_ *= multiplier;
}

unsigned int LinearWDE::GetEmpiricalCoefficientsNumber() const {
  return GetEmpiricalScalingCoefficientsNumber() + GetEmpiricalWaveletCoefficientsNumber();
}

WaveletDensityEstimator *LinearWDE::Merge(WaveletDensityEstimator *other_wde) const {

  vector<EmpiricalCoefficientData> merged_coefficients = {};

  for(auto coefficient_data : empirical_scaling_coefficients_){
    merged_coefficients.push_back(coefficient_data);
  }

  for(int i = 0; i < merged_coefficients.size(); ++i) {
    auto weighted_coefficient = merged_coefficients[i].coefficient_ * weight_;
    merged_coefficients[i].coefficient_ = weighted_coefficient;
  }

  auto other_coefficients = other_wde->GetEmpiricalScalingCoefficients();
  auto other_weight = other_wde->GetWeight();

  unsigned int i = 0;

  for(auto coefficient_data : other_coefficients){
    while(merged_coefficients[i].k_ < coefficient_data.k_){
      ++i;
    }
    if(merged_coefficients[i].k_ == coefficient_data.k_){
      merged_coefficients[i].coefficient_ += coefficient_data.coefficient_ * other_weight;
    } else {
      if(merged_coefficients.size() <= i){
        merged_coefficients.push_back(coefficient_data);
      } else {
        merged_coefficients.insert(merged_coefficients.begin() + i, coefficient_data);
      }
      merged_coefficients[i].coefficient_ *= other_weight;
    }
  }

  return new LinearWDE(merged_coefficients);
}

vector<EmpiricalCoefficientData> LinearWDE::GetEmpiricalWaveletCoefficients() const {
  return vector<EmpiricalCoefficientData>();
}

unsigned int LinearWDE::GetEmpiricalScalingCoefficientsNumber() const {
  return empirical_scaling_coefficients_.size();
}

unsigned int LinearWDE::GetEmpiricalWaveletCoefficientsNumber() const {
  return 0;
}

void LinearWDE::RemoveSmallestWaveletCoefficients(const unsigned int &number_of_coefficients) {
  return;
}

vector<EmpiricalCoefficientData> LinearWDE::ComputeLowerResolutionScalingCoefficients(
    const vector<EmpiricalCoefficientData> &coefficients) const {

  // All coefficients should have same resolution.
  int new_resolution = coefficients[0].j_ - 1;

  auto filter_coefficients = translated_dilated_scaling_function_.GetFilterCoefficients();

  auto lower_resolution_translations = ComputeLowerResolutionTranslations(filter_coefficients.size() - 1,
                                                                          coefficients);

  vector<EmpiricalCoefficientData> LowerResolutionEmpiricalCoefficients = {};

  for(int k = lower_resolution_translations[0]; k <= lower_resolution_translations[1]; ++k){

    double lower_res_coefficient = 0;

    for(auto val : coefficients){
      int l = val.k_;
      int filter_coefficient_index = l - 2 * k;

      if(filter_coefficient_index < 0 || filter_coefficient_index >= filter_coefficients.size()){
        continue;
      }

      lower_res_coefficient += val.coefficient_ * filter_coefficients[filter_coefficient_index];
    }

    auto data = EmpiricalCoefficientData();
    data.coefficient_ = lower_res_coefficient;
    data.j_ = new_resolution;
    data.k_ = k;
    LowerResolutionEmpiricalCoefficients.push_back(data);
  }

  return LowerResolutionEmpiricalCoefficients;
}






