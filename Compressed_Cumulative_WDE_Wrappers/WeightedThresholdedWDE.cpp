//
// Created by Tomev on 03/02/2021.
//

#include "WeightedThresholdedWDE.h"

TranslatedDilatedWaveletFunction WeightedThresholdedWDE::translated_dilated_wavelet_function_ = TranslatedDilatedWaveletFunction(0, 0);

WeightedThresholdedWDE::WeightedThresholdedWDE(ThresholdingStrategyPtr thresholding_strategy, const double &threshold, const double &weight)
{
  thresholding_strategy_ = thresholding_strategy;
  weight_ = 0.95;
}

void WeightedThresholdedWDE::UpdateWDEData(const vector<double> &values_block) {
  auto prepared_block = PrepareBlockData(values_block);
  WeightedLinearWDE::UpdateWDEData(values_block);
  ComputeEmpiricalWaveletCoefficients(prepared_block);
  thresholding_strategy_->UpdateThresholds(empirical_wavelet_coefficients_);
}

void WeightedThresholdedWDE::ComputeEmpiricalWaveletCoefficients(const vector<StreamElementData> &values) {

  empirical_wavelet_coefficients_ = {};

  if(values.empty()){
    return;
  }

  for(int j = scaling_function_resolution_index_; j <= resolution_index_1_; ++j) {
    for(int k = k_min_; k <= k_max_; ++k) {

      double coefficient = 0;

      translated_dilated_wavelet_function_.UpdateIndices(j, k);

      for(auto val : values) {
        coefficient += translated_dilated_wavelet_function_.GetValue(val.value_) * val.weight_;
      }

      coefficient /= weights_sum_;

      EmpiricalCoefficientData data;
      data.coefficient_ = coefficient;
      data.j_ = j;
      data.k_ = k;
      empirical_wavelet_coefficients_.push_back(data);
    }
  }
}

void WeightedThresholdedWDE::RemoveSmallestWaveletCoefficients(const unsigned int &number_of_coefficients) {

  if(number_of_coefficients == 0 || number_of_coefficients >= empirical_wavelet_coefficients_.size()){
    empirical_wavelet_coefficients_.clear();
    return;
  }

  std::sort(empirical_wavelet_coefficients_.begin(), empirical_wavelet_coefficients_.end());
  for(unsigned int i = 0; i < number_of_coefficients; ++i){
    empirical_wavelet_coefficients_.erase(empirical_wavelet_coefficients_.begin());
  }
}

void WeightedThresholdedWDE::LowerCoefficientsResolution() {

  auto filter_coefficients = translated_dilated_scaling_function_.GetFilterCoefficients();

  auto lower_resolution_translations = ComputeLowerResolutionTranslations(filter_coefficients.size() - 1,
                                                                          empirical_scaling_coefficients_);

  k_min_ = lower_resolution_translations[0];
  k_max_ = lower_resolution_translations[1];

  empirical_wavelet_coefficients_ = ComputeLowerResolutionWaveletCoefficients(empirical_scaling_coefficients_);
  empirical_scaling_coefficients_ = ComputeLowerResolutionScalingCoefficients(empirical_scaling_coefficients_);

  scaling_function_resolution_index_ = resolution_index_1_ = scaling_function_resolution_index_ - 1;
}

unsigned int WeightedThresholdedWDE::GetEmpiricalWaveletCoefficientsNumber() const {
  return empirical_wavelet_coefficients_.size();
}

void WeightedThresholdedWDE::ComputeOptimalResolutionIndex(const vector<StreamElementData> &values_block) {

  double r = translated_dilated_scaling_function_.GetRegularity();
  int n = values_block.size();

  scaling_function_resolution_index_ = ceil(log2(n) / (2 * r + 1));
  resolution_index_1_ = floor(log2(n) - log2(log(n)));
}

vector<EmpiricalCoefficientData> WeightedThresholdedWDE::ComputeLowerResolutionWaveletCoefficients(
    const vector<EmpiricalCoefficientData> &coefficients) const {
  // All coefficients should have same resolution.
  int new_resolution = coefficients[0].j_ - 1;

  auto filter_coefficients = translated_dilated_scaling_function_.GetFilterCoefficients();

  vector<EmpiricalCoefficientData> LowerResolutionEmpiricalCoefficients = {};

  for(int k = k_min_; k <= k_max_; ++k){

    double lower_res_coefficient = 0;

    for(auto val : coefficients){
      int l = val.k_;
      int filter_coefficient_index = -l + 2 * k + 1;

      if(filter_coefficient_index < 0 || filter_coefficient_index >= filter_coefficients.size()){
        continue;
      }

      lower_res_coefficient += pow(-1, l) * val.coefficient_ * filter_coefficients[filter_coefficient_index];
    }

    auto data = EmpiricalCoefficientData();
    data.coefficient_ = lower_res_coefficient;
    data.j_ = new_resolution;
    data.k_ = k;
    LowerResolutionEmpiricalCoefficients.push_back(data);
  }

  return LowerResolutionEmpiricalCoefficients;
}

double WeightedThresholdedWDE::GetValue(const double &x) const {

  double val = WeightedLinearWDE::GetValue(x);

  val += ComputeWaveletAddend(x);

  return val;
}

double WeightedThresholdedWDE::ComputeWaveletAddend(const double &x) const {

  double wavelet_addend = 0;

  for(auto val : empirical_wavelet_coefficients_){
    translated_dilated_wavelet_function_.UpdateIndices(val.j_, val.k_);
    wavelet_addend += thresholding_strategy_->ComputeThresholdedCoefficient(val) * translated_dilated_wavelet_function_.GetValue(x);
  }

  return wavelet_addend;
}