//
// Created by Tomev on 26/01/2021.
//

#include <algorithm>

#include "weightedLinearWde.h"
#include "math_helpers.h"

WeightedLinearWDE::WeightedLinearWDE(vector<EmpiricalCoefficientData> empirical_scaling_coefficients,
                                     const double &threshold) : LinearWDE(empirical_scaling_coefficients, threshold) {
  weight_ = 0.99;
}

WeightedLinearWDE::WeightedLinearWDE(const double &threshold) : LinearWDE(threshold) {
  weight_ = 0.99;
};


vector<StreamElementData> WeightedLinearWDE::PrepareBlockData(const vector<double> &values_block) {
  if(values_block.size() != elements_weights_.size()){
    ComputeElementsWeights(values_block.size());
  }

  auto prepared_block_data = vector<StreamElementData>();

  for(unsigned int i = 0; i < values_block.size(); ++i){
    prepared_block_data.emplace_back(StreamElementData());
    prepared_block_data.back().value_ = values_block[i];
    prepared_block_data.back().weight_ = elements_weights_[i];
  }

  std::sort(prepared_block_data.begin(), prepared_block_data.end());

  return prepared_block_data;
}

void WeightedLinearWDE::ComputeOptimalResolutionIndex(const vector<StreamElementData> &values_block) {
  std::vector<double> unweighted_values = {};
  std::vector<double> weights = {};

  for(auto val : values_block){
    unweighted_values.push_back(val.value_);
    weights.push_back(val.weight_);
  }

  resolution_index_ = log2(values_block.size()) / 3.0 - 2.0 - log2(WeightedStDev(unweighted_values, weights));
}

void WeightedLinearWDE::ComputeEmpiricalScalingCoefficients(const vector<StreamElementData> &values) {

  if(values.size() != elements_weights_.size()) {
    ComputeElementsWeights(values.size());
  }

  empirical_scaling_coefficients_ = {};

  if(values.empty()){
    return;
  }

  for(int k = k_min_; k <= k_max_; ++k){

    double coefficient = 0;

    translated_dilated_scaling_function_.UpdateIndices(resolution_index_, k);

    for(auto val : values){
      coefficient += translated_dilated_scaling_function_.GetValue(val.value_) * val.weight_;
    }

    coefficient /= weights_sum_;

    EmpiricalCoefficientData data;
    data.coefficient_ = coefficient;
    data.j_ = resolution_index_;
    data.k_ = k;
    empirical_scaling_coefficients_.push_back(data);
  }

}

void WeightedLinearWDE::ComputeElementsWeights(const unsigned int &elements_number) {
  weights_sum_ = 0;
  elements_weights_ = {};

  for(unsigned int i = 0; i < elements_number; ++i){
    double element_weight = pow(weight_, i);
    elements_weights_.push_back(element_weight);
    weights_sum_ += element_weight;
  }

  std::reverse(elements_weights_.begin(), elements_weights_.end());
}

void WeightedLinearWDE::SetWeight(const double &new_weight) {
  LinearWDE::SetWeight(new_weight);
  ComputeElementsWeights(elements_weights_.size());
}







