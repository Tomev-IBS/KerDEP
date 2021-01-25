//
// Created by Tomev on 26/01/2021.
//

#include "weightedLinearWde.h"

WeightedLinearWDE::WeightedLinearWDE(vector<EmpiricalCoefficientData> empirical_scaling_coefficients,
                                     const double &threshold) : LinearWDE(empirical_scaling_coefficients, threshold) {
  weight_ = 0.99;
}

void WeightedLinearWDE::ComputeEmpiricalScalingCoefficients(const vector<double> &values) {

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

    for(unsigned int i = 0; i < values.size(); ++i){
      auto val = values[i];
      auto weight = elements_weights_[i];
      coefficient += translated_dilated_scaling_function_.GetValue(val) * weight;
    }

    coefficient /= weights_sum_;

    //if(fabs(coefficient) > coefficient_threshold_){
    EmpiricalCoefficientData data;
    data.coefficient_ = coefficient;
    data.j_ = resolution_index_;
    data.k_ = k;
    empirical_scaling_coefficients_.push_back(data);
    //}
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
}

void WeightedLinearWDE::SetWeight(const double &new_weight) {
  LinearWDE::SetWeight(new_weight);
  ComputeElementsWeights(elements_weights_.size());
}

WeightedLinearWDE::WeightedLinearWDE(const double &threshold) : LinearWDE(threshold) {};


