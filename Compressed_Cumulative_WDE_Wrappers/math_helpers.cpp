//
// Created by Tomev on 23/01/2021.
//

#include "math_helpers.h"

#include <cmath>
#include <numeric>
using std::inner_product;
using std::accumulate;

/** Computes standard deviation estimator from given vector.
 * @brief Computes standard deviation estimator from given vector.
 * @param values - Vector of values.
 * @return Standard deviation estimator of values.
 */
double StDev(const vector<double> &values){

  if(values.size() < 2){
    return 0;
  }

  auto mean = accumulate(values.begin(), values.end(), 0.0) / values.size();
  auto squares_sum = inner_product(values.begin(), values.end(), values.begin(), 0.0);

  return sqrt(squares_sum / values.size() - pow(mean, 2));
}

double WeightedStDev(const vector<double> &values, const vector<double> &weights){

  if(values.size() != weights.size() || values.size() < 2){
    return StDev(values);
  }

  double s0 = accumulate(weights.begin(), weights.end(), 0.0);
  double s1 = 0;
  double s2 = 0;

  for(int i  = 0; i < values.size(); ++i){
    s1 += values[i] * weights[i];
    s2 += values[i] * values[i] * weights[i];
  }

  double stDev = s0 * s2 - s1 * s1;
  return sqrt(stDev) / s0;
}