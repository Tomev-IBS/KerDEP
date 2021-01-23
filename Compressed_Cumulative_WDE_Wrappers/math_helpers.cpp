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
double stDev(const vector<double> &values){

  if(values.size() < 2){
    return 0;
  }

  auto mean = accumulate(values.begin(), values.end(), 0.0) / values.size();
  auto squares_sum = inner_product(values.begin(), values.end(), values.begin(), 0.0);

  return sqrt(squares_sum / values.size() - pow(mean, 2));
}