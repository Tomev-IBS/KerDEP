#include "errorsCalculator.h"

#include <math.h>

double ErrorsCalculator::CalculateL1Error() {
  /** This method calculates l1 error between estimated values and model values.
   * @brief This method calculates l1 error between estimated values and model values.
   */
  // Require that # of model and estimated values are equal.
  if(model_function_values_->size() != estimator_values_->size()){
    return -1;
  }
  // First calculate sum of absolute values of differences between estimations and model on the domain.
  double sum_of_absolute_differences = 0;
  for(auto i = 0; i < model_function_values_->size(); ++i){
    sum_of_absolute_differences += fabs((*model_function_values_)[i] - (*estimator_values_)[i]);
  }
  // Then, get its average.
  auto average_of_absolute_differences = sum_of_absolute_differences / model_function_values_->size();
  // Finally multiply it by domains length, area or so...
  return average_of_absolute_differences * *domain_quantity_;
}

double ErrorsCalculator::CalculateL2Error() {
  /** This method calculates l2 error between estimated values and model values.
   * @brief This method calculates l2 error between estimated values and model values.
   */
  // Require that # of model and estimated values are equal.
  if(model_function_values_->size() != estimator_values_->size()){
    return -1;
  }
  // First calculate sum of squared values of differences between estimations and model on the domain.
  double sum_of_absolute_differences = 0;
  for(auto i = 0; i < model_function_values_->size(); ++i){
    sum_of_absolute_differences += pow((*model_function_values_)[i] - (*estimator_values_)[i], 2);
  }
  // Then, get its average.
  auto average_of_absolute_differences = sum_of_absolute_differences / model_function_values_->size();
  // Finally multiply it by domains length, area or so, and square it all.
  return pow(average_of_absolute_differences * *domain_quantity_, 0.5);
}

double ErrorsCalculator::CalculateSupError() {
  /** This method calculates sup error between estimated values and model values.
   * @brief This method calculates sup error between estimated values and model values.
   */

  // Initialize the sup.
  double sup = fabs((*model_function_values_)[0] - (*estimator_values_)[0]);

  // Find the greatest absolute difference.
  for(int i = 0; i < model_function_values_->size(); ++i) {
    double new_sup = fabs((*model_function_values_)[i] - (*estimator_values_)[i]);
    sup = new_sup > sup ? new_sup : sup;
  }

  // Return it.
  return sup;
}

double ErrorsCalculator::CalculateModError() {
  return 0;
}

double ErrorsCalculator::FindMaxValueIndex(const std::vector<double> &values) {
  /** Finds the index of maximal value in given vector. It checks one value at the time.
   * @brief Finds the index of maximal value in given vector.
   */
   // If there are no values in the vector return -1;

   // Initialize index and maximal value with the first value of vector.
}
