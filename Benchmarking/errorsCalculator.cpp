#include "errorsCalculator.h"

#include <math.h>

ErrorsCalculator::ErrorsCalculator(std::vector<double> *model_function_values,
                                   std::vector<double> *estimator_values,
                                   std::vector<
                                       std::vector<double>> *error_domain,
                                   double *domain_quantity)
 : model_function_values_(model_function_values),
   estimator_values_(estimator_values), error_domain_(error_domain),
   domain_quantity_(domain_quantity)
{}

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
  /** Calculates mod error between model and estimated values.
   * @brief Calculates mod error between model and estimated values.
   */
  // Ensure, that both estimated and model values have the same dimension
  // as domain. In this case return -1.
  if( error_domain_->size() != estimator_values_->size() ||
      error_domain_->size() != model_function_values_->size()){
    return -1;
  }
  // Find greatest value of model and estimated values.
  auto max_model_value_index = FindMaxValueIndex(*model_function_values_);
  auto max_estimator_value_index =
      FindMaxValueIndex(*estimator_values_);
  auto model_max_value_point = (*error_domain_)[max_model_value_index];
  auto estimator_max_value_point = (*error_domain_)[max_estimator_value_index];
  // Return distance between selected points.
  return CalculateEuclideanDistance(model_max_value_point,
                                    estimator_max_value_point);
}

double ErrorsCalculator::FindMaxValueIndex(const std::vector<double> &values) {
  /** Finds the index of maximal value in given vector. It checks one value at the time.
   * @brief Finds the index of maximal value in given vector.
   */
   // If there are no values in the vector return -1;
   if(values.empty()){
     return -1;
   }
   // Initialize index and maximal value with the first value of vector.
   int max_value_index = 0;
   double max_value = values[0];

   // Check if there are no greater values then initialized one.
   for(auto i = 0; i < values.size(); ++i){
     if(max_value < values[i]){
       max_value = values[i];
       max_value_index = i;
     }
   }

   // Return value found.
   return max_value_index;
}

double ErrorsCalculator::CalculateEuclideanDistance(
    const std::vector<double> &point_1, const std::vector<double> &point_2) {
  /** Calculates euclidean distance between two given points. For this method
   * to work correctly both points have to have same dimension.
   * @brief Calculates euclidean distance between two given points.
   */
  // Ensure, that dimension of both points i equal. If not return -1.
  if(point_1.size() != point_2.size()){
    return -1;
  }
  // Calculate sum of squares of differences between the points.
  double sum = 0;
  for(auto i = 0; i < point_2.size(); ++i){
    sum += pow(point_1[i] - point_2[i], 2);
  }
  // Return square root of this sum (the euclidean distance).
  return pow(sum, 0.5);
}


