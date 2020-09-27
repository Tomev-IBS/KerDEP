#ifndef KERDEP_ERRORSCALCULATOR_H
#define KERDEP_ERRORSCALCULATOR_H

#include <vector>

class ErrorsCalculator {

    /** Error calculator class is designed to calculate the usual notions of describing differences between two
     * probability density functions. It includes methods for L1, L2, sup and mod error calculations.
     *
     * @brief Error calculator class is designed to calculate the usual notions of describing differences between two
     * probability density functions.
     */
  public:
    ErrorsCalculator(std::vector<double> *model_function_values,
                     std::vector<double> *estimator_values,
                     std::vector<std::vector<double>> *error_domain,
                     double *domain_quantity);
    double CalculateL1Error();
    double CalculateL2Error();
    double CalculateSupError();
    double CalculateModError();
  private:
    std::vector<double> *model_function_values_;
    std::vector<double> *estimator_values_;
    std::vector<std::vector<double>> *error_domain_;
    double *domain_quantity_; // For one dimension it should be length, in two it should be area, and so...

    double FindMaxValueIndex(const std::vector<double> &values);
    double CalculateEuclideanDistance(const std::vector<double> &point_1, const std::vector<double> &point_2);
};

#endif //KERDEP_ERRORSCALCULATOR_H
