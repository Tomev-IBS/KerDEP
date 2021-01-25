//
// Created by Tomev on 26/01/2021.
//

#ifndef KERDEP_WEIGHTEDLINEARWDE_H
#define KERDEP_WEIGHTEDLINEARWDE_H

#include "LinearWDE.h"

class WeightedLinearWDE : public LinearWDE {

  public:
    explicit WeightedLinearWDE(const double &threshold = 1e-5);
    WeightedLinearWDE(vector<EmpiricalCoefficientData> empirical_scaling_coefficients,
    const double &threshold = 1e-5);

    void SetWeight(const double &new_weight) override;

  protected:
    virtual void ComputeEmpiricalScalingCoefficients(const vector<double> &values) override;
    void ComputeElementsWeights(const unsigned int &elements_number);

    vector<double> elements_weights_;
    double weights_sum_;

};

#endif //KERDEP_WEIGHTEDLINEARWDE_H
