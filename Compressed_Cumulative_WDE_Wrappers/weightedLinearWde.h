//
// Created by Tomev on 26/01/2021.
//

#ifndef DEDSTA_WEIGHTEDLINEARWDE_H
#define DEDSTA_WEIGHTEDLINEARWDE_H

#include "LinearWDE.h"

class WeightedLinearWDE : public LinearWDE {

  public:
    explicit WeightedLinearWDE(const double &threshold = 1e-5);
    WeightedLinearWDE(vector<EmpiricalCoefficientData> empirical_scaling_coefficients,
    const double &threshold = 1e-5);

    void SetWeight(const double &new_weight) override;

  protected:

    virtual vector<StreamElementData> PrepareBlockData(const vector<double> &values_block);
    virtual void ComputeOptimalResolutionIndex(const vector<StreamElementData> &values_block);
    void ComputeEmpiricalScalingCoefficients(const vector<StreamElementData> &values) override;
    void ComputeElementsWeights(const unsigned int &elements_number);

    vector<double> elements_weights_;
    double weights_sum_;

};

#endif //DEDSTA_WEIGHTEDLINEARWDE_H
