//
// Created by Tomev on 19/01/2021.
//
// Header file for Linear Wavelet Density Estimation Strategy. It's meant for 1D data.
//

#ifndef COMPRESSED_CUMULATIVE_WDE_OVER_STREAM_LINEARWDE_H
#define COMPRESSED_CUMULATIVE_WDE_OVER_STREAM_LINEARWDE_H

#include "WaveletDensityEstimator.h"
#include "TranslatedDilatedScalingFunction.h"

#include <vector>

using std::vector;

class LinearWDE : public WaveletDensityEstimator {
  public:
    explicit LinearWDE(const double &threshold = 1e-5);

    double GetValue(const double &x) const override;

    void UpdateWDEData(const vector<double> &values_block) override;
    void LowerCoefficientsResolution() override;

    int GetResolutionIndex() const override;
    vector<EmpiricalCoefficientData> GetEmpiricalCoefficients() const override;

    double GetWeight() const override;
    void SetWeight(const double &new_weight) override;
    void MultiplyWeight(const double &multiplier) override;


  protected:

    int resolution_index_ = 0;
    int k_min_ = 0;
    int k_max_ = 0;
    vector<EmpiricalCoefficientData> empirical_scaling_coefficients_ = {};
    double weight_ = 1;

    double coefficient_threshold_ = 1e-5; // Denotes if the coefficient should be included.

    void ComputeOptimalResolutionIndex(const vector<double> &values_block);
    void ComputeTranslations(const vector<double> &values_block);
    void ComputeEmpiricalScalingCoefficients(const vector<double> &values);
    void ComputeLowerResolutionTranslations(const int &number_of_filter_coefficients);
};

#endif //COMPRESSED_CUMULATIVE_WDE_OVER_STREAM_LINEARWDE_H
