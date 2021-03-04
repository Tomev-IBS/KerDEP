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

struct StreamElementData {
    double value_ = 0;
    double weight_ = 1;

  bool operator >(const StreamElementData &other_stream_element) const {
    return this->value_ > other_stream_element.value_;
  }

  bool operator <(const StreamElementData &other_stream_element) const {
    return this->value_ < other_stream_element.value_;
  }

};

class LinearWDE : public WaveletDensityEstimator {
  public:
    explicit LinearWDE(const double &threshold = 1e-5);
    LinearWDE(vector<EmpiricalCoefficientData> empirical_scaling_coefficients,
              const double &threshold = 1e-5);

    double GetValue(const double &x) const override;

    void UpdateWDEData(const vector<double> &values_block) override;
    void LowerCoefficientsResolution() override;
    void RemoveSmallestWaveletCoefficients(const unsigned int &number_of_coefficients) override;

    int GetResolutionIndex() const override;
    vector<EmpiricalCoefficientData> GetEmpiricalScalingCoefficients() const override;
    vector<EmpiricalCoefficientData> GetEmpiricalWaveletCoefficients() const override;
    unsigned int GetEmpiricalScalingCoefficientsNumber() const override;
    unsigned int GetEmpiricalWaveletCoefficientsNumber() const override;
    unsigned int GetEmpiricalCoefficientsNumber() const override;

    double GetWeight() const override;
    void SetWeight(const double &new_weight) override;
    void MultiplyWeight(const double &multiplier) override;

    WaveletDensityEstimator* Merge(WaveletDensityEstimator *other_wde) const override;

  protected:

    static TranslatedDilatedScalingFunction translated_dilated_scaling_function_;

    int scaling_function_resolution_index_ = 0;
    int k_min_ = 0;
    int k_max_ = 0;
    vector<EmpiricalCoefficientData> empirical_scaling_coefficients_ = {};
    double weight_ = 1;

    double coefficient_threshold_ = 1e-5; // Denotes if the coefficient should be included.

    virtual vector<StreamElementData> PrepareBlockData(const vector<double> &values_block);
    virtual void ComputeOptimalResolutionIndex(const vector<StreamElementData> &values_block);
    virtual void ComputeTranslations(const vector<StreamElementData> &values_block);
    virtual void ComputeEmpiricalScalingCoefficients(const vector<StreamElementData> &values);
    vector<int> ComputeLowerResolutionTranslations(const int &number_of_filter_coefficients,
                                                   const vector<EmpiricalCoefficientData> &scaling_coefficients) const;
    vector<EmpiricalCoefficientData> ComputeLowerResolutionScalingCoefficients(
        const vector<EmpiricalCoefficientData> &coefficients) const;
};

#endif //COMPRESSED_CUMULATIVE_WDE_OVER_STREAM_LINEARWDE_H
