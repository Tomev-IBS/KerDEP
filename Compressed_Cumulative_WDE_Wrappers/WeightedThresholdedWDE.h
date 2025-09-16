//
// Created by Tomev on 03/02/2021.
//

#ifndef DEDSTA_WEIGHTEDTHRESHOLDEDWDE_H
#define DEDSTA_WEIGHTEDTHRESHOLDEDWDE_H

#include "weightedLinearWde.h"
#include "TranslatedDilatedWaveletFunction.h"

#include "ThresholdingStrategies/ThresholdingStrategyInterface.h"

class WeightedThresholdedWDE : public WeightedLinearWDE {
  public:

    explicit WeightedThresholdedWDE(ThresholdingStrategyPtr thresholding_strategy, const double &weight = 1,
                                    const double &threshold = 1e-5);


    double GetValue(const double &x) const override;
    void UpdateWDEData(const vector<double> &values_block) override;
    void LowerCoefficientsResolution() override;
    unsigned int GetEmpiricalWaveletCoefficientsNumber() const override;
    void RemoveSmallestWaveletCoefficients(const unsigned int &number_of_coefficients) override;

  protected:

    static TranslatedDilatedWaveletFunction translated_dilated_wavelet_function_;
    ThresholdingStrategyPtr thresholding_strategy_;
    vector<EmpiricalCoefficientData> empirical_wavelet_coefficients_ = {};

    int resolution_index_1_ = 0;

    void ComputeOptimalResolutionIndex(const vector<StreamElementData> &values_block) override;
    virtual void ComputeEmpiricalWaveletCoefficients(const vector<StreamElementData> &values);
    double ComputeWaveletAddend(const double &x) const;
    vector<EmpiricalCoefficientData> ComputeLowerResolutionWaveletCoefficients(
        const vector<EmpiricalCoefficientData> &coefficients) const;

};

#endif //DEDSTA_WEIGHTEDTHRESHOLDEDWDE_H
