//
// Created by Tomev on 23/01/2021.
//

#ifndef DEDSTA_DEDSTACCWDE_H
#define DEDSTA_DEDSTACCWDE_H

#include "CompressedCumulativeWaveletDensityEstimator.h"

class DEDSTA_CC_WDE : protected CompressedCumulativeWaveletDensityEstimator {
  public:
    DEDSTA_CC_WDE(const unsigned int &maximal_number_of_empirical_coefficients,
                  const double &weights_modifier_,
                  WaveletDensityEstimator* (*wde_factory_method)(const vector<double> &values_block),
                  const unsigned int &block_size);

    virtual void PerformStep(point *pt);
    vector<point> GetErrorDomain() const;
    vector<point> GetErrorDomainFromBlock() const;
    vector<point> GetErrorDomainFromCoefficients() const;
    std::pair<double, double> GetEstimatorSupport() const;
    std::vector<double> GetEstimatorValuesOnDomain(std::vector<point> domain) const;
    unsigned int GetCurrentCoefficientsNumber() const;


  protected:
    vector<double> block;
};

#endif //DEDSTA_DEDSTACCWDE_H
