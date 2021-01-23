//
// Created by Tomev on 23/01/2021.
//

#ifndef KERDEP_KERDEPCCWDE_H
#define KERDEP_KERDEPCCWDE_H

#include "CompressedCumulativeWaveletDensityEstimator.h"

class KerDEP_CC_WDE : protected CompressedCumulativeWaveletDensityEstimator {
  public:
    KerDEP_CC_WDE(const unsigned int &maximal_number_of_empirical_coefficients,
                  const double &weights_modifier_,
                  WaveletDensityEstimator* (*wde_factory_method)(const vector<double> &values_block),
                  const unsigned int &block_size);

    void PerformStep(point *pt);
    vector<point> GetErrorDomain() const;
    std::vector<double> GetEstimatorValuesOnDomain(std::vector<point> domain) const;
    unsigned int GetCurrentCoefficientsNumber() const;


  protected:
    vector<double> block;
};

#endif //KERDEP_KERDEPCCWDE_H
