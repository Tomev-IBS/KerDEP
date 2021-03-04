//
// Created by Tomev on 25/01/2021.
//

#include "kerDepWindowedWde.h"

KerDEPWindowedWDE::KerDEPWindowedWDE(const unsigned int &maximal_number_of_empirical_coefficients,
                                     const double &weights_modifier_,
                                     WaveletDensityEstimator *(*wde_factory_method)(const vector<double> &),
                                     const unsigned int &block_size) : KerDEP_CC_WDE(
    maximal_number_of_empirical_coefficients, weights_modifier_, wde_factory_method, block_size) {}

void KerDEPWindowedWDE::PerformStep(point *pt) {

  block.push_back((*pt)[0]);

  if(block_size_ < block.size()){
    block.erase(block.begin());
    estimators_.clear();
    UpdateEstimator(block);
  }
}
