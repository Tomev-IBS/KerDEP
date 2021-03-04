//
// Created by Tomev on 25/01/2021.
//

#ifndef KERDEP_KERDEPWINDOWEDWDE_H
#define KERDEP_KERDEPWINDOWEDWDE_H

#include "kerDepCcWde.h"

class KerDEPWindowedWDE : public KerDEP_CC_WDE {
  public:
    KerDEPWindowedWDE(const unsigned int &maximal_number_of_empirical_coefficients,
                      const double &weights_modifier_,
                      WaveletDensityEstimator* (*wde_factory_method)(const vector<double> &values_block),
                      const unsigned int &block_size);

    virtual void PerformStep(point *pt) override;

};

typedef KerDEPWindowedWDE Windowed_WDE;

#endif //KERDEP_KERDEPWINDOWEDWDE_H
