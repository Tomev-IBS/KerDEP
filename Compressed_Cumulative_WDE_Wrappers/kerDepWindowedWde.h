//
// Created by Tomev on 25/01/2021.
//

#ifndef DEDSTA_DEDSTAWINDOWEDWDE_H
#define DEDSTA_DEDSTAWINDOWEDWDE_H

#include "DEDSTACcWde.h"

class DEDSTAWindowedWDE : public DEDSTA_CC_WDE {
  public:
    DEDSTAWindowedWDE(const unsigned int &maximal_number_of_empirical_coefficients,
                      const double &weights_modifier_,
                      WaveletDensityEstimator* (*wde_factory_method)(const vector<double> &values_block),
                      const unsigned int &block_size);

    virtual void PerformStep(point *pt) override;

};

typedef DEDSTAWindowedWDE Windowed_WDE;

#endif //DEDSTA_DEDSTAWINDOWEDWDE_H
