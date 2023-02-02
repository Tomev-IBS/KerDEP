//
// Created by tomev on 1/24/2023.
//

#ifndef KERDEP_WEIGHTEDCVBANDWIDTHSELECTOR_H
#define KERDEP_WEIGHTEDCVBANDWIDTHSELECTOR_H

#include <armadillo>

using namespace arma;

class WeightedCVBandwidthSelector {
  public:
    WeightedCVBandwidthSelector();
    mat compute_bandwidth (const mat &data, const vec &weights) const;
  protected:

};

#endif //KERDEP_WEIGHTEDCVBANDWIDTHSELECTOR_H
