//
// Created by tomev on 1/24/2023.
//

#ifndef KERDEP_WEIGHTEDCVBANDWIDTHSELECTOR_H
#define KERDEP_WEIGHTEDCVBANDWIDTHSELECTOR_H

#include <armadillo>
#include <vector>

using namespace arma;
using std::vector;

class WeightedCVBandwidthSelector {
  // Weighted version of the simplified version from Kulczycki's book.
  public:
    WeightedCVBandwidthSelector();
    double compute_bandwidth (const vector<vec> &data, const vector<double> &weights);
  protected:
    vector<double> _tested_bandwidths = {};

    vector<vec> _data;
    vector<double> _weights;

    double _k0 = 0;

    double g(const double &h) const;
    static double G(const vec &x);
    static double normal_kernel_value(const vec &x);
    static double normal_kernel_squared_convolution_value(const vec &x);
};

#endif //KERDEP_WEIGHTEDCVBANDWIDTHSELECTOR_H
