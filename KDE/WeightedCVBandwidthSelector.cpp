#include "WeightedCVBandwidthSelector.h"

#include <cmath>
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif


WeightedCVBandwidthSelector::WeightedCVBandwidthSelector() {
  double h = 0.05;
  double h_step = 0.05;

  while(h < 1.3){
    _tested_bandwidths.push_back(h);
    h += h_step;
  }

  _k0 = normal_kernel_value(vec(2, fill::zeros));
}

double WeightedCVBandwidthSelector::compute_bandwidth(const vector<vec> &data, const vector<double> &weights) {
  _data = data;
  _weights = weights;

  double min_g_val = g(_tested_bandwidths[0]);
  double h_min = _tested_bandwidths[0];

  for(auto h : _tested_bandwidths){
      double g_val = g(h);

      if(g_val < min_g_val){
        min_g_val = g_val;
        h_min = h;
      }
  }

  return h_min;
}

double WeightedCVBandwidthSelector::g(const double &h) const{
  int m = _data.size();
  int n = _data[0].size(); // No check here.

  double g_value = 0;

  for(int i = 0; i < m; ++i){
    for(int j = 0; j < m; ++j){
      vec argument = _data[j] - _data[i];
      argument /= h;
      g_value += _weights[i] * _weights[j] * G(argument);
    }
  }

  g_value /= m * m;

  g_value += 2 * _k0 / m;

  return g_value / pow(h, n);
}

double WeightedCVBandwidthSelector::G(const vec &x){
  return normal_kernel_squared_convolution_value(x) - 2 * normal_kernel_value(x);
}

double WeightedCVBandwidthSelector::normal_kernel_value(const vec &x) {
  return pow(2 * M_PI, - x.size() / 2) * exp(-as_scalar(x.t() * x) / 2);
}

double WeightedCVBandwidthSelector::normal_kernel_squared_convolution_value(const vec &x) {
  return pow(4 * M_PI, - x.size() / 2) * exp(- as_scalar(x.t() * x) / 4);
}