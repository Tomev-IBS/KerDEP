#include "enhancedClusterKernelAlgorithm.h"

#include <cmath>

EnhancedClusterKernelAlgorithm::EnhancedClusterKernelAlgorithm(const int &m,
                                                               ClusterKernel *(*clusterKernelFactoryMethod)(ClusterKernelStreamElement *stream_element))
 : WeightedUnivariateListBasedClusterKernelAlgorithm(m, clusterKernelFactoryMethod)
{}

std::vector<Point> EnhancedClusterKernelAlgorithm::GetErrorDomain(const int &dimension) {
  double minimal_value = FindMinimalValueOnDimension(dimension);
  double maximal_value = FindMaximalValueOnDimension(dimension);
  // Similarly to 5 sigma.
  minimal_value -= 5 * bandwidth_[dimension];
  maximal_value += 5 * bandwidth_[dimension];
  // We want hard 1000 points, so:
  double step_size = (maximal_value - minimal_value) / 1000;
  auto domain = std::vector<Point>();

  for(auto current_value = minimal_value; current_value < maximal_value; current_value += step_size){
    domain.push_back({current_value});
  }

  return domain;
}

Point EnhancedClusterKernelAlgorithm::GetKDEValuesOnDomain(std::vector<Point> domain) {
  Point kde_values_on_domain = {};

  for(const Point& pt : domain){
    // Note that KDE returned values should be one dimensional
    double value_on_point = GetValue(pt)[0];
    kde_values_on_domain.push_back(value_on_point);
  }

  return kde_values_on_domain;
}

double EnhancedClusterKernelAlgorithm::GetStandardDeviation() {
  return sqrt(variation_estimator_[0]);
}

double EnhancedClusterKernelAlgorithm::GetBandwidth() {
  return bandwidth_[0];
}
