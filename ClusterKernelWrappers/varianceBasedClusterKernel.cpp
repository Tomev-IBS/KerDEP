#include <cmath>

#include "varianceBasedClusterKernel.h"


namespace VarianceBasedClusterKernelUtilities{
  Point ElementwiseVectorsMultiplication(const Point &first_vector, const Point &second_vector){
    if(first_vector.size() != second_vector.size()){
      return {};
    }
    Point multiplied_vector = {};
    for(auto i = 0; i < first_vector.size(); ++i){
      multiplied_vector.push_back(first_vector[i] * second_vector[i]);
    }
    return multiplied_vector;
  }

  Point SumVectors(const Point &first_vector, const Point &second_vector){
    if(first_vector.size() != second_vector.size()){
      return {};
    }
    Point sum = {};

    for(auto i = 0; i < first_vector.size(); ++i){
      sum.push_back(first_vector[i] + second_vector[i]);
    }

    return sum;
  }

  Point VectorScalarMultiplication(const Point &vector, const double &scalar){
    Point multiplied_vector = {};

    for(auto value : vector){
      multiplied_vector.push_back(value * scalar);
    }

    return multiplied_vector;
  }
};

using namespace VarianceBasedClusterKernelUtilities;

VarianceBasedClusterKernel::VarianceBasedClusterKernel(ClusterKernelStreamElement *stream_element, RealValuedFunction *kernel) {
  if(!kernel){
    kernel = new EpanecznikowKernelRealValuedFunction();
  }
  cardinality_ = 1;
  elements_sum_ = stream_element->GetMean();
  elements_squared_sum_ = ElementwiseVectorsMultiplication(elements_sum_, elements_sum_);
  kernel_ = RealValuedFunctionPtr(kernel);
}

Point VarianceBasedClusterKernel::GetMean() {
  return VectorScalarMultiplication(elements_sum_, 1.0 / double(cardinality_));
}

void VarianceBasedClusterKernel::Update(ClusterKernelStreamElement *stream_element) {
  ++cardinality_;
  Point elements_mean = stream_element->GetMean();
  elements_sum_ = SumVectors(elements_sum_, elements_mean);
  elements_squared_sum_ = SumVectors(
    elements_squared_sum_, ElementwiseVectorsMultiplication(elements_mean, elements_mean)
  );
}

ClusterKernel *VarianceBasedClusterKernel::Merge(ClusterKernel *other_cluster_kernel) {
  auto other_variance_based_cluster_kernel = dynamic_cast<VarianceBasedClusterKernel*>(other_cluster_kernel);
  auto merged_cluster = new VarianceBasedClusterKernel(
      cardinality_ + other_variance_based_cluster_kernel->GetCardinality(),
      SumVectors(elements_sum_, other_variance_based_cluster_kernel->GetElementsSum()),
      SumVectors(elements_squared_sum_, other_variance_based_cluster_kernel->GetElementsSquaredSum())
  );
  auto merged_cluster_weight = weight_ * cardinality_;
  merged_cluster_weight += other_cluster_kernel->GetWeight() * other_cluster_kernel->GetCardinality();
  merged_cluster_weight /= merged_cluster->GetCardinality();
  merged_cluster->RescaleWeight(merged_cluster_weight); // Weight is initialized with 1.
  return merged_cluster;
}

Point VarianceBasedClusterKernel::GetValue(const Point &pt) {
  if(pt.size() != bandwidth_.size()){
    return {};
  }

  // There are many resampling methods. Version below is called mean-var-resampling with uniformly distributed values.
  // It's for univariate data only!
  auto mean = GetMean()[0];
  // I need standard deviation.
  auto standard_deviation = 0;
  if(cardinality_ > 1) {
    standard_deviation = sqrt((elements_squared_sum_[0] - pow(elements_sum_[0], 2) / cardinality_) / (cardinality_ - 1));
  }
  // Then I can calculate the domain (or perform the resampling).
  double min_domain_value = mean - 2 * standard_deviation;
  double max_domain_value = mean + 2 * standard_deviation;
  std::vector<Point> domain = {};
  double step_size = (max_domain_value - min_domain_value) / cardinality_;
  double current_point_value = min_domain_value;
  while(domain.size() < cardinality_){
    domain.push_back({(pt[0] - current_point_value) / bandwidth_[0]});
    current_point_value += step_size;
  }
  // Now calculate KDE value from this cluster kernel.
  double kde_value = 0;
  for(const auto& point : domain){
    kde_value += kernel_->GetValue(point)[0];
  }
  kde_value /= bandwidth_[0];
  return {kde_value};

  /*
  // There are many versions of resampling. Here, I'll use the simplest one-value-resampling approach.
  auto point_for_kenrel = SumVectors(pt, VectorScalarMultiplication(GetMean(), -1));
  for(auto i = 0; i < point_for_kenrel.size(); ++i){
    point_for_kenrel[i] = point_for_kenrel[i] / bandwidth_[i];
  }

  // Note, that this should only have one dimension.
  double kernel_returned_value = kernel_->GetValue(point_for_kenrel)[0];
  for(auto value : bandwidth_){
    kernel_returned_value /= value;
  }

  // Multiplying times cardinality instead of sum calculation.
  return {cardinality_ * kernel_returned_value};
 */
}

double VarianceBasedClusterKernel::GetWeight() {
  return weight_;
}

unsigned int VarianceBasedClusterKernel::GetCardinality() {
  return cardinality_;
}

VarianceBasedClusterKernel::VarianceBasedClusterKernel(unsigned int cardinality, Point elements_sum,
                                                      Point elements_squared_sum, RealValuedFunction *kernel) {
  if(!kernel){
    kernel = new EpanecznikowKernelRealValuedFunction();
  }
  cardinality_ = cardinality;
  elements_sum_ = elements_sum;
  elements_squared_sum_ = elements_squared_sum;
  kernel_ = RealValuedFunctionPtr(kernel);
}

Point VarianceBasedClusterKernel::GetElementsSum() {
  return elements_sum_;
}

Point VarianceBasedClusterKernel::GetElementsSquaredSum() {
  return elements_squared_sum_;
}

void VarianceBasedClusterKernel::SetBandwidth(const Point &bandwidth) {
  bandwidth_ = bandwidth;
}

void VarianceBasedClusterKernel::RescaleWeight(const double &modifier) {
  weight_ *= modifier;
}
