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
  return VectorScalarMultiplication(elements_sum_, 1.0 / cardinality_);
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
  return new VarianceBasedClusterKernel(
    cardinality_ + other_variance_based_cluster_kernel->GetCardinality(),
    SumVectors(elements_sum_, other_variance_based_cluster_kernel->GetElementsSum()),
    SumVectors(elements_squared_sum_, other_variance_based_cluster_kernel->GetElementsSquaredSum())
  );
}

Point VarianceBasedClusterKernel::GetValue(const Point &pt) {
  if(pt.size() != bandwidth_.size()){
    return {};
  }

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
}

double VarianceBasedClusterKernel::GetWeight() {
  return GetCardinality(); // This is an unweighted version.
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
