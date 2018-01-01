#include "velocityDensityEstimator.h"

velocityDensityEstimator::velocityDensityEstimator(std::shared_ptr<kernelDensityEstimator> kde)
{
  KDE = kde;
}

double velocityDensityEstimator::countVelocitiDensityFromClusters(std::shared_ptr<cluster> clusters)
{

}
