#include "velocityDensityEstimator.h"

velocityDensityEstimator::velocityDensityEstimator(
    std::shared_ptr<kernelDensityEstimator> kde,
    long time)
{
  KDE.reset(kde);
  this->time = time;
}

double velocityDensityEstimator::countVelocitiDensityFromClusters(
    std::vector<std::shared_ptr<cluster> > clusters)
{
  countTemporalWindowFromClusters(clusters);

  double result = 0.0;

  return result;
}

long velocityDensityEstimator::setTime(long time)
{
  this->time = time;
}

long velocityDensityEstimator::countTemporalWindowFromClusters(std::vector<std::shared_ptr<cluster> > clusters)
{
  long largestTimestamp, smallestTimestamp = largestTimestamp = clusters[0]->getTimestamp();
  long clustersTimestamp;

  for(std::shared_ptr<cluster> c : clusters)
  {
    clustersTimestamp = c->getTimestamp();

    if(clustersTimestamp > largestTimestamp)
      largestTimestamp = clustersTimestamp;

    if(clustersTimestamp < smallestTimestamp)
      smallestTimestamp = clustersTimestamp;
  }

  temporalWindow = largestTimestamp - smallestTimestamp;

  return this->temporalWindow;
}

double velocityDensityEstimator::countForwardTimeSliceDensityFromClusters(
    std::vector<std::shared_ptr<cluster> > clusters)
{
  double result = 0.0;

  for(std::shared_ptr<cluster> c: clusters)
  {

  }

  return result;
}

double velocityDensityEstimator::countSpatiotemporalKernelValueFromCluster(
    std::shared_ptr<cluster> c, long moment)
{
  double result = 0.0;



  result *= (1.0 - (double) moment / (double) temporalWindow);

  return result;
}
