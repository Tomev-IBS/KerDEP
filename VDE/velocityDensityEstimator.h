#ifndef VELOCITYDENSITYESTIMATOR_H
#define VELOCITYDENSITYESTIMATOR_H

#include "../KDE/kerneldensityestimator.h"
#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"

#include <memory>
#include <vector>


class velocityDensityEstimator
{
  public:
    velocityDensityEstimator(std::shared_ptr<kernelDensityEstimator> kde);

    double countVelocitiDensityFromClusters(std::shared_ptr<cluster> clusters);

  protected:

    std::shared_ptr<kernelDensityEstimator> KDE;

    double countForwardTimeSliceDensityFromClusters(std::shared_ptr<cluster> clusters);
    double countReverseTimeSliceDensityFromClusters(std::shared_ptr<cluster> clusters);
};

#endif // VELOCITYDENSITYESTIMATOR_H
