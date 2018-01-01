#ifndef VELOCITYDENSITYESTIMATOR_H
#define VELOCITYDENSITYESTIMATOR_H

#include "../KDE/kerneldensityestimator.h"
#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"

#include <memory>
#include <vector>

class velocityDensityEstimator
{
  public:

    velocityDensityEstimator(kernelDensityEstimator *kde, long time);

    double countVelocitiDensityFromClusters(std::vector<std::shared_ptr<cluster>> clusters);

    long setTime(long time);

  protected:

    long temporalWindow;
    long time;

    std::shared_ptr<kernelDensityEstimator> KDE;

    long countTemporalWindowFromClusters(
      std::vector<std::shared_ptr<cluster>> clusters
    );
    double countForwardTimeSliceDensityFromClusters(
      std::vector<std::shared_ptr<cluster>> clusters
    );
      double countSpatiotemporalKernelValueFromCluster(
        std::shared_ptr<cluster> c, long moment);
    double countReverseTimeSliceDensityFromClusters(
      std::vector<std::shared_ptr<cluster>> clusters
    );
};

#endif // VELOCITYDENSITYESTIMATOR_H
