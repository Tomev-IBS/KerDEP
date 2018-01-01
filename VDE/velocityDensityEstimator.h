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

    double countVelocitiDensityFromClusters(std::vector<std::shared_ptr<cluster>> clusters);

  protected:

    long temporalWindow;

    std::shared_ptr<kernelDensityEstimator> KDE;

    long countTemporalWindowFromClusters(
      std::vector<std::shared_ptr<cluster>> clusters
    );
    double countForwardTimeSliceDensityFromClusters(
      std::vector<std::shared_ptr<cluster>> clusters
    );
      double countSpatiotemporalKernelValueFromCluster(
        std::shared_ptr<cluster> c
      );
    double countReverseTimeSliceDensityFromClusters(
      std::vector<std::shared_ptr<cluster>> clusters
    );
};

#endif // VELOCITYDENSITYESTIMATOR_H
