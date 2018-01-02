#ifndef VELOCITYDENSITYESTIMATOR_H
#define VELOCITYDENSITYESTIMATOR_H

#include "../KDE/kerneldensityestimator.h"
#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"

#include <memory>
#include <vector>
#include <map>

class velocityDensityEstimator
{
  public:

    velocityDensityEstimator(kernelDensityEstimator *kde, long time);

    double countTemporalVelocityDensityProfileFromClusters(
        std::vector<std::shared_ptr<cluster>> clusters);

    long setTime(long time);

    QVector<std::shared_ptr<point>>* getDomainPtr();

    std::map<long, std::map<point, double> > *getTemporalVelocityDensityProfilePtr();

  protected:

    long temporalWindow, time;

    QVector<std::shared_ptr<point>> domain;

    std::map<long, std::map<point, double>> temporalVelocityDensityProfile;

    std::shared_ptr<kernelDensityEstimator> KDE;

    long countTemporalWindowFromClusters
      (std::vector<std::shared_ptr<cluster>> clusters);

    double countForwardTimeSliceDensityInPoint(std::vector<std::shared_ptr<cluster>> clusters, std::shared_ptr<point> pt);

      QVector<qreal> countSpatialLocationForTimeSliceComputation
        (std::shared_ptr<point> pt, std::shared_ptr<cluster> c);

      double countSpatiotemporalKernelValue(std::shared_ptr<point> pt, long moment);

    double countReverseTimeSliceDensityInPoint
      (std::vector<std::shared_ptr<cluster>> clusters, std::shared_ptr<point> pt);
};

#endif // VELOCITYDENSITYESTIMATOR_H
