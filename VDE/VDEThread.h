#ifndef VDETHREAD_H
#define VDETHREAD_H

#include <QThread>

#include "VDE/velocityDensityEstimator.h"

class VDEThread : public QThread
{
  public:
    VDEThread(kernelDensityEstimator *kde,
              std::map<long, std::map<point, double>>* temporalVelocityDensityProfile,
              std::vector<std::shared_ptr<cluster>> clusters);

    void run();

    long setTime(long time);

    QVector<std::shared_ptr<point>>* getDomainPtr();

  protected:

    std::shared_ptr<velocityDensityEstimator> VDE;
    std::vector<std::shared_ptr<cluster>> clusters;

};

#endif // VDETHREAD_H
