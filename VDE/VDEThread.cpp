#include "VDEThread.h"

VDEThread::VDEThread
  (kernelDensityEstimator *kde,
   std::map<long, std::map<point, double> > *temporalVelocityDensityProfile,
   std::vector<std::shared_ptr<cluster> > clusters)
{
  VDE.reset(new velocityDensityEstimator(kde, temporalVelocityDensityProfile));
  this->clusters = clusters;
}

void VDEThread::run()
{
  VDE->countTemporalVelocityDensityProfileFromClusters(clusters);
}

long VDEThread::setTime(long time)
{
  return VDE->setTime(time);
}

QVector<std::shared_ptr<point> > *VDEThread::getDomainPtr()
{
  return VDE->getDomainPtr();
}
