#ifndef DESDARESERVOIR_H
#define DESDARESERVOIR_H

#include "KDE/kerneldensityestimator.h"

using std::vector;

class DESDAReservoir
{
  public:
    DESDAReservoir();

  private:
    unsigned int _minimalReservoirSize = 20; // m_Min
    unsigned int _maximalReservoirSize = 1000; // m_Max
    unsigned int _currentStep = 0;

    // k new clusters, that appeared before algorithm step
    vector<clusterPtr> _newClusters;
    // main reservoir
    vector<clusterPtr> _clusters;

};

#endif // DESDARESERVOIR_H
