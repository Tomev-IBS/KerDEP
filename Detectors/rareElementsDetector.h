#ifndef RAREELEMENTSDETECTOR_H
#define RAREELEMENTSDETECTOR_H

#include "./KDE/kerneldensityestimator.h"

typedef QVector<qreal> point;

class rareElementsDetector
{
  public:
    rareElementsDetector(kernelDensityEstimator *KDE, qreal r);

    int findAtypicalElementsInDomain(QVector<point *> *domain,
                                     QVector<int> *target);

  private:
    kernelDensityEstimator *KDE;
    qreal r;

    int getEstimatedValuesFromDomain(QVector<point *> *domain,
                                     QVector<qreal> *estimatedValues);

    int countI(int domainSize);
    qreal countQ(int i, QVector<qreal> *estimatedValues);



};

#endif // RAREELEMENTSDETECTOR_H
