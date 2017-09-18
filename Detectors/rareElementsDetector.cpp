#include "rareElementsDetector.h"
#include "qdebug.h"

#include <QtMath>

rareElementsDetector::rareElementsDetector(kernelDensityEstimator *KDE,
                                           qreal r)
{
  this->KDE = KDE;
  this->r = r;
}

int rareElementsDetector::findAtypicalElementsInDomain(QVector<point *>* domain,
                                                       QVector<int>* target)
{
  // It's better to store estimated values in order not to count them
  // again when comparing to q.

  QVector<qreal> estimatedValues, sortedValues;

  getEstimatedValuesFromDomain(domain, &estimatedValues);

  sortedValues = estimatedValues;
  qSort(sortedValues);

  int i = countI(domain->size());
  qreal q = countQ(i, &estimatedValues);

  target->clear();

  for(int index = 0; index < estimatedValues.size(); ++index)
  {
    if(estimatedValues.at(index) <= q)
      target->append(index);
  }

  return target->size();
}

int rareElementsDetector::getEstimatedValuesFromDomain(QVector<point *> *domain,
                                                       QVector<qreal> *estimatedValues)
{
  estimatedValues->clear();

  foreach(auto arg, *domain)
    estimatedValues->append(KDE->getValue(arg));

  return estimatedValues->size();
}

int rareElementsDetector::countI(int domainSize)
{
  return qFloor(domainSize * r + 0.5);
}

qreal rareElementsDetector::countQ(int i, QVector<qreal> *estimatedValues)
{
  qreal mr = estimatedValues->size() * r;

  if(mr <= 0.5) return estimatedValues->at(i);

  qreal q = (0.5 + i - mr) * estimatedValues->at(i);
  q += (0.5 - i + mr) * estimatedValues->at(i+1);

  return q;
}
