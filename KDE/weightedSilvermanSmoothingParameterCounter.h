#ifndef WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H
#define WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H

#include "smoothingParameterCounter.h"

#include <QVector>

class weightedSilvermanSmoothingParameterCounter
    : public smoothingParameterCounter
{
  public:

    weightedSilvermanSmoothingParameterCounter(QVector<qreal> *samples,
                                               QVector<int> *weights);

    double countSmoothingParameterValue();

  protected:

    QVector<qreal> *samples;
    QVector<int> *weights;

    double countWeightedStandardDeviation();

};

#endif // WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H
