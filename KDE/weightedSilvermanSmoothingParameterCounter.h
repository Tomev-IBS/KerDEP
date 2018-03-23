#ifndef WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H
#define WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H

#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"

#include "smoothingParameterCounter.h"

#include <QVector>

class weightedSilvermanSmoothingParameterCounter
    : public smoothingParameterCounter
{
  public:

    weightedSilvermanSmoothingParameterCounter(QVector<qreal> *samples,
                                               QVector<double> *weights);

    weightedSilvermanSmoothingParameterCounter(std::vector<std::shared_ptr<cluster> > *clusters,
                                               int dimension);

    double countSmoothingParameterValue();

  protected:

    QVector<qreal> *samples;
    QVector<double> *weights;

    double countWeightedStandardDeviation();

};

#endif // WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H
