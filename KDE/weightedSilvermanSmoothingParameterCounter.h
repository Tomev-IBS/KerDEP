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

    void setClusters(std::vector<std::shared_ptr<cluster> > *clusters, int dimension);

    ~weightedSilvermanSmoothingParameterCounter();

    double countSmoothingParameterValue() override;

  protected:

    QVector<qreal> *samples;
    QVector<double> *weights;

    double countWeightedStandardDeviation();

};

#endif // WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H
