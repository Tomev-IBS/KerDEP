#ifndef WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H
#define WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H

#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"

#include "smoothingParameterCounter.h"

#include <QVector>

class weightedSilvermanSmoothingParameterCounter
    : public smoothingParameterCounter
{
  public:

    weightedSilvermanSmoothingParameterCounter(QVector<qreal> *_samples,
                                               QVector<double> *_weights);

    weightedSilvermanSmoothingParameterCounter(std::vector<std::shared_ptr<cluster> > *clusters,
                                               int dimension);

    void setClusters(std::vector<std::shared_ptr<cluster> > *clusters, int dimension);


    // NEW
    double getSmoothingParameterValue();
    void updateSmoothingParameterValue(double weightModifier, double newSample);


    ~weightedSilvermanSmoothingParameterCounter();

    double _stDev = 0.0;

    double countSmoothingParameterValue() override;

  protected:

    QVector<qreal> *_samples;
    QVector<double> *_weights;

    // NEW
    double _m = 0.0;
    double _sampleSum = 0.0;
    double _squaredSampleSum = 0.0;

    double _h = 0.0;

    void recountWeightedStandardDeviation();
    void countStandardDeviation();
    double updateWeightedStandardDeviation();

};

#endif // WEIGHTEDSILVERMANSMOOTHINGPARAMETERCOUNTER_H
