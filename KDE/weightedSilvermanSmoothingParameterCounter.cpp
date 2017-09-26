#include "weightedSilvermanSmoothingParameterCounter.h"
#include "math.h"

#include "iostream"

weightedSilvermanSmoothingParameterCounter::weightedSilvermanSmoothingParameterCounter
(QVector<qreal> *samples, QVector<int> *weights) : samples(samples), weights(weights)
{}

double weightedSilvermanSmoothingParameterCounter::countSmoothingParameterValue()
{
  double result = 0.9 * pow(samples->size(), -0.2);

  result *= countWeightedStandardDeviation();

  return result;
}

double weightedSilvermanSmoothingParameterCounter::countWeightedStandardDeviation()
{
  if(weights == NULL)
  {
    std::cerr  << "Weighted Silverman smoothing parameter counter: "
          << "Weights pointer is null. Cannot perform operation." << std::endl;

    return -1.0;
  }

  if(samples->size() != weights->size())
  {
    std::cerr << "Weighted Silverman smoothing parameter counter: "
         << "Weights and samples are not of the same size." << std::endl;

    return -2.0;
  }

  double stDev = 0;
  int weightsSum = 0;

  for(int i = 0; i < samples->size(); ++i)
  {
    stDev += weights->at(i) * pow(samples->at(i), 2);
    weightsSum += weights->at(i);
  }

  stDev /= weightsSum - 1;

  double substract = 0;

  for(int i = 0; i < samples->size(); ++i)
    substract += pow(weights->at(i) * samples->at(i), 2);

  substract /= weightsSum * (weightsSum - 1);

  stDev -= substract;

  stDev = pow(stDev, 0.5);

  return stDev;
}
