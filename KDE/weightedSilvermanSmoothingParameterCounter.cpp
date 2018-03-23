#include "weightedSilvermanSmoothingParameterCounter.h"
#include "math.h"

#include "iostream"

weightedSilvermanSmoothingParameterCounter::weightedSilvermanSmoothingParameterCounter
(QVector<qreal> *samples, QVector<double> *weights) : samples(samples), weights(weights)
{}

weightedSilvermanSmoothingParameterCounter::weightedSilvermanSmoothingParameterCounter(std::vector<std::shared_ptr<cluster>> *clusters, int dimension)
{
  samples = new QVector<qreal>();
  weights = new QVector<double>();

  for(std::shared_ptr<cluster> c : *clusters)
  {
    int counter = 0;

    for(auto attrVal : c.get()->getObject().get()->attributesValues)
    {
      if(counter == dimension) samples->append(stod(attrVal.second));
    }

    weights->append(c.get()->getWeight());
  }
}

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
  double weightsSum = 0;

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
