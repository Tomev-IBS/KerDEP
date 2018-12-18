#define _USE_MATH_DEFINES

#include "weightedSilvermanSmoothingParameterCounter.h"
#include <cmath>
#include <qdebug.h>

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

void weightedSilvermanSmoothingParameterCounter::setClusters(std::vector<std::shared_ptr<cluster> > *clusters, int dimension)
{
  samples->clear();
  weights->clear();

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

weightedSilvermanSmoothingParameterCounter::~weightedSilvermanSmoothingParameterCounter()
{
  delete samples;
  delete weights;
}

double weightedSilvermanSmoothingParameterCounter::countSmoothingParameterValue()
{
  //double result = 0.9 * pow(samples->size(), -0.2);
  double result = pow(4.0 / (3 * samples->size()), 0.2);

  result *= countWeightedStandardDeviation();

  return result;
}

double weightedSilvermanSmoothingParameterCounter::countWeightedStandardDeviation()
{
  if(weights == nullptr)
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

  // Weight's has to be normalized to n = number_of_samples, so...
  double weightsNormalizationFactor = samples->size();
  std::vector<double> normalizedWeights;
  double weightsSum = 0;

  for(auto weight : *weights)
    weightsSum += weight;

  weightsNormalizationFactor /= weightsSum;
  qDebug() << "factor = " << weightsNormalizationFactor;


  // Normalizing all weights.
  for(long i = 0; i < weights->size(); ++i)
    normalizedWeights.push_back((*weights)[i] * weightsNormalizationFactor);

  // Count weighted stDev
  double stDev = 0;

  for(long i = 0; i < samples->size(); ++i)
    stDev += pow(normalizedWeights[i] * (*samples)[i], 2);

  stDev /= samples->size() - 1;

  double substract = 0;

  for(long i = 0; i < samples->size(); ++i)
    substract += normalizedWeights[i] * (*samples)[i];

  substract = pow(substract, 2);

  substract /= samples->size() * (samples->size() - 1);

  stDev -= substract;

  stDev = pow(stDev, 0.5);

  return stDev;
}
