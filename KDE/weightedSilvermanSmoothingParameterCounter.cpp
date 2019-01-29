#define _USE_MATH_DEFINES

#include "weightedSilvermanSmoothingParameterCounter.h"
#include <cmath>
#include <qdebug.h>

#include "iostream"

weightedSilvermanSmoothingParameterCounter::weightedSilvermanSmoothingParameterCounter
(QVector<qreal> *samples, QVector<double> *weights) : _samples(samples), _weights(weights)
{}

weightedSilvermanSmoothingParameterCounter::weightedSilvermanSmoothingParameterCounter(std::vector<std::shared_ptr<cluster>> *clusters, int dimension)
{
  _samples = new QVector<qreal>();
  _weights = new QVector<double>();

  for(std::shared_ptr<cluster> c : *clusters)
  {
    int counter = 0;

    for(auto attrVal : c.get()->getObject().get()->attributesValues)
    {
      if(counter == dimension) _samples->append(stod(attrVal.second));
    }

    _weights->append(c.get()->getWeight());
  }
}

void weightedSilvermanSmoothingParameterCounter::setClusters(std::vector<std::shared_ptr<cluster> > *clusters, int dimension)
{
  _samples->clear();
  _weights->clear();

  for(std::shared_ptr<cluster> c : *clusters)
  {
    int counter = 0;

    for(auto attrVal : c.get()->getObject().get()->attributesValues)
    {
      if(counter == dimension) _samples->append(stod(attrVal.second));
    }

    _weights->append(c.get()->getWeight());
  }
}

weightedSilvermanSmoothingParameterCounter::~weightedSilvermanSmoothingParameterCounter()
{
  delete _samples;
  delete _weights;
}

double weightedSilvermanSmoothingParameterCounter::countSmoothingParameterValue()
{
  //double result = pow(4.0 / (3 * _samples->size()), 0.2);
  //countStandardDeviation();

  double weightsSum = 0.0;
  for (auto w : *_weights) weightsSum += w;

  double result = pow(4.0 / (3 * weightsSum), 0.2);
  recountWeightedStandardDeviation();

  return result * _stDev;
}

void weightedSilvermanSmoothingParameterCounter::recountWeightedStandardDeviation()
{
  if(_samples->size() == 1)
  {
    _stDev = 1.0;
    return;
  }

  if(_weights == nullptr)
  {
    std::cerr  << "Weighted Silverman smoothing parameter counter: "
          << "Weights pointer is null. Cannot perform operation." << std::endl;

    return;
  }

  if(_samples->size() != _weights->size())
  {
    std::cerr << "Weighted Silverman smoothing parameter counter: "
         << "Weights and samples are not of the same size." << std::endl;

    return;
  }

  double weightsSum = 0.0;

  for(auto weight : *_weights)
    weightsSum += weight;

  double weightedMean = 0.0;

  for(int i = 0; i < _samples->size(); ++i)
    weightedMean += (*_samples)[i] * (*_weights)[i];

  weightedMean /= weightsSum;

  double var = 0.0;

  for(int i = 0; i < _samples->size(); ++i)
    var += (*_weights)[i] * pow((*_samples)[i] - weightedMean, 2.0);

  double N = _samples->size(); // So I don't have to cast during division

  var /= weightsSum * (N - 1) / N;

  _stDev = pow(var, 0.5);
}

void weightedSilvermanSmoothingParameterCounter::countStandardDeviation()
{
  double m = _samples->size();
  _stDev = 0.0;

  double substract = 0.0;

  for(auto sample : *_samples)
  {
    _stDev += pow(sample, 2.0);
    substract += sample;
  }

  _stDev /= (m - 1);
  substract = pow(substract, 2) / (m * (m - 1));

  _stDev = pow(_stDev - substract, 0.5);
}

void weightedSilvermanSmoothingParameterCounter::updateSmoothingParameterValue(double weightModifier, double newSample)
{
  if(_stDev - 1e-10 < 0)
  {
    _stDev = 1.0;
    _m = 1.0;
    _sampleSum = newSample;
    _squaredSampleSum = pow(newSample, 2.0);
    _h = pow(4.0 / (3 * _m), 0.2) * _stDev;
    return;
  }

  _m = _m * weightModifier + 1.0;
  _squaredSampleSum = _squaredSampleSum * weightModifier + pow(newSample, 2.0);
  _sampleSum = _sampleSum * weightModifier + newSample;
  _stDev = _squaredSampleSum / (_m - weightModifier);
  _stDev -= pow(_sampleSum, 2.0) / (_m * (_m - weightModifier));
  _stDev = pow(_stDev, 0.5);

  _h = pow(4.0 / (3 * _m), 0.2) * _stDev;
}

double weightedSilvermanSmoothingParameterCounter::getSmoothingParameterValue()
{
  return _h;
}
