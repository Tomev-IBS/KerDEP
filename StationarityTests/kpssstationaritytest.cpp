#include "kpssstationaritytest.h"

#include <cmath>
#include <QDebug>

KPSSStationarityTest::KPSSStationarityTest(int maxM)
  : _maxM(maxM)
{}

// It's eta with hat from 1992 KPSS work.
double KPSSStationarityTest::getTestsValue()
{
  int T = _samples.size();

  if(T < 2) return 0;

  double testValue = 0.0;

  _l = round(4 * pow(T / 100, 0.25));

  calculateRegressionRests();
  testValue = getSumOfSquaredRegressionRests();
  testValue /= T * T;

  testValue /= getLongRunVarianceEstimator();

  return testValue;
}

void KPSSStationarityTest::addNewSample(double sample)
{
  _samples.push_back(sample);

  while(_samples.size() > _maxM)
    _samples.erase(_samples.begin(), _samples.begin() + 1);
}

void KPSSStationarityTest::setSampleSize(int newSize)
{
  _maxM = newSize;
}

void KPSSStationarityTest::calculateRegressionRests()
{
  _regressionRests.clear();

  double averageSampleValue = 0;

  for(auto val : _samples)
    averageSampleValue += val;

  averageSampleValue /= _samples.size();

  for(auto val : _samples)
    _regressionRests.push_back(val - averageSampleValue);
}

// It's sum of S_t^2
double KPSSStationarityTest::getSumOfSquaredRegressionRests()
{
  double ithSum = 0;
  double sumOfSquaredRests = 0;

  for(auto val : _regressionRests) {
    ithSum += val;
    sumOfSquaredRests += ithSum * ithSum;
  }

  return sumOfSquaredRests;
}

// It's s^2(l) from 1992 KPSS work.
double KPSSStationarityTest::getLongRunVarianceEstimator()
{
  double estimator = 0;

  int T = _samples.size();

  for(auto val : _regressionRests)
    estimator += val * val;

  double esSum = 0;

  for(auto s = 0; s < _l; ++s){

    esSum = 0;

    for(auto t = s + 1; t < T; ++t)
      esSum += (_regressionRests[t]) * (_regressionRests[t - s]);

    estimator += 2 * esSum * getBarlettWindow(s, _l);
  }

  estimator /= T;

  return estimator;
}

// It's w(s, l) form 1992 KPSS work.
double KPSSStationarityTest::getBarlettWindow(int s, int l)
{
  double window = 1.0;

  window -= static_cast<double>(s) / (l + 1);

  return window;
}
