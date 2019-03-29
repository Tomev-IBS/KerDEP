#include "kpssstationaritytest.h"

#include <cmath>

KPSSStationarityTest::KPSSStationarityTest(int maxM, double &avg, int l)
  : _avg(avg), _l(l), _maxM(maxM)
{}

// It's eta with hat from 1992 KPSS work.
double KPSSStationarityTest::getTestsValue()
{
  if(_regressionRests.size() == 0)
    return 0;

  double testValue = 0.0;
  int m = std::min(_maxM, static_cast<int>(_regressionRests.size()));

  testValue = getSumOfRegressionRests();
  testValue /=  getLongRunVarianceEstimator();
  testValue /= m * m;

  return testValue;
}

double KPSSStationarityTest::addNewSample(double sample)
{
  _regressionRests.push_back(sample - _avg);

  if(_regressionRests.size() == _maxM){
    _regressionRests.erase(_regressionRests.begin(), _regressionRests.begin() + 1);
  }
}

// It's sum of S_t^2
double KPSSStationarityTest::getSumOfRegressionRests()
{
  double ithSum = 0;
  double sumOfSquaredSums = 0;

  for(auto val : _regressionRests) {
    ithSum += val;
    sumOfSquaredSums += ithSum * ithSum;
  }

  return sumOfSquaredSums;
}

// It's s^2(l) from 1992 KPSS work.
double KPSSStationarityTest::getLongRunVarianceEstimator()
{
  double estimator = 0;

  int m = std::min(_maxM, static_cast<int>(_regressionRests.size()));

  for(auto val : _regressionRests)
    estimator += val * val;

  double esSum = 0;

  for(auto s = 0; s < _l; ++s){

    esSum = 0;

    for(auto i = s + 1; i < m; ++i)
      esSum += (_regressionRests[i]) * (_regressionRests[i - s]);

    estimator += 2 * esSum * getBarlettWindow(s, _l);
  }

  estimator /= m;

  return estimator;
}

// It's w(s, l) form 1992 KPSS work.
double KPSSStationarityTest::getBarlettWindow(int s, int l)
{
  double window = 1.0;

  window -= static_cast<double>(s) / (l + 1);

  return window;
}
