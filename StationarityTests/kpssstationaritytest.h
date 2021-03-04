#ifndef KPSSSTATIONARITYTEST_H
#define KPSSSTATIONARITYTEST_H

#include<vector>

#include "i_stationaritytest.h"

class KPSSStationarityTest : public i_stationarityTest
{
  public:
    KPSSStationarityTest(int maxM);
    double getTestsValue();
    void addNewSample(double sample);
    void setSampleSize(int newSize);
  private:
    int _l = 0;
    int _maxM = 0;
    std::vector<double> _samples = {};
    std::vector<double> _regressionRests = {};

    void calculateRegressionRests();
    double getSumOfSquaredRegressionRests();
    double getLongRunVarianceEstimator();
    double getBarlettWindow(int s, int l);
};

#endif // KPSSSTATIONARITYTEST_H
