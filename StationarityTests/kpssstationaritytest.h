#ifndef KPSSSTATIONARITYTEST_H
#define KPSSSTATIONARITYTEST_H

#include<vector>

#include "i_stationaritytest.h"

class KPSSStationarityTest : public i_stationarityTest
{
  public:
    KPSSStationarityTest(int maxM, double &avg, int l = 0);
    double getTestsValue();
    void addNewSample(double sample);
    void setSampleSize(int newSize);
  private:
    int _l = 0;
    int _maxM = 0;
    double & _avg;
    std::vector<double> _regressionRests = {};

    double getSumOfRegressionRests();
    double getLongRunVarianceEstimator();
    double getBarlettWindow(int s, int l);
};

#endif // KPSSSTATIONARITYTEST_H
