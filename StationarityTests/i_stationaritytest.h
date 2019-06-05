#ifndef I_STATIONARITYTEST_H
#define I_STATIONARITYTEST_H

#include <memory>

class i_stationarityTest
{
  public:
    virtual double getTestsValue() = 0;
    virtual double addNewSample(double sample) = 0;
    virtual void setSampleSize(int newSize) = 0;
};

typedef std::shared_ptr<i_stationarityTest> stationarityTestPtr;

#endif // I_STATIONARITYTEST_H
