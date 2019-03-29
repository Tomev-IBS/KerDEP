#ifndef I_STATIONARITYTEST_H
#define I_STATIONARITYTEST_H

class i_stationarityTest
{
  public:
    virtual double getTestsValue() = 0;
    virtual void addNewItem(double item) = 0;
};

#endif // I_STATIONARITYTEST_H
