#ifndef PKSAMPLINGALGORITHM_H
#define PKSAMPLINGALGORITHM_H

#include <random>
#include "reservoirSamplingAlgorithm.h"

class PKSamplingAlgorithm : public reservoirSamplingAlgorithm
{
  public:
    PKSamplingAlgorithm(dataReader *reader, dataParser* parser, int reservoirMaxSize, int stepsNumber);

    void fillReservoir(std::vector<std::shared_ptr<sample>> *reservoir);
    void performSingleStep(std::vector<std::shared_ptr<sample>> *reservoir, int stepNumber);
    unsigned int getReservoidMaxSize();
    void changeReservoirMaxSize(unsigned int newMaxSize);

  private:

    unsigned int reservoirMaxSize = 1000;
    void addDatumToReservoir(std::vector<std::shared_ptr<sample> > *reservoir);
};

#endif // PKSAMPLINGALGORITHM_H
