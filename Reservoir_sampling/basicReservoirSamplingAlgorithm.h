#ifndef RESERVOIRALGORITHM_BASICRESERVOIRSAMPLINGALGORITHM_H
#define RESERVOIRALGORITHM_BASICRESERVOIRSAMPLINGALGORITHM_H

#include <random>

#include "reservoirSamplingAlgorithm.h"

class basicReservoirSamplingAlgorithm : public reservoirSamplingAlgorithm
{
  public:

    basicReservoirSamplingAlgorithm(dataReader *reader, dataParser* parser, int reservoirMaxSize, int stepsNumber);

    void fillReservoir(std::vector<std::shared_ptr<sample>> *reservoir);
    void performSingleStep(std::vector<std::shared_ptr<sample>> *reservoir, int stepNumber);
    unsigned int getReservoidMaxSize();
    void changeReservoirMaxSize(unsigned int newMaxSize);
    std::vector<std::string> *getAttributesList();

  private:

    unsigned int reservoirMaxSize = 1000;

    void addDatumToReservoir(std::vector<std::shared_ptr<sample> > *reservoir);
    bool shouldDatumBeAdded(int stepNumber);
    void addDatumOnRandomPosition(std::vector<std::shared_ptr<sample>> *reservoir);
};


#endif //RESERVOIRALGORITHM_BASICRESERVOIRSAMPLINGALGORITHM_H
