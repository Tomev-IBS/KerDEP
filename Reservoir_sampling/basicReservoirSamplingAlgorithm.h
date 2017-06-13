//
// Created by Tomev on 29.05.2017.
//

#ifndef RESERVOIRALGORITHM_BASICRESERVOIRSAMPLINGALGORITHM_H
#define RESERVOIRALGORITHM_BASICRESERVOIRSAMPLINGALGORITHM_H

#include <random>

#include "reservoirSamplingAlgorithm.h"

class basicReservoirSamplingAlgorithm : public reservoirSamplingAlgorithm
{

  public:

    basicReservoirSamplingAlgorithm(dataReader *reader, dataParser* parser, int reservoirMaxSize, int stepsNumber);

    void fillReservoir(std::vector<sample*> *reservoir);
    void performSingleStep(std::vector<sample*> *reservoir, int stepNumber);

  private:

    int reservoirMaxSize = 1000;

    void initializeReservoir(std::vector<sample*> *reservoir);

    void addDatumToReservoir(std::vector<sample *> *reservoir);
    bool shouldDatumBeAdded(int stepNumber);
    void addDatumOnRandomPosition(std::vector<sample *> *reservoir);
};


#endif //RESERVOIRALGORITHM_BASICRESERVOIRSAMPLINGALGORITHM_H
