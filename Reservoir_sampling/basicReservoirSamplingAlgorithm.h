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

  basicReservoirSamplingAlgorithm(dataReader *reader, dataParser* parser, int reservoirSize, int stepsNumber);
    void fillReservoir(void *reservoir);

  private:

    int reservoirSize = 1000;

    void initializeReservoir(void *reservoir);

};


#endif //RESERVOIRALGORITHM_BASICRESERVOIRSAMPLINGALGORITHM_H
