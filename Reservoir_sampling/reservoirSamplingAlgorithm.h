//
// Created by Tomev on 29.05.2017.
//

#ifndef RESERVOIRALGORITHM_RESERVOIRSAMPLINGALGORITHM_H
#define RESERVOIRALGORITHM_RESERVOIRSAMPLINGALGORITHM_H

#include "dataReader.h"
#include "dataParser.h"

#include <random>

class reservoirSamplingAlgorithm
{

  public:

    reservoirSamplingAlgorithm(dataReader *reader, dataParser* parser);
    void fillReservoir(void *reservoir);

  private:
    const int RESERVOIR_SIZE = 1000;
    const int STEPS_NUMBER = 20000;

    dataReader *reader;
    dataParser *parser;


    void initializeReservoir(void *reservoir);

};


#endif //RESERVOIRALGORITHM_RESERVOIRSAMPLINGALGORITHM_H
