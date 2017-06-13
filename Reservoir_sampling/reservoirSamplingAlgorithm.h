//
// Created by Tomev on 05.06.2017.
//

#ifndef RESERVOIRALGORITHM_RESERVOIRSAMPLINGALGORITHM_H
#define RESERVOIRALGORITHM_RESERVOIRSAMPLINGALGORITHM_H

#include <vector>

#include "dataParser.h"
#include "dataReader.h"
#include "sample.h"

class reservoirSamplingAlgorithm
{
  public:
    virtual void fillReservoir(std::vector<sample*> *reservoir) = 0;
    virtual void performSingleStep(std::vector<sample*> *reservoir, int stepNumber) = 0;

  protected:
    dataParser *parser;
    dataReader *reader;

    int stepsNumber = 10000;
};

#endif //RESERVOIRALGORITHM_RESERVOIRSAMPLINGALGORITHM_H
