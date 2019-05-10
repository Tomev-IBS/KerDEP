#ifndef RESERVOIRALGORITHM_RESERVOIRSAMPLINGALGORITHM_H
#define RESERVOIRALGORITHM_RESERVOIRSAMPLINGALGORITHM_H

#include <vector>
#include <memory>

#include "dataParser.h"
#include "dataReader.h"
#include "sample.h"

class reservoirSamplingAlgorithm
{
  public:
    virtual void fillReservoir(std::vector<std::shared_ptr<sample>> *reservoir) = 0;
    virtual void performSingleStep(std::vector<std::shared_ptr<sample>> *reservoir, int stepNumber) = 0;
    virtual unsigned int getReservoidMaxSize() = 0;
    virtual void changeReservoirMaxSize(unsigned int newMaxSize) = 0;

  protected:
    dataParser *parser;
    dataReader *reader;

    int stepsNumber = 10000;
};

#endif //RESERVOIRALGORITHM_RESERVOIRSAMPLINGALGORITHM_H
