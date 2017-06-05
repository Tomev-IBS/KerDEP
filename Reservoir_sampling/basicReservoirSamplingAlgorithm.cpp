//
// Created by Tomev on 29.05.2017.
//

#include "basicReservoirSamplingAlgorithm.h"

basicReservoirSamplingAlgorithm::basicReservoirSamplingAlgorithm(
        dataReader *reader, dataParser *parser, int reservoirSize, int stepsNumber) :
        reservoirSize(reservoirSize)
{
  this->reader = reader;
  this->parser = parser;
  this->stepsNumber = stepsNumber;
}

void basicReservoirSamplingAlgorithm::fillReservoir(void *reservoir)
{
  initializeReservoir(reservoir);

  std::string rawData;

  double addChance;

  for(int step = reservoirSize+1; step < reservoirSize; ++step)
  {
    addChance = (double) reservoirSize/step;
    reader->getNextRawDatum(&rawData);

    if(addChance > ((double) rand() / (RAND_MAX)))
    {
      // Adding to reservoir
      int idxToDelete = (((double) rand() / (RAND_MAX)) * reservoirSize);

      parser->writeDatumOnPosition(&rawData, reservoir, idxToDelete);
    }
  }
}

void basicReservoirSamplingAlgorithm::initializeReservoir(void *reservoir)
{
    std::string rawData;

    for(int step = 0; step < reservoirSize; ++step)
    {
        reader->getNextRawDatum(&rawData);

        parser->addDatumToContainer(reservoir);

        parser->writeDatumOnPosition(&rawData, reservoir, step);
    }
}
