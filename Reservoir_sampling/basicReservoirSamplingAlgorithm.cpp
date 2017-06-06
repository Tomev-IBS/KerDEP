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

  double addChance;

  for(int step = reservoirSize+1; step < stepsNumber; ++step)
  {
    addChance = (double) reservoirSize/step;

    reader->getNextRawDatum(parser->buffor);

    if(addChance > ((double) rand() / (RAND_MAX)))
    {
      // Adding to reservoir
      int idxToDelete = (((double) rand() / (RAND_MAX)) * (reservoirSize - 1));

      parser->writeDatumOnPosition(reservoir, idxToDelete);
    }
  }
}

void basicReservoirSamplingAlgorithm::initializeReservoir(void *reservoir)
{
    for(int step = 0; step < reservoirSize; ++step)
    {
        reader->getNextRawDatum(parser->buffor);

        parser->addDatumToContainer(reservoir);

        parser->writeDatumOnPosition(reservoir, step);
    }
}
