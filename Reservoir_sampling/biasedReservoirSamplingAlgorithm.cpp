//
// Created by Tomev on 03.06.2017.
//

#include "biasedReservoirSamplingAlgorithm.h"

biasedReservoirSamplingAlgorithm::biasedReservoirSamplingAlgorithm(
        dataReader *reader, dataParser *parser, int reservoirSize, int stepsNumber):
        biasRate(1.0/reservoirSize)
{
  this->reader = reader;
  this->parser = parser;
  this->stepsNumber = stepsNumber;
}

void biasedReservoirSamplingAlgorithm::fillReservoir(void *reservoir)
{
    int indexOfSampleToWriteOn, currentReservoirSize = 0;

    // For each incoming data sample
    for(int step = 0; step < stepsNumber; ++step)
    {
      updateFractionOfReservoirFilled(currentReservoirSize);

      // Check if a sample should be removed
      if(fractionOfReservoirFilled >= ((double) rand() / (RAND_MAX)))
      {
          // If so randomly choose a sample
          indexOfSampleToWriteOn = (((double) rand() / (RAND_MAX)) * (currentReservoirSize - 1));
      }
      else
      {
        // If not add new sample to reservoir
        currentReservoirSize = parser->addDatumToContainer(reservoir);
        indexOfSampleToWriteOn  = currentReservoirSize - 1;
      }

      // Add new data to the reservoir
      reader->getNextRawDatum(parser->buffor);
      parser->writeDatumOnPosition(reservoir, indexOfSampleToWriteOn);
    }
}

void biasedReservoirSamplingAlgorithm::updateFractionOfReservoirFilled(int currentReservoirSize)
{
    // Bias rate is equal to 1 over reservoir capacity, thus:
    this->fractionOfReservoirFilled = currentReservoirSize * biasRate;
}


