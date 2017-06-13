//
// Created by Tomev on 03.06.2017.
//

#include "biasedReservoirSamplingAlgorithm.h"

biasedReservoirSamplingAlgorithm::biasedReservoirSamplingAlgorithm(
        dataReader *reader, dataParser *parser, int reservoirMaxSize, int stepsNumber):
        biasRate(1.0/reservoirMaxSize)
{
  this->reader = reader;
  this->parser = parser;
  this->stepsNumber = stepsNumber;
}

void biasedReservoirSamplingAlgorithm::fillReservoir(std::vector<sample*> *reservoir)
{
    // For each incoming data sample
    for(int step = 0; step < stepsNumber; ++step) performSingleStep(reservoir, stepsNumber);
}

void biasedReservoirSamplingAlgorithm::updateFractionOfReservoirFilled(int currentReservoirSize)
{
    // Bias rate is equal to 1 over reservoir capacity, thus:
    this->fractionOfReservoirFilled = currentReservoirSize * biasRate;
}

void biasedReservoirSamplingAlgorithm::performSingleStep(std::vector<sample *> *reservoir, int stepNumber)
{
  int indexOfSampleToWriteOn;

  updateFractionOfReservoirFilled(reservoir->size());

  // Check if a sample should be removed
  if(fractionOfReservoirFilled >= ((double) rand() / (RAND_MAX)))
  {
    // If so randomly choose a sample
    indexOfSampleToWriteOn = (((double) rand() / (RAND_MAX)) * (reservoir->size() - 1));
  }
  else
  {
    // If not add new sample to reservoir
    indexOfSampleToWriteOn  = reservoir->size();
    parser->addDatumToContainer(reservoir);
  }

  // Add new data to the reservoir
  reader->getNextRawDatum(parser->buffer);
  parser->writeDatumOnPosition(reservoir, indexOfSampleToWriteOn);
}


