//
// Created by Tomev on 29.05.2017.
//

#include "basicReservoirSamplingAlgorithm.h"

basicReservoirSamplingAlgorithm::basicReservoirSamplingAlgorithm(
        dataReader *reader, dataParser *parser, int reservoirMaxSize, int stepsNumber) :
        reservoirMaxSize(reservoirMaxSize)
{
  this->reader = reader;
  this->parser = parser;
  this->stepsNumber = stepsNumber;
}

void basicReservoirSamplingAlgorithm::fillReservoir(std::vector<sample*> *reservoir)
{
  initializeReservoir(reservoir);

  for(int step = reservoirMaxSize+1; step < reservoirMaxSize; ++step)
    if(shouldDatumBeAdded(step)) addDatumOnRandomPosition (reservoir);
}

void basicReservoirSamplingAlgorithm::initializeReservoir(std::vector<sample*> *reservoir)
{
  while(reservoir->size() < reservoirMaxSize) addDatumToReservoir(reservoir);
}

void basicReservoirSamplingAlgorithm::performSingleStep(std::vector<sample*> *reservoir, int stepNumber)
{
  if(reservoir->size() < reservoirMaxSize) addDatumToReservoir(reservoir);
  else if(shouldDatumBeAdded(stepNumber)) addDatumOnRandomPosition(reservoir);
}

void basicReservoirSamplingAlgorithm::addDatumToReservoir(std::vector<sample *> *reservoir)
{
  reader->getNextRawDatum(parser->buffer);

  parser->addDatumToContainer(reservoir);

  parser->writeDatumOnPosition(reservoir, reservoir->size()-1);
}

void basicReservoirSamplingAlgorithm::addDatumOnRandomPosition(std::vector<sample *> *reservoir)
{
  int idxToDelete = (((double) rand() / (RAND_MAX)) * reservoir->size());

  parser->writeDatumOnPosition(reservoir, idxToDelete);
}

bool basicReservoirSamplingAlgorithm::shouldDatumBeAdded(int stepNumber)
{
  reader->getNextRawDatum(parser->buffer);

  double addChance = (double) reservoirMaxSize/stepNumber;

  return (addChance > ((double) rand() / (RAND_MAX)));
}
