#include "basicReservoirSamplingAlgorithm.h"

basicReservoirSamplingAlgorithm::basicReservoirSamplingAlgorithm(
        dataReader *reader, dataParser *parser, int reservoirMaxSize, int stepsNumber) :
        reservoirMaxSize(reservoirMaxSize)
{
  this->reader = reader;
  this->parser = parser;
  this->stepsNumber = stepsNumber;
}

void basicReservoirSamplingAlgorithm::fillReservoir(std::vector<std::shared_ptr<sample> > *reservoir)
{
  for(int step = 0; step < stepsNumber; ++step)
    performSingleStep(reservoir, step);
}

void basicReservoirSamplingAlgorithm::performSingleStep(std::vector<std::shared_ptr<sample> > *reservoir, int stepNumber)
{
  if(reservoir->size() < reservoirMaxSize) addDatumToReservoir(reservoir);
  else if(shouldDatumBeAdded(stepNumber)) addDatumOnRandomPosition(reservoir);
}

unsigned int basicReservoirSamplingAlgorithm::getReservoidMaxSize()
{
  return this->reservoirMaxSize;
}

void basicReservoirSamplingAlgorithm::changeReservoirMaxSize(unsigned int newMaxSize)
{
  reservoirMaxSize = newMaxSize;
}

std::vector<std::string>* basicReservoirSamplingAlgorithm::getAttributesList()
{
  return reader->getAttributesOrder();
}

void basicReservoirSamplingAlgorithm::addDatumToReservoir(std::vector<std::shared_ptr<sample>> *reservoir)
{
  reader->getNextRawDatum(parser->buffer);

  parser->addDatumToContainer(reservoir);

  parser->writeDatumOnPosition(reservoir, reservoir->size() - 1);
}

void basicReservoirSamplingAlgorithm::addDatumOnRandomPosition(std::vector<std::shared_ptr<sample>> *reservoir)
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
