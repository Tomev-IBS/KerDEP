#include "PKSamplingAlgorithm.h"

PKSamplingAlgorithm::PKSamplingAlgorithm(
        dataReader *reader, dataParser *parser, int reservoirMaxSize, int stepsNumber) :
        reservoirMaxSize(reservoirMaxSize)
{
  this->reader = reader;
  this->parser = parser;
  this->stepsNumber = stepsNumber;
}

void PKSamplingAlgorithm::fillReservoir(std::vector<std::shared_ptr<sample> > *reservoir)
{
  for(int step = 0; step < stepsNumber; ++step)
    performSingleStep(reservoir, step);
}

void PKSamplingAlgorithm::performSingleStep(std::vector<std::shared_ptr<sample> > *reservoir, int stepNumber)
{
  while(reservoir->size() >= reservoirMaxSize){
      int idxToDelete = (((double) rand() / (RAND_MAX)) * reservoir->size());
      reservoir->erase(reservoir->begin() + idxToDelete);
  }

  addDatumToReservoir(reservoir);
}

unsigned int PKSamplingAlgorithm::getReservoidMaxSize()
{
  return this->reservoirMaxSize;
}

void PKSamplingAlgorithm::changeReservoirMaxSize(unsigned int newMaxSize)
{
  reservoirMaxSize = newMaxSize;
}

void PKSamplingAlgorithm::addDatumToReservoir(std::vector<std::shared_ptr<sample>> *reservoir)
{
  reader->getNextRawDatum(parser->buffer);

  parser->addDatumToContainer(reservoir);

  parser->writeDatumOnPosition(reservoir, reservoir->size() - 1);
}
