#include "progressivedistributiondatareader.h"
#include "../groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

#include <string>
#include <unordered_map>
#include <cmath>

#include <QDebug>

progressiveDistributionDataReader::progressiveDistributionDataReader(distribution *source, double progressionSize, int delay, distribution *alternativeSource) :
    sourceDistribution(source), progressionSize(progressionSize), _delay(delay)

{
  _alternativeDistribution.reset(alternativeSource);
}

void progressiveDistributionDataReader::getNextRawDatum(void *target)
{
    vector<double>* targetPtr = static_cast<vector<double>*>(target);
    targetPtr->clear();

    //sourceDistribution->getValue(targetPtr);
    _alternativeDistribution->getValue(targetPtr);

    /*
    for(auto attributeName : attributesOrder)
      (*attrs_ptr)[attributeName];
    */

    int bimodalMean = 5;
    int trimodalMean = -5;

    /*
    // Bimodal selection scheme.
    if((_currentIteration - 1) % 10 == 1 || (_currentIteration - 1) % 10 == 4 ||
       (_currentIteration - 1) % 10 == 7 || (_currentIteration - 1) % 10 == 9){
      (*targetPtr)[0] += bimodalMean;
    }
    //*/

    //*
    // Trimodal selection scheme.
    if((_currentIteration - 1) % 10 == 1 || (_currentIteration - 1) % 10 == 4 ||
       (_currentIteration - 1) % 10 == 7){
      (*targetPtr)[0] += bimodalMean;
    }

    if((_currentIteration - 1) % 10 == 2 || (_currentIteration - 1) % 10 == 5 ||
       (_currentIteration - 1) % 10 == 8){
      (*targetPtr)[0] += trimodalMean;
    }
    //*/

    // 26 III 2020 article formula
    // Stops at 0.2 + 30 + 3 + 1 = 34.2 without offset. Set maxX = 40.
    /*
    switch(_currentIteration - 1){ // For exps with seed, remove later
      case 0:
        progressionSize = 0.0001;
        break;
      case 2000:
        progressionSize = 0.01;
        break;
      case 5000:
        progressionSize = 0.001;
        break;
      case 8000:
        progressionSize = 1;
        break;
      case 8001:
        progressionSize = 0;
        break;
    }
    /**/

    //*
    switch(_currentIteration - 1){ // For exps with seed, remove later
      case 0:
        progressionSize = 0;
        break;
      case 2000:
        progressionSize = 0.005;
        break;
      case 6000:
        progressionSize = 0;
        break;
    }
    /**/

    /*
    switch(_currentIteration - 1){ // For exps with seed, remove later
      case 0:
        progressionSize = 0;
        break;
      case 2000:
        progressionSize = 1;
        break;
      case 2001:
        progressionSize = 0;
        break;
      case 4000:
        progressionSize = -1;
        break;
      case 4001:
        progressionSize = 0;
        break;
      case 6000:
        progressionSize = 1;
        break;
      case 6001:
        progressionSize = 0;
        break;
    }
    /**/

    sourceDistribution->increaseMeans(progressionSize);
    _alternativeDistribution->increaseMeans(progressionSize);

    ++_currentIteration;
}

void progressiveDistributionDataReader::gatherAttributesData(void *attributes)
{
  std::unordered_map<std::string, attributeData*>* attrs_ptr =
  static_cast<std::unordered_map<std::string, attributeData*>*>(attributes);

  // There are no attributes in distribution, just numbers;
  vector<double>* dummyValue = new vector<double>;
  sourceDistribution->getValue(dummyValue);

  for(size_t i = 0; i < dummyValue->size(); ++i)
  {
    std::string attrName = "Val"+std::to_string(i);

    attributesOrder.push_back(attrName);

    (*attrs_ptr)[attrName] = new numericalAttributeData(attrName);
  }

  return;
}

bool progressiveDistributionDataReader::hasMoreData()
{
    // One can always generate more data from distribution.
  return true;
}

std::vector<std::__cxx11::string> *progressiveDistributionDataReader::getAttributesOrder()
{
  return &attributesOrder;
}

void progressiveDistributionDataReader::setNewSource(distribution *source)
{
  this->sourceDistribution = source;
}
