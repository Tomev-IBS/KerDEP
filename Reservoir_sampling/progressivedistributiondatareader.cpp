#include "progressivedistributiondatareader.h"
#include "../groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

#include <string>
#include <unordered_map>
#include <cmath>

#include <QDebug>

progressiveDistributionDataReader::progressiveDistributionDataReader(distribution *source, double progressionSize,
                                                                     int delay, distribution *alternativeSource) :
    sourceDistribution(source), x_progression_size(progressionSize), _delay(delay) {
  _alternativeDistribution.reset(alternativeSource);
}

void progressiveDistributionDataReader::getNextRawDatum(void *target) {
  vector<double> *targetPtr = static_cast<vector<double> *>(target);
  targetPtr->clear();

  sourceDistribution->getValue(targetPtr);
  //_alternativeDistribution->getValue(targetPtr);

  //qDebug() << "Sample: " << (*targetPtr)[0] << ", " << (*targetPtr)[1];

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

  /*
  // Symmetric trimodal selection scheme.
  if((_currentIteration - 1) % 10 == 1 || (_currentIteration - 1) % 10 == 4 ||
     (_currentIteration - 1) % 10 == 7){
    (*targetPtr)[0] += bimodalMean;
  }

  if((_currentIteration - 1) % 10 == 2 || (_currentIteration - 1) % 10 == 5 ||
     (_currentIteration - 1) % 10 == 8){
    (*targetPtr)[0] += trimodalMean;
  }
  //*/

  /*
  // Asymmetric trimodal selection scheme.
  if((_currentIteration - 1) % 10 == 1 || (_currentIteration - 1) % 10 == 4 ||
     (_currentIteration - 1) % 10 == 7 || (_currentIteration - 1) % 10 == 8){
    (*targetPtr)[0] += bimodalMean;
  }

  if((_currentIteration - 1) % 10 == 2 || (_currentIteration - 1) % 10 == 5){
    (*targetPtr)[0] += trimodalMean;
  }
  //*/

  // 26 III 2020 article formula
  // Stops at 0.2 + 30 + 3 + 1 = 34.2 without offset. Set maxX = 40.
  switch(_currentIteration - 1) { // For exps with seed, remove later
    case 0:
      x_progression_size = 0.0001;
      break;
    //case 200: // Added for faster q test
      //x_progression_size = 0.1;
    case 2000:
    //case 5000:
      x_progression_size = 0.01;
      break;
    case 5000:
    //case 8000:
      x_progression_size = 0.001;
      break;
    case 8000:
    //case 11000:
      x_progression_size = 1;
      break;
    case 8001:
    //case 11001:
      x_progression_size = 0;
      break;
  }

  switch(_currentIteration - 1) { // For exps with seed, remove later
    case 0:
      y_progression_size = 0;
      break;
    case 7000:
      y_progression_size = 0.01;
      break;
    case 9000:
      y_progression_size = 0;
      break;
    case 11000:
      y_progression_size = 1;
      break;
    case 11001:
      y_progression_size = 0;
      break;
  }

  // SECOND DIMENSION
  /**/

  /*
  switch(_currentIteration - 1){ // For exps with seed, remove later
    case 0:
      x_progression_size = 0;
      break;
    case 2000:
      x_progression_size = 0.005;
      break;
    case 6000:
      x_progression_size = 0;
      break;
  }
  /**/

  /*
  switch(_currentIteration - 1){ // For exps with seed, remove later
    case 0:
      x_progression_size = 0;
      break;
    case 2000:
      x_progression_size = 1;
      break;
    case 2001:
      x_progression_size = 0;
      break;
    case 4000:
      x_progression_size = -1;
      break;
    case 4001:
      x_progression_size = 0;
      break;
    case 6000:
      x_progression_size = 1;
      break;
    case 6001:
      x_progression_size = 0;
      break;
  }
  /**/

  sourceDistribution->increaseMeans(0, 0);
  //sourceDistribution->increaseMeans(x_progression_size, 0);
  //_alternativeDistribution->increaseMeans(x_progression_size, 0);
  //sourceDistribution->increaseMeans(y_progression_size, 1);
  //_alternativeDistribution->increaseMeans(y_progression_size, 1);

  ++_currentIteration;
}

void progressiveDistributionDataReader::gatherAttributesData(void *attributes) {
  std::unordered_map<std::string, attributeData *> *attrs_ptr =
      static_cast<std::unordered_map<std::string, attributeData *> *>(attributes);

  // There are no attributes in distribution, just numbers;
  vector<double> *dummyValue = new vector<double>;
  sourceDistribution->getValue(dummyValue);

  for(size_t i = 0; i < dummyValue->size(); ++i) {
    std::string attrName = "Val" + std::to_string(i);

    attributesOrder.push_back(attrName);

    (*attrs_ptr)[attrName] = new numericalAttributeData(attrName);
  }

  return;
}

bool progressiveDistributionDataReader::hasMoreData() {
  // One can always generate more data from distribution.
  return true;
}

std::vector<std::__cxx11::string> *progressiveDistributionDataReader::getAttributesOrder() {
  return &attributesOrder;
}

void progressiveDistributionDataReader::setNewSource(distribution *source) {
  this->sourceDistribution = source;
}
