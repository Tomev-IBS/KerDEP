#include "progressivedistributiondatareader.h"
#include "../groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

#include <string>
#include <unordered_map>
#include <cmath>

#include <QDebug>

progressiveDistributionDataReader::progressiveDistributionDataReader(distribution *source, double progressionSize,
                                                                     int delay, distribution *alternativeSource,
                                                                     const double &d2_speed_multiplier,
                                                                     const double &d3_speed_multiplier) :
    sourceDistribution(source), x_progression_size(progressionSize), _delay(delay),
    d2_speed_multiplier_(d2_speed_multiplier), d3_speed_multiplier_(d3_speed_multiplier) {
  _alternativeDistribution.reset(alternativeSource);
}

void progressiveDistributionDataReader::getNextRawDatum(void *target) {
  vector<double> *targetPtr = static_cast<vector<double> *>(target);
  targetPtr->clear();

  sourceDistribution->getValue(targetPtr);

  /*
  switch(_currentIteration - 1) { // For exps with seed, remove later
    case 0:
      x_progression_size = 0;
      //x_progression_size = 0.0005;
      break;
      // case 300: // Test
    case 1000: // 1D
      //case 5000: // 2D
      x_progression_size = 0.005; // + 40
      //x_progression_size = 0.05; // + 400
      // x_progression_size = 0.1; // + 800
      break;
    case 9000:
      //case 11001: // 2D
      x_progression_size = 0;
      break;
  }
  //*/
  // 26 III 2020 article formula
  // Stops at 0.5 + 6 + 3 + 1 = 10.5 without offset. Set maxX = 15.
  // 0.2 + 30 + 3 + 1 = 34.2 // Dla klasycznej. maxX = 38
  //*
  switch(_currentIteration - 1) { // For exps with seed, remove later
    case 0:
      x_progression_size = 0.001;
      //x_progression_size = 0.0005;
      break;
    // case 300: // Test
    case 2000: // 1D
    //case 5000: // 2D
      // x_progression_size = 0.1; // Test
      x_progression_size = 0.01; // Klasyczna + 30
      //x_progression_size = 0.005; // Spowolniona + 15
      // x_progression_size = 0.002; //  Leniwa + 6
      break;
    case 6000:
    //case 8000: // 2D
      //x_progression_size = 0.005;
      x_progression_size = 0;
      break;
    case 8000:
    //case 11000: // 2D
      x_progression_size = 1;
      break;
    case 8001:
    //case 11001: // 2D
      x_progression_size = 0;
      break;
  }
  //*/

  // Ziegler-Nichols 0-1-0-1 jumps
  /*
  x_progression_size = 0;
  switch(_currentIteration - 1) { // For exps with seed, remove later
    case 2500:
      x_progression_size = 1;
      break;
    case 2501:
      x_progression_size = 0;
      break;
    case 5000:
      x_progression_size = -1;
      break;
    case 5001:
      x_progression_size = 0;
      break;
    case 7500:
      x_progression_size = 1;
      break;
    case 7501:
      x_progression_size = 0;
      break;
  }
  //*/

  sourceDistribution->increaseMeans(x_progression_size, 0);
  _alternativeDistribution->increaseMeans(x_progression_size, 0);
  if(fabs(d2_speed_multiplier_) > 1e-15) {
    sourceDistribution->increaseMeans(x_progression_size * d2_speed_multiplier_, 1);
    _alternativeDistribution->increaseMeans(x_progression_size * d2_speed_multiplier_, 1);
  }

  if(fabs(d3_speed_multiplier_) > 1e-15) {
    sourceDistribution->increaseMeans(x_progression_size * d3_speed_multiplier_, 2);
    _alternativeDistribution->increaseMeans(x_progression_size * d3_speed_multiplier_, 2);
  }

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
