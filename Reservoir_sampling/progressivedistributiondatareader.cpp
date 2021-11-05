#include "progressivedistributiondatareader.h"
#include "../groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

#include <string>
#include <unordered_map>
#include <cmath>

#include <QDebug>

progressiveDistributionDataReader::progressiveDistributionDataReader(distribution *source, double progressionSize,
                                                                     int delay, distribution *alternativeSource,
                                                                     const double &d2_speed_multiplier) :
    sourceDistribution(source), x_progression_size(progressionSize), _delay(delay),
    d2_speed_multiplier_(d2_speed_multiplier) {
  _alternativeDistribution.reset(alternativeSource);
}

void progressiveDistributionDataReader::getNextRawDatum(void *target) {
  vector<double> *targetPtr = static_cast<vector<double> *>(target);
  targetPtr->clear();

  sourceDistribution->getValue(targetPtr);

  int bimodalMean = 5;
  int trimodalMean = -5;

  // 26 III 2020 article formula
  // Stops at 0.5 + 6 + 3 + 1 = 10.5 without offset. Set maxX = 15.
  // 0.2 + 30 + 3 + 1 = 34.2 // Dla klasycznej. maxX = 38
  switch(_currentIteration - 1) { // For exps with seed, remove later
    case 0:
      //x_progression_size = 0.0005;
      x_progression_size = 0.001;
      break;
    //case 200: // Added for faster q test
      //x_progression_size = 0.1;
    case 2000:
    //case 5000: // 2D
      x_progression_size = 0.001; // Klasyczna + 30
      // x_progression_size = 0.005; // Spowolniona + 15
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

  sourceDistribution->increaseMeans(x_progression_size, 0);
  _alternativeDistribution->increaseMeans(x_progression_size, 0);
  if(fabs(d2_speed_multiplier_) > 1e-15) {
    sourceDistribution->increaseMeans(x_progression_size * d2_speed_multiplier_, 1);
    _alternativeDistribution->increaseMeans(x_progression_size * d2_speed_multiplier_, 1);
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
