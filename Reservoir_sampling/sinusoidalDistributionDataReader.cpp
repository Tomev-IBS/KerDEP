#include "sinusoidalDistributionDataReader.h"
#include "../groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

#include <string>
#include <unordered_map>
#include <cmath>

#include <QDebug>

sinusoidalDistributionDataReader::sinusoidalDistributionDataReader(distribution *source, double periodsNumber,
                                                                     int stepsNumber, distribution *alternativeSource, const double &d2_speed_multiplier,
                                                                     int sineStart, int sineStop) :
    sourceDistribution(source), _stepsNumber(stepsNumber), _periodsNumber(periodsNumber), _d2_speed_multiplier_(d2_speed_multiplier),
    _sineStart(sineStart), _sineStop(sineStop){
  _alternativeDistribution.reset(alternativeSource);
  _sineSteps = _sineStop - _sineStart;
}

void sinusoidalDistributionDataReader::getNextRawDatum(void *target) {
  vector<double> *targetPtr = static_cast<vector<double> *>(target);
  targetPtr->clear();

  sourceDistribution->getValue(targetPtr);
  if(_currentIteration > _sineStart && _currentIteration < _sineStop && _sineSteps > 0) {
    double newMean = sin(_periodsNumber * (_currentIteration - _sineStart) * 2 * M_PI / _sineSteps);

    sourceDistribution->setMeans(newMean, 0);
    _alternativeDistribution->setMeans(newMean, 0);

    if(fabs(_d2_speed_multiplier_) > 1e-15) {
      sourceDistribution->setMeans(newMean * _d2_speed_multiplier_, 1);
      _alternativeDistribution->setMeans(newMean * _d2_speed_multiplier_, 1);
    }
  }



  ++_currentIteration;
}

void sinusoidalDistributionDataReader::gatherAttributesData(void *attributes) {
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

bool sinusoidalDistributionDataReader::hasMoreData() {
  // One can always generate more data from distribution.
  return true;
}

std::vector<std::__cxx11::string> *sinusoidalDistributionDataReader::getAttributesOrder() {
  return &attributesOrder;
}

void sinusoidalDistributionDataReader::setNewSource(distribution *source) {
  this->sourceDistribution = source;
}
