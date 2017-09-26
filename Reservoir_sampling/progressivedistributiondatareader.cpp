#include "progressivedistributiondatareader.h"
#include "../groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

#include <string>
#include <unordered_map>

#include <QDebug>

progressiveDistributionDataReader::progressiveDistributionDataReader(distribution *source, qreal progressionSize) :
    sourceDistribution(source), progressionSize(progressionSize)
{}

void progressiveDistributionDataReader::getNextRawDatum(void *target)
{
    QVector<qreal>* targetPtr = static_cast<QVector<qreal>*>(target);
    targetPtr->clear();

    sourceDistribution->getValue(targetPtr);

    sourceDistribution->increaseMeans(progressionSize);
}

void progressiveDistributionDataReader::gatherAttributesData(void *attributes)
{
  std::unordered_map<std::string, attributeData*>* attrs_ptr =
  static_cast<std::unordered_map<std::string, attributeData*>*>(attributes);

  // There are no attributes in distribution, just numbers;
  QVector<qreal>* dummyValue = new QVector<qreal>;
  sourceDistribution->getValue(dummyValue);

  for(int i = 0; i < dummyValue->size(); ++i)
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
