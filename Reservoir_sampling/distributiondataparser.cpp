#include <QDebug>

#include "distributiondataparser.h"
#include "../groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

distributionDataParser::distributionDataParser(std::unordered_map<std::string, attributeData *> *attributesData)
{
    buffer = new QVector<qreal>();
    this->attributesData = attributesData;
}

void distributionDataParser::parseData(void *target)
{
  QVector<qreal> *data            = static_cast<QVector<qreal>*>(buffer);
  distributionDataSample *sample  = static_cast<distributionDataSample*>(target);

  sample->attributesValues.clear();

  for(int i = 0; i < data->size(); ++i)
    sample->attributesValues[attributesOrder->at(i)] = std::to_string(data->at(i));

  sample->attributesData = attributesData;

  updateAttributesData(sample);
}

int distributionDataParser::addDatumToContainer(std::vector<std::shared_ptr<sample>> *container)
{
    container->push_back(std::shared_ptr<sample>(new distributionDataSample()));

    return container->size();
}

void distributionDataParser::writeDatumOnPosition(std::vector<std::shared_ptr<sample>> *container, int position)
{
  parseData(container->at(position).get());
}

void distributionDataParser::setAttributesOrder(std::vector<std::string> *attributesOrder)
{
  this->attributesOrder = attributesOrder;
}

void distributionDataParser::updateAttributesData(distributionDataSample *sample)
{
  for(auto kv : sample->attributesValues)
    {
      // For now only consider numerical data
      if(attributesData->at(kv.first)->getType() == "numerical")
      {
        numericalAttributeData *numAttribute = static_cast<numericalAttributeData*>(attributesData->at(kv.first));
        double value = stod(kv.second);

        numAttribute->setMaximalValue(value);
        numAttribute->setMinimalValue(value);

        numAttribute->attributeOccured();
      }
    }
}
