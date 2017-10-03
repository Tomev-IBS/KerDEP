#ifndef DISTRIBUTIONDATAPARSER_H
#define DISTRIBUTIONDATAPARSER_H

#include <QVector>
#include <unordered_map>
#include <vector>
#include <string>

#include "dataParser.h"
#include "distributionDataSample.h"
#include "../groupingThread/kMedoidsAlgorithm/attributeData.h"

class distributionDataParser: public dataParser
{
    public:
    distributionDataParser(std::unordered_map<std::string, attributeData*> *attributesData);
    void parseData(void *target);
    int addDatumToContainer(std::vector<std::shared_ptr<sample>> *container);
    void writeDatumOnPosition(std::vector<std::shared_ptr<sample>> *container, int position);

    void setAttributesOrder(std::vector<std::string> *attributesOrder);

  protected:
    std::vector<std::string> *attributesOrder;
    std::unordered_map<std::string, attributeData*> *attributesData;

    void updateAttributesData(distributionDataSample *sample);
};

#endif // DISTRIBUTIONDATAPARSER_H
