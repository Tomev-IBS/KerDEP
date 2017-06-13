#ifndef DISTRIBUTIONDATAPARSER_H
#define DISTRIBUTIONDATAPARSER_H

#include <QVector>

#include "dataParser.h"
#include "distributionDataSample.h"

class distributionDataParser: public dataParser
{
    public:
        distributionDataParser();
        void parseData(void *target);
        int addDatumToContainer(std::vector<sample*> *container);
        void writeDatumOnPosition(std::vector<sample*> *container, int position);
};

#endif // DISTRIBUTIONDATAPARSER_H
