#ifndef DISTRIBUTIONDATAPARSER_H
#define DISTRIBUTIONDATAPARSER_H

#include "dataParser.h"

#include <QVector>

class distributionDataParser: public dataParser
{
    public:
        distributionDataParser();
        void parseData(void *target);
        int addDatumToContainer(void *container);
        void writeDatumOnPosition(void *container, int position);
};

#endif // DISTRIBUTIONDATAPARSER_H
