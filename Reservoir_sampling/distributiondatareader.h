#ifndef DISTRIBUTIONDATAREADER_H
#define DISTRIBUTIONDATAREADER_H

#include "Distributions/distribution.h"
#include "dataReader.h"

#include <QVector>

class distributionDataReader : public dataReader
{
    public:

        distributionDataReader(distribution* sourceDistribution);

        void getNextRawDatum(void *target);
        void gatherAttributesData(void *attributes);
        bool hasMoreData();

    protected:

        distribution *sourceDistribution;
};

#endif // DISTRIBUTIONDATAREADER_H
