#ifndef PROGRESSIVEDISTRIBUTIONDATAREADER_H
#define PROGRESSIVEDISTRIBUTIONDATAREADER_H

#include "Reservoir_sampling/dataReader.h"
#include "Distributions/distribution.h"

#include <QVector>
#include <vector>

class progressiveDistributionDataReader : public dataReader
{
    public:

        progressiveDistributionDataReader(distribution *source, qreal progressionSize);

        void getNextRawDatum(void *target);
        void gatherAttributesData(void *attributes);
        bool hasMoreData();

        std::vector<std::string>* getAttributesOrder();

    protected:
        distribution *sourceDistribution;
        qreal progressionSize = 0.0;

        std::vector<std::string> attributesOrder;

};

#endif // PROGRESSIVEDISTRIBUTIONDATAREADER_H
