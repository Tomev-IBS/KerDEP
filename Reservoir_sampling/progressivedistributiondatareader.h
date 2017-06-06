#ifndef PROGRESSIVEDISTRIBUTIONDATAREADER_H
#define PROGRESSIVEDISTRIBUTIONDATAREADER_H

#include "Reservoir_sampling/dataReader.h"
#include "Distributions/distribution.h"

#include <QVector>

class progressiveDistributionDataReader : public dataReader
{
    public:
        progressiveDistributionDataReader(distribution *source, qreal progressionSize);
        void getNextRawDatum(void *target);

    protected:
        distribution *sourceDistribution;
        qreal progressionSize;

};

#endif // PROGRESSIVEDISTRIBUTIONDATAREADER_H
