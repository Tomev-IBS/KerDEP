#ifndef PROGRESSIVEDISTRIBUTIONDATAREADER_H
#define PROGRESSIVEDISTRIBUTIONDATAREADER_H

#include "Reservoir_sampling/dataReader.h"
#include "Distributions/distribution.h"

#include <QVector>
#include <vector>

class progressiveDistributionDataReader : public dataReader
{
    public:

        progressiveDistributionDataReader(distribution *source, qreal progressionSize, int delay, bool shouldJump = false);

        void getNextRawDatum(void *target);
        void gatherAttributesData(void *attributes);
        bool hasMoreData();

        std::vector<std::string>* getAttributesOrder();

        void setNewSource(distribution *source);

    protected:
        distribution *sourceDistribution;
        qreal progressionSize = 0.0;
        int _delay = 0;
        int _currentIteration = 1; // PK always starts from 1.
        bool _hasJumped = false;
        bool _shouldJump = false;

        std::vector<std::string> attributesOrder;



};

#endif // PROGRESSIVEDISTRIBUTIONDATAREADER_H
