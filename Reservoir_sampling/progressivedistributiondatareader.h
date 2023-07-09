#ifndef PROGRESSIVEDISTRIBUTIONDATAREADER_H
#define PROGRESSIVEDISTRIBUTIONDATAREADER_H

#include "Reservoir_sampling/dataReader.h"
#include "Distributions/distribution.h"

#include <QVector>
#include <vector>
#include <memory>

class progressiveDistributionDataReader : public dataReader
{
    public:

        progressiveDistributionDataReader(distribution *source, qreal progressionSize, int delay,
                                          distribution* alternativeSource, const double &d2_speed_multiplier = 0,
                                          const double &d3_speed_multiplier = 0);

        void getNextRawDatum(void *target);
        void gatherAttributesData(void *attributes);
        bool hasMoreData();

        std::vector<std::string>* getAttributesOrder();

        void setNewSource(distribution *source);

    protected:
        distribution *sourceDistribution;
        qreal x_progression_size = 0.0;
        int _currentIteration = 1; // PK always starts from 1.
        double d2_speed_multiplier_;
        double d3_speed_multiplier_;
        int _delay;

        std::vector<std::string> attributesOrder;
        std::shared_ptr<distribution> _alternativeDistribution; // For non-random multimodals
};

#endif // PROGRESSIVEDISTRIBUTIONDATAREADER_H
