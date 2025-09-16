//
// Created by tomev on 1/4/2023.
//

#ifndef DEDSTA_SINUSOIDALDISTRIBUTIONDATAREADER_H
#define DEDSTA_SINUSOIDALDISTRIBUTIONDATAREADER_H

#include "Reservoir_sampling/dataReader.h"
#include "Distributions/distribution.h"

#include <QVector>
#include <vector>
#include <memory>

class sinusoidalDistributionDataReader : public dataReader
{
    public:

      sinusoidalDistributionDataReader(distribution *source, qreal periodsNumber, int stepsNumber,
                                        distribution* alternativeSource, const double &d2_speed_multiplier = 0,
                                        int sineStart=1000, int sineStop=9000);

        void getNextRawDatum(void *target);
        void gatherAttributesData(void *attributes);
        bool hasMoreData();

        std::vector<std::string>* getAttributesOrder();

        void setNewSource(distribution *source);

    protected:
        distribution *sourceDistribution;
        qreal _periodsNumber = 0.0;
        int _currentIteration = 1; // PK always starts from 1.
        int _stepsNumber;
        float _lastMeansIncrease = 0;
    double _d2_speed_multiplier_;

      int _sineStart = 1000;
      int _sineStop = 9000;
      int _sineSteps = 1;

        std::vector<std::string> attributesOrder;
        std::shared_ptr<distribution> _alternativeDistribution; // For non-random multimodals
};


#endif //DEDSTA_SINUSOIDALDISTRIBUTIONDATAREADER_H
