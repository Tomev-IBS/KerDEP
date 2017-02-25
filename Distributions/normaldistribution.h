#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H

#include "distribution.h"

#include "random"

class normalDistribution : public distribution
{
    public:
        normalDistribution(int seed, qreal mean = 0.0, qreal standardDeviation = 1.0);

        qreal getValue();

    private:
        std::normal_distribution<qreal>* distribution;

};

#endif // NORMALDISTRIBUTION_H
