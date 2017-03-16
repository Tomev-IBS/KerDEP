#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H

#include "distribution.h"

#include "random"

class normalDistribution : public distribution
{
    public:
        normalDistribution(int seed);

        qreal getValue();

    private:
        std::normal_distribution<qreal>* distribution;

};

#endif // NORMALDISTRIBUTION_H
