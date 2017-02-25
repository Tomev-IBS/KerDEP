#include "normaldistribution.h"

#include "QDebug"

normalDistribution::normalDistribution(int seed, qreal mean, qreal standardDeviation)
{
    generator = std::default_random_engine(seed);
    distribution = new std::normal_distribution<qreal>(mean, standardDeviation);
}

qreal normalDistribution::getValue()
{
    return (*distribution)(generator);
}
