#include "normaldistribution.h"

#include "QDebug"

normalDistribution::normalDistribution(int seed)
{
    generator = std::default_random_engine(seed);
    distribution = new std::normal_distribution<qreal>(0.0, 1.0);
}

qreal normalDistribution::getValue()
{
    return (*distribution)(generator);
}
