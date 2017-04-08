#include "complexdistribution.h"

complexDistribution::complexDistribution(int seed, QVector<distribution *> *elementalDistributions, QVector<qreal> *contributions) :
    elementalDistributions(elementalDistributions), contributions(contributions), uniformDistribution(0.0, 1.0)
{
    generator = std::default_random_engine(seed);
}

void complexDistribution::getValue(QVector<qreal> *result)
{
    int distributionIndex = randomizeDistributionIndex();

    elementalDistributions->at(distributionIndex)->getValue(result);
}

int complexDistribution::randomizeDistributionIndex()
{
    std::uniform_real_distribution<qreal> uniformDistribution(0,1);

    qreal threshold = uniformDistribution(generator) * 100,
          checker = contributions->at(0);

    int distributionIndex = 0;

    while(checker < threshold)
        checker += contributions->at(++distributionIndex);

    return distributionIndex;
}
