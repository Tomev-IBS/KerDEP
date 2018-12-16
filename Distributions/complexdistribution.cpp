#include "complexdistribution.h"

complexDistribution::complexDistribution(int seed, QVector<std::shared_ptr<distribution>> *elementalDistributions, QVector<qreal> *contributions) :
    uniformDistribution(0.0, 1.0)
{
    this->elementalDistributions = QVector<std::shared_ptr<distribution>>(*elementalDistributions);
    this->contributions = QVector<qreal>(*contributions);

    generator = std::default_random_engine(seed);
}

void complexDistribution::getValue(QVector<qreal> *result)
{
    int distributionIndex = randomizeDistributionIndex();

    elementalDistributions.at(distributionIndex)->getValue(result);
}

void complexDistribution::increaseMeans(qreal addend)
{
    foreach (std::shared_ptr<distribution> elementalDistribution, elementalDistributions)
    {
        elementalDistribution->increaseMeans(addend);
    }
}

int complexDistribution::randomizeDistributionIndex()
{
    std::uniform_real_distribution<qreal> uniformDistribution(0,1);

    qreal threshold = uniformDistribution(generator) * 100;
    qreal checker = contributions.at(0);

    int distributionIndex = 0;

    while(checker < threshold)
        checker += contributions.at(++distributionIndex);

    return distributionIndex;
}
