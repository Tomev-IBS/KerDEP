#include "complexdistribution.h"

complexDistribution::complexDistribution(int seed, vector<std::shared_ptr<distribution>> *elementalDistributions, vector<double> *contributions) :
    uniformDistribution(0.0, 1.0)
{
    this->elementalDistributions = vector<std::shared_ptr<distribution>>(*elementalDistributions);
    this->contributions = vector<double>(*contributions);

    generator = std::default_random_engine(seed);
}

void complexDistribution::getValue(vector<double> *result)
{
    int distributionIndex = randomizeDistributionIndex();

    elementalDistributions.at(distributionIndex)->getValue(result);
}

void complexDistribution::increaseMeans(double addend)
{
    for(std::shared_ptr<distribution> elementalDistribution : elementalDistributions)
    {
        elementalDistribution->increaseMeans(addend);
    }
}

int complexDistribution::randomizeDistributionIndex()
{
    std::uniform_real_distribution<double> uniformDistribution(0,1);

    double threshold = uniformDistribution(generator) * 100;
    double checker = contributions.at(0);

    int distributionIndex = 0;

    while(checker < threshold)
        checker += contributions.at(++distributionIndex);

    return distributionIndex;
}
