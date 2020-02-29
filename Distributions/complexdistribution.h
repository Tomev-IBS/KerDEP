#ifndef COMPLEXDISTRIBUTION_H
#define COMPLEXDISTRIBUTION_H

#include "distribution.h"

#include <memory>

class complexDistribution : public distribution
{
    public:

        complexDistribution(int seed, vector<std::shared_ptr<distribution>>* elementalDistributions, vector<double>* contributions);

        void getValue(vector<double>* result);
        void increaseMeans(double addend);

    private:
        vector<std::shared_ptr<distribution>> elementalDistributions;
        vector<double> contributions;

        const std::uniform_real_distribution<double> uniformDistribution;

        int randomizeDistributionIndex();
};

#endif // COMPLEXDISTRIBUTION_H
