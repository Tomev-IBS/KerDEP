#ifndef COMPLEXDISTRIBUTION_H
#define COMPLEXDISTRIBUTION_H

#include "distribution.h"

#include <memory>

class complexDistribution : public distribution
{
    public:

        complexDistribution(int seed, vector<std::shared_ptr<distribution>>* elementalDistributions, vector<double>* contributions);

        void getValue(vector<double>* result);
        virtual void increaseMeans(double addend, int index=-1);
        void setMeans(double newMean, int index=-1);

    protected:
        vector<std::shared_ptr<distribution>> elementalDistributions;
        vector<double> contributions;

        const std::uniform_real_distribution<double> uniformDistribution;

        virtual int randomizeDistributionIndex();
};

#endif // COMPLEXDISTRIBUTION_H
