#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H

#include "distribution.h"
#include "../Libraries/matrixoperationslibrary.h"
#include "random"

class normalDistribution : public distribution
{
    public:
        normalDistribution(int seed, vector<double>* means, vector<double>* stDevs);

        void getValue(vector<double>* result);
        void increaseMeans(double addend, int index=-1);

    private:      
        vector<double>* means;
        vector<double>* stDevs;
        vector<std::normal_distribution<double>*> distributions;
        vector<std::default_random_engine> generators;

        matrix A;

        void fillA(vector<vector<double> *> *covarianceMatrix);
};

#endif // NORMALDISTRIBUTION_H
