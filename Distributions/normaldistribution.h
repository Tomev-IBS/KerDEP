#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H

#include "distribution.h"
#include "../Libraries/matrixoperationslibrary.h"
#include "random"

class normalDistribution : public distribution
{
    public:
        normalDistribution(int seed, vector<double>* means, vector<double>* stDevs, double maxMean);

        void getValue(vector<double>* result);
        void increaseMeans(double addend);

    private:      
        vector<double>* means;
        vector<double>* stDevs;
        vector<std::normal_distribution<double>*> distributions;

        matrix A;

        double _maxMean = 0;

        void fillA(vector<vector<double> *> *covarianceMatrix);
};

#endif // NORMALDISTRIBUTION_H
