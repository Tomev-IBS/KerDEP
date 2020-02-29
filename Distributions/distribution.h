#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>
#include <random>

using std::vector;

class distribution
{
    public:
        virtual void getValue(vector<double>* result) = 0;
        virtual void increaseMeans(double addend) = 0;

    protected:

        std::default_random_engine generator;
};

#endif // DISTRIBUTION_H
