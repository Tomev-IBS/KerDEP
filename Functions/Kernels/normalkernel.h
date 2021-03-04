#ifndef NORMALKERNEL_H
#define NORMALKERNEL_H

#include "kernel.h"
#include <cmath>
#define _USE_MATH_DEFINES

class normalKernel : public kernel
{
    public:
        normalKernel();

        double getValue(vector<double>* arguments); // Only one argument should be passed
        static double getW(){ return 1.0 / (2.0 * 3.14); }
        static double getU(){ return 1.0; }

    private:
        double  mean = 0.0,
                standardDeviation = 1.0;

};

#endif // NORMALKERNEL_H
