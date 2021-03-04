#ifndef DULLKERNEL_H
#define DULLKERNEL_H

#include "kernel.h"

class dullKernel : public kernel
{
public:
    dullKernel();

    double getValue(vector<double>* arguments); // Only one argument should be passed
    static double getW(){ return 1.0 / 3.0; }
    static double getU(){ return 0.5; }
};

#endif // DULLKERNEL_H
