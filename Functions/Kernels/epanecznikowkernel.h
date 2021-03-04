#ifndef EPANECZNIKOWKERNEL_H
#define EPANECZNIKOWKERNEL_H

#include "kernel.h"

class epanecznikowKernel : public kernel
{
public:
    epanecznikowKernel();

    double getValue(vector<double>* arguments); // Only one argument should be passed
    static double getW(){ return 0.2; }
    static double getU(){ return 0.6; }
};

#endif // EPANECZNIKOWKERNEL_H
