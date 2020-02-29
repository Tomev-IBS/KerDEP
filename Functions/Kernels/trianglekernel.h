#ifndef TRIANGLEKERNEL_H
#define TRIANGLEKERNEL_H

#include "kernel.h"

class triangleKernel : public kernel
{
public:
    triangleKernel();

    double getValue(vector<double>* arguments); // Only one argument should be passed
    static double getW(){ return 2.0/3.0; }
    static double getU(){ return 1.0/6.0; }
};

#endif // TRIANGLEKERNEL_H
