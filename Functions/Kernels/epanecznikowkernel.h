#ifndef EPANECZNIKOWKERNEL_H
#define EPANECZNIKOWKERNEL_H

#include "kernel.h"

class epanecznikowKernel : public kernel
{
public:
    epanecznikowKernel();

    qreal getValue(QVector<qreal>* arguments); // Only one argument should be passed
    static qreal getW(){ return 0.2; }
    static qreal getU(){ return 0.6; }
};

#endif // EPANECZNIKOWKERNEL_H
