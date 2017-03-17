#ifndef TRIANGLEKERNEL_H
#define TRIANGLEKERNEL_H

#include "kernel.h"

class triangleKernel : public kernel
{
public:
    triangleKernel();

    qreal getValue(QVector<qreal>* arguments); // Only one argument should be passed
    static qreal getW(){ return 2.0/3.0; }
    static qreal getU(){ return 1.0/6.0; }
};

#endif // TRIANGLEKERNEL_H
