#ifndef DULLKERNEL_H
#define DULLKERNEL_H

#include "kernel.h"

class dullKernel : public kernel
{
public:
    dullKernel();

    qreal getValue(QVector<qreal>* arguments); // Only one argument should be passed
    qreal getW(){ return 1.0 / 3.0; }
    qreal getU(){ return 0.5; }
};

#endif // DULLKERNEL_H
