#ifndef NORMALKERNEL_H
#define NORMALKERNEL_H

#include "kernel.h"

class normalKernel : public kernel
{
    public:
        normalKernel();

        qreal getValue(QVector<qreal>* arguments); // Only one argument should be passed
        static qreal getW(){ return 1.0 / (2.0 * M_PI); }
        static qreal getU(){ return 1.0; }

    private:
        qreal   mean = 0.0,
                standardDeviation = 1.0;

};

#endif // NORMALKERNEL_H
