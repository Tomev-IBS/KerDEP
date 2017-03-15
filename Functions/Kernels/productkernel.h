#ifndef PRODUCTKERNEL_H
#define PRODUCTKERNEL_H

#include "kernels.h"

class productKernel : public kernel
{
    public:
        productKernel(QVector<int>* idsOfKernels);

        qreal getValue(QVector<qreal>* arguments);

        qreal getU(){ return -1.0; }
        qreal getW(){ return -1.0; }

   private:

        QVector<kernelPtr> kernels;
};

#endif // PRODUCTKERNEL_H
