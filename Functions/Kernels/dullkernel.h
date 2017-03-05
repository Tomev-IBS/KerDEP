#ifndef DULLKERNELFUNCTION_H
#define DULLKERNELFUNCTION_H

#include "function.h"

class dullKernelFunction : public function
{
public:
    dullKernelFunction();

    qreal getValue(QVector<qreal>* arguments); // Only one argument should be passed
};

#endif // DULLKERNELFUNCTION_H
