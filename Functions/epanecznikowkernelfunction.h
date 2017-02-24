#ifndef EPANECZNIKOWKERNELFUNCTION_H
#define EPANECZNIKOWKERNELFUNCTION_H

#include "function.h"

class epanecznikowKernelFunction : public function
{
public:
    epanecznikowKernelFunction();

    qreal getValue(QVector<qreal>* arguments); // Only one argument should be passed
};

#endif // EPANECZNIKOWKERNELFUNCTION_H
