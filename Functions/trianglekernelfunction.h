#ifndef TRIANGLEKERNELFUNCTION_H
#define TRIANGLEKERNELFUNCTION_H

#include "function.h"

class triangleKernelFunction : public function
{
public:
    triangleKernelFunction();

    qreal getValue(QVector<qreal>* arguments); // Only one argument should be passed
};

#endif // TRIANGLEKERNELFUNCTION_H
