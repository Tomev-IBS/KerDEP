#ifndef KERNEL_H
#define KERNEL_H

#include "function.h"

class kernel : public function
{
    public:
        virtual qreal getU() = 0;
        virtual qreal getW() = 0;
};

#endif // KERNEL_H
