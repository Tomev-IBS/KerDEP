#ifndef KERNEL_H
#define KERNEL_H

#include "../function.h"
#include "memory"

class kernel : public function
{
    public:
        static qreal getU(){return 0.0;}
        static qreal getW(){return 0.0;}
};

typedef std::shared_ptr<kernel> kernelPtr;

#endif // KERNEL_H
