#ifndef KERNEL_H
#define KERNEL_H

#include "../function.h"
#include "memory"

class kernel : public function
{
    public:
        static double getU(){return 0.0;}
        static double getW(){return 0.0;}
};

typedef std::shared_ptr<kernel> kernelPtr;

#endif // KERNEL_H
