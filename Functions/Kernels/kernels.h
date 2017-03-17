#ifndef KERNELS_H
#define KERNELS_H

#include "normalkernel.h"
#include "trianglekernel.h"
#include "epanecznikowkernel.h"
#include "dullkernel.h"

enum kernelsIDs
{
    NORMAL          =   0,
    TRIANGLE        =   1,
    EPANECZNIKOW    =   2,
    DULL            =   3
};

#endif // KERNELS_H
