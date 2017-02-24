#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "gaussianprobabilitydensityfunction.h"
#include "trianglekernelfunction.h"
#include "epanecznikowkernelfunction.h"
#include "dullkernelfunction.h"

enum functionsIDs
{
    NORMAL      = 0,
    TRIANGLE    = 1,
    EPACZNIKOW  = 2,
    DULL        = 3,
};

#endif // FUNCTIONS_H
