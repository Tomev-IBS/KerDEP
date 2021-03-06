#include "trianglekernel.h"

#include <cmath>

triangleKernel::triangleKernel()
{}

double triangleKernel::getValue(vector<double> *arguments)
{
    // TODO: Consider using asserts of some sort

    // Check for nullpointer
    if(!arguments)
    {
        //qDebug() << "Nullpointer in gaussianProbabilityDensityFunction";
        return -1;
    }

    // Check if arguments length is equal to 1
    if(arguments->size() != 1)
    {
        // If not return -1 with wrong arguments size.
        //qDebug() << "Wrong arguments size.";
        return -1;
    }

    double x = fabs(arguments->at(0));

    // For |x| < 1 return 1 - |x|
    if(x < 1)
    {
        return 1.0 - x;
    }

    // Otherwise return 0

    return 0;
}
