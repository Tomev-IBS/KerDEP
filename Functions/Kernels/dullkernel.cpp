#include "dullkernel.h"
#include <cmath>
#include <iostream>

dullKernel::dullKernel()
{}

double dullKernel::getValue(vector<double> *arguments)
{
    // TODO: Consider using asserts of some kind

    // Check for nullpointer
    if(!arguments)
    {
        std::cout << "Nullpointer in gaussianProbabilityDensityFunction";
        return -1;
    }

    // Check if arguments length is equal to 1
    if(arguments->size() != 1)
    {
        // If not return -1 with wrong arguments size.
        std::cout << "Wrong arguments size.";
        return -1;
    }

    double x = arguments->at(0);

    // If |x| <= 1 return 0.5
    if(fabs(x) <= 1)
    {
        return 0.5;
    }

    return 0;
}
