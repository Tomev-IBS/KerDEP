#define _USE_MATH_DEFINES
#include <cmath>
#include "normalkernel.h"
#include "QDebug"

normalKernel::normalKernel()
{}

double normalKernel::getValue(vector<double>* arguments)
{
    // TODO: Consider using asserts of some kind

    // Check for nullpointer
    if(!arguments)
    {
        qDebug() << "Nullpointer in gaussianProbabilityDensityFunction";
        return -1;
    }

    // Check if arguments length is equal to 1
    if(arguments->size() != 1)
    {
        // If not return -1 with wrong arguments size.
        qDebug() << "Wrong arguments size.";
        qDebug() << "Got: " << arguments->size() << " Expected: " << 1;
        return -1;
    }

    double x = arguments->at(0);

    double result = exp(- pow((x - mean), 2) / (2 * pow(standardDeviation, 2)));
    result /= (standardDeviation * sqrt(2 * M_PI));

    return result;
}
