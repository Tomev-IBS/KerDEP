#include "normalkernel.h"

#include "QDebug"

normalKernel::normalKernel()
{}

qreal normalKernel::getValue(QVector<qreal>* arguments)
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
        return -1;
    }

    qreal x = arguments->at(0);

    qreal result = exp(- pow((x - mean), 2) / (2 * pow(standardDeviation, 2)));
    result /= (standardDeviation * sqrt(2 * M_PI));

    return result;
}
