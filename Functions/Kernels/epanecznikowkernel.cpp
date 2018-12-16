#include "epanecznikowkernel.h"
#include "QDebug"
#include <cmath>

epanecznikowKernel::epanecznikowKernel()
{}

qreal epanecznikowKernel::getValue(QVector<qreal> *arguments)
{
    // TODO: Consider using asserts of some sort

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

    // For |x| <= 1 return 3/4(1-x^2)
    if(qAbs(x) <= 1)
    {
        return 3.0/4.0 * (1 - pow(x,2));
    }

    return 0;
}
