#include "dullkernel.h"

#include "QDebug"

dullKernel::dullKernel()
{}

qreal dullKernel::getValue(QVector<qreal> *arguments)
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

    // If |x| <= 1 return 0.5
    if(qAbs(x) <= 1)
    {
        return 0.5;
    }

    return 0;
}
