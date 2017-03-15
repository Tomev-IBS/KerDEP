#include "productkernel.h"

#include "QDebug"

productKernel::productKernel(QVector<int> *idsOfKernels)
{
    // Check for null pointer
    if(idsOfKernels == NULL)
    {
        // Log and return if occured
        qDebug() << "Kernels IDs vector pointer is null.";
        return;
    }

    // Check if vector size is lower than 1
    if(idsOfKernels->size() < 1)
    {
        // Log and return if so
        qDebug() << "Kernels IDs vector has no entries.";
        return;
    }

    // For each entry in this vector create a kernel and add it to list
    foreach(const int kernelID, *idsOfKernels)
    {
        switch (kernelID)
        {
            case NORMAL:
                kernels.append(kernelPtr(new normalKernel(0, 1.0)));
            break;
            case TRIANGLE:
                kernels.append(kernelPtr(new triangleKernel()));
            break;
            case EPANECZNIKOW:
                kernels.append(kernelPtr(new epanecznikowKernel()));
            break;
            case DULL:
            default:
                kernels.append(kernelPtr(new dullKernel()));
            break;
        }
    }
}

qreal productKernel::getValue(QVector<qreal> *arguments)
{
    // Check for null in arguments pointer
    if(arguments == NULL)
    {
        qDebug() << "Null arguments pointer in kerenels getValue.";
        return -1.0;
    }

    // Check if arguments size is equal to kernels size
    if(arguments->size() != kernels.size())
    {
        qDebug() << "Kernel's dimension and argument's dimension aren't equal.";
        return -1.0;
    }

    // Initiate result with value from first kernel

    QVector<qreal> *tempValueHolder = new QVector<qreal>();
    qreal result = 1.0;

    // For each other kernel multiply result
    for(int kernelNum = 0; kernelNum < kernels.size(); ++kernelNum)
    {
        tempValueHolder->clear();
        tempValueHolder->append(arguments->at(kernelNum));

        result *= kernels.at(kernelNum)->getValue(arguments);
    }


    return result;
}
