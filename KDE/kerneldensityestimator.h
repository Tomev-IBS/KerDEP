#ifndef KERNELDENSITYESTIMATOR_H
#define KERNELDENSITYESTIMATOR_H

#include "../Functions/Kernels/kernels.h"
#include "../Distributions/distribution.h"

#include <QObject>

enum estimatorsKernelsType
{
    PRODUCT =   0,
    RADIAL  =   1
};

class kernelDensityEstimator
{
    public:
        kernelDensityEstimator(QVector<qreal>* samples, QVector<qreal>* smoothingParameter, int kernelType, QVector<int>* kernelsIDs);

        qreal getValue(QVector<qreal>* x);

    private:

        int kernelType;

        QVector<qreal>* samples;
        QVector<kernelPtr> kernels;
        QVector<qreal>* smoothingParameters;

        qreal getProductKernelValue(QVector<qreal>* x);
        qreal getRadialKernelValue(QVector<qreal>* x);

        void fillKernelsList(QVector<int>* kernelsIDs);
            void addProductKernelsToTheList(QVector<int>* kernelsIDs);
            void addRadialKernelsToTheList(QVector<int>* kernelsIDs);



};

#endif // KERNELDENSITYESTIMATOR_H
