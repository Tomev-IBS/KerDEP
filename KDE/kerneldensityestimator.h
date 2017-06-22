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
        kernelDensityEstimator(QVector<QVector<qreal>*>* samples, QVector<qreal>* smoothingParameter, QVector<QString>* carriersRestrictions, int kernelType, QVector<int>* kernelsIDs);

        void setSamples(QVector<QVector<qreal>*>* samples);

        qreal getValue(QVector<qreal>* x);

    private:

        int kernelType;

        QVector<QVector<qreal>*>    samples;
        QVector<kernelPtr>          kernels;
        QVector<qreal>              smoothingParameters;
        QVector<QString>            carriersRestrictions;

        qreal getProductKernelValue(QVector<qreal>* x);
        qreal getRadialKernelValue(QVector<qreal>* x);

        void fillKernelsList(QVector<int>* kernelsIDs);
            void addProductKernelsToTheList(QVector<int>* kernelsIDs);
            void addRadialKernelsToTheList(QVector<int>* kernelsIDs);

        int partitionCharacteristicFunction(qreal carrier, qreal restriction);
};

#endif // KERNELDENSITYESTIMATOR_H
