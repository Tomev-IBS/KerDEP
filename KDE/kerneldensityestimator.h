#ifndef KERNELDENSITYESTIMATOR_H
#define KERNELDENSITYESTIMATOR_H

#include "../Functions/gaussianprobabilitydensityfunction.h"

#include <QObject>

class kernelDensityEstimator
{
    public:
        kernelDensityEstimator(int seed, int sampleSize, qreal mean, qreal standardDeviation, function* kernel, qreal smoothingParameter);

        qreal getValue(qreal x);

    private:

        QVector<qreal> samples;

        qreal smoothingParameter;

        function* kernel;
};

#endif // KERNELDENSITYESTIMATOR_H
