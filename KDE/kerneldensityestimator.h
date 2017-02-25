#ifndef KERNELDENSITYESTIMATOR_H
#define KERNELDENSITYESTIMATOR_H

#include "../Functions/gaussianprobabilitydensityfunction.h"
#include "../Distributions/distribution.h"

#include <QObject>

class kernelDensityEstimator
{
    public:
        kernelDensityEstimator(int sampleSize, qreal smoothingParameter, function* kernel, distribution* targetDistribution);

        qreal getValue(qreal x);

    private:

        QVector<qreal> samples;

        qreal smoothingParameter;

        function* kernel;
};

#endif // KERNELDENSITYESTIMATOR_H
