#include "kerneldensityestimator.h"

#include "QDebug"

kernelDensityEstimator::kernelDensityEstimator(int sampleSize, qreal smoothingParameter, function* kernel, distribution* targetDistribution)
    : kernel(kernel), smoothingParameter(smoothingParameter)
{
    for(int sampleNumber = 0; sampleNumber < sampleSize; ++sampleNumber)
        samples.append(targetDistribution->getValue());
}

qreal kernelDensityEstimator::getValue(qreal x)
{
    qreal result = 0.0;

    QVector<qreal>* tempValueHolder = new QVector<qreal>();

    foreach(qreal sample, samples)
    {
        tempValueHolder->clear();
        tempValueHolder->append((x-sample)/smoothingParameter);
        result += kernel->getValue(tempValueHolder);
    }

    result /= samples.size() * smoothingParameter;

    return result;
}
