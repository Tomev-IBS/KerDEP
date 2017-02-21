#include "kerneldensityestimator.h"

kernelDensityEstimator::kernelDensityEstimator(int seed, int sampleSize, qreal mean, qreal standardDeviation, function* kernel, qreal smoothingParameter)
    : kernel(kernel), smoothingParameter(smoothingParameter)
{
    std::default_random_engine generator(seed);

    std::normal_distribution<qreal> distribution(mean, standardDeviation);

    for(int sampleNumber = 0; sampleNumber < sampleSize; ++sampleNumber)
        samples.append(distribution(generator));
}

qreal kernelDensityEstimator::getValue(qreal x)
{
    qreal result = 0.0;

    foreach(qreal sample, samples)
        result += kernel->getValue(new QVector<qreal>{(x-sample)/smoothingParameter});

    result /= samples.size() * smoothingParameter;

    return result;
}
