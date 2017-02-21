#ifndef GAUSSIANPROBABILITYDENSITYFUNCTION_H
#define GAUSSIANPROBABILITYDENSITYFUNCTION_H

#include "function.h"

class gaussianProbabilityDensityFunction : public function
{
    public:
        gaussianProbabilityDensityFunction();
        gaussianProbabilityDensityFunction(qreal mean, qreal standardDeviation);

        qreal getValue(QVector<qreal>* arguments); // Only one argument should be passed

    private:
        qreal   mean = 0.0,
                standardDeviation = 1.0;

};

#endif // GAUSSIANPROBABILITYDENSITYFUNCTION_H
