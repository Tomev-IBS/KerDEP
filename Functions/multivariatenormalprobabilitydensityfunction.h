#ifndef MULTIVARIATENORMALPROBABILITYDENSITYFUNCTION_H
#define MULTIVARIATENORMALPROBABILITYDENSITYFUNCTION_H

#include <QObject>
#include "function.h"
#include "../Libraries/matrixoperationslibrary.h"

class multivariateNormalProbabilityDensityFunction : public function
{
    public:
        multivariateNormalProbabilityDensityFunction(QVector<qreal>* means, QVector<qreal>* stDevs);

        qreal getValue(point* arguments);

    private:

        matrix inverseCovarianceMatrix;
        QVector<qreal> means;
        qreal covarianceMatrixDeterminant;
};

#endif // MULTIVARIATENORMALPROBABILITYDENSITYFUNCTION_H
