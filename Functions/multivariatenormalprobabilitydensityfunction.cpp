#include "multivariatenormalprobabilitydensityfunction.h"

multivariateNormalProbabilityDensityFunction::multivariateNormalProbabilityDensityFunction(QVector<qreal> *means, QVector<qreal> *stDevs)
    : means(means)
{
    qreal correlationCoefficient = 0.5;
    QVector<QVector<qreal>*> covarianceMatrix;

    fillCovarianceMatrix(correlationCoefficient, stDevs, &covarianceMatrix);

    covarianceMatrixDeterminant = countMatrixDeterminantRecursively(&covarianceMatrix);

    fillInverseMatrix(&covarianceMatrix, &inverseCovarianceMatrix);
}

qreal multivariateNormalProbabilityDensityFunction::getValue(point *arguments)
{
    return 0;
}
