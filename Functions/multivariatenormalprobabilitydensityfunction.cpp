#include "multivariatenormalprobabilitydensityfunction.h"

#include "QtMath"
#include <QDebug>

multivariateNormalProbabilityDensityFunction::multivariateNormalProbabilityDensityFunction(QVector<qreal> *means, QVector<qreal> *stDevs)
{
    qreal correlationCoefficient = 0.5;
    matrix covarianceMatrix;

    this->means = QVector<qreal>(*means);

    fillCovarianceMatrix(correlationCoefficient, stDevs, &covarianceMatrix);

    covarianceMatrixDeterminant = countMatrixDeterminantRecursively(&covarianceMatrix);

    fillInverseMatrix(&covarianceMatrix, &inverseCovarianceMatrix);
}

qreal multivariateNormalProbabilityDensityFunction::getValue(point *arguments)
{
    if(covarianceMatrixDeterminant == 0)
    {
        qDebug() << "Determinant is 0.";
        return -1;
    }

    QVector<qreal> vectorMatrixProduct;
    qreal value, result = 0;

    for(int rowIndex = 0; rowIndex < inverseCovarianceMatrix.size(); ++rowIndex)
    {
        value = 0;

        for(int columnIndex = 0; columnIndex < inverseCovarianceMatrix.at(rowIndex)->size(); ++columnIndex)
        {
            value += inverseCovarianceMatrix.at(rowIndex)->at(columnIndex)
                    * (arguments->at(columnIndex) - means.at(columnIndex));
        }

        vectorMatrixProduct.append(value);
    }

    for(int i = 0; i < vectorMatrixProduct.size(); ++i)
        result += vectorMatrixProduct.at(i) * (arguments->at(i) - means.at(i));

    result /= -2;

    result = exp(result);
    result /= qSqrt(qPow(2*M_PI, arguments->size()) * covarianceMatrixDeterminant);

    return result;
}
