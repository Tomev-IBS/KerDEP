#include "multivariatenormalprobabilitydensityfunction.h"

#include "QtMath"
#include <QDebug>

multivariateNormalProbabilityDensityFunction::multivariateNormalProbabilityDensityFunction(vector<double> *means, vector<double> *stDevs, int covarianceCoefficient)
{
    double correlationCoefficient = covarianceCoefficient;
    matrix covarianceMatrix;

    this->means = vector<double>(*means);

    fillCovarianceMatrix(correlationCoefficient, stDevs, &covarianceMatrix);

    covarianceMatrixDeterminant = countMatrixDeterminantRecursively(&covarianceMatrix);

    fillInverseMatrix(&covarianceMatrix, &inverseCovarianceMatrix);
}

double multivariateNormalProbabilityDensityFunction::getValue(point *arguments)
{
    if(covarianceMatrixDeterminant == 0)
    {
        qDebug() << "Determinant is 0.";
        return -1;
    }

    vector<double> vectorMatrixProduct;
    double value, result = 0;

    for(size_t rowIndex = 0; rowIndex < inverseCovarianceMatrix.size(); ++rowIndex)
    {
        value = 0;

        for(size_t columnIndex = 0; columnIndex < inverseCovarianceMatrix.at(rowIndex)->size(); ++columnIndex)
        {
            value += inverseCovarianceMatrix.at(rowIndex)->at(columnIndex)
                    * (arguments->at(columnIndex) - means.at(columnIndex));
        }

        vectorMatrixProduct.push_back(value);
    }

    for(size_t i = 0; i < vectorMatrixProduct.size(); ++i)
        result += vectorMatrixProduct.at(i) * (arguments->at(i) - means.at(i));

    result /= -2;

    result = exp(result);
    result /= qPow(2 * M_PI, arguments->size() / 2.0);
    result /= qSqrt(covarianceMatrixDeterminant);

    return result;
}

void multivariateNormalProbabilityDensityFunction::setMeans(const vector<double> &newMeans)
{
  means = newMeans;
}
