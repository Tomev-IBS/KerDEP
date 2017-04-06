#include "matrixoperationslibrary.h"
#include "QtMath"



void fillCovarianceMatrix(qreal correlationCoefficient, QVector<qreal> *stDevs, QVector<QVector<qreal> *> *covarianceMatrix)
{
    for(int i = 0; i < stDevs->size(); ++i)
    {
        covarianceMatrix->append(new QVector<qreal>());

        for(int j = 0; j < stDevs->size(); ++j)
        {
            if(i==j)
                covarianceMatrix->at(i)->append(stDevs->at(i) * stDevs->at(i));
            else
                covarianceMatrix->at(i)->append(stDevs->at(i) * stDevs->at(j) * correlationCoefficient);
        }
    }
}

void fillCholeskyDecompositionMatrix(matrixPtr matrix, matrixPtr decomposedMatrix)
{
    // Decomposing covariance matrix using Cholesky decomposition
    // Note that it was used only for normal distributions covariance matrices.
    // Warning: May note work correctly.

    qreal value;
    decomposedMatrix->clear();

    for(int rowNum = 0; rowNum < matrix->size(); ++rowNum)
    {
        decomposedMatrix->append(new QVector<qreal>());

        for(int columnNum = 0; columnNum < matrix->size(); ++columnNum)
        {
            if(columnNum > rowNum)
            {
                value = 0;
            }
            else if(rowNum == columnNum)
            {
                value = matrix->at(columnNum)->at(rowNum);

                for(int k = 0; k < rowNum -1; ++k)
                    value -= qPow(decomposedMatrix->at(k)->at(rowNum), 2);

                value = qSqrt(value);
            }
            else
            {
                value = matrix->at(rowNum)->at(columnNum);

                for(int k = 0; k < rowNum; ++k)
                    value -= decomposedMatrix->at(k)->at(rowNum) * decomposedMatrix->at(k)->at(columnNum);

                value /= decomposedMatrix->at(columnNum)->at(columnNum);
            }

            decomposedMatrix->at(rowNum)->append(value);
        }
    }
}

qreal countMatrixDeterminantRecursively(matrixPtr matrix)
{
    return 0;
}

void fillInverseMatrix(matrixPtr matrix, matrixPtr inverseMatrix)
{

}

void fillTransposedMatrix(matrixPtr matrix, matrixPtr transposedMatrix)
{

}
