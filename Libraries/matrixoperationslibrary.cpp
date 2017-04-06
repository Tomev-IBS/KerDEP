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

void fillCholeskyDecompositionMatrix(matrixPtr baseMatrix, matrixPtr decomposedMatrix)
{
    // Decomposing covariance matrix using Cholesky decomposition
    // Note that it was used only for normal distributions covariance matrices.
    // Warning: May note work correctly.

    qreal value;
    decomposedMatrix->clear();

    for(int rowNum = 0; rowNum < baseMatrix->size(); ++rowNum)
    {
        decomposedMatrix->append(new QVector<qreal>());

        for(int columnNum = 0; columnNum < baseMatrix->size(); ++columnNum)
        {
            if(columnNum > rowNum)
            {
                value = 0;
            }
            else if(rowNum == columnNum)
            {
                value = baseMatrix->at(columnNum)->at(rowNum);

                for(int k = 0; k < rowNum -1; ++k)
                    value -= qPow(decomposedMatrix->at(k)->at(rowNum), 2);

                value = qSqrt(value);
            }
            else
            {
                value = baseMatrix->at(rowNum)->at(columnNum);

                for(int k = 0; k < rowNum; ++k)
                    value -= decomposedMatrix->at(k)->at(rowNum) * decomposedMatrix->at(k)->at(columnNum);

                value /= decomposedMatrix->at(columnNum)->at(columnNum);
            }

            decomposedMatrix->at(rowNum)->append(value);
        }
    }
}

qreal countMatrixDeterminantRecursively(matrixPtr baseMatrix)
{
    // Recursive definition from https://pl.wikipedia.org/wiki/Wyznacznik
    // Column index is always = 0

    // Warning: not checking if matrix is N x N

    if(baseMatrix->size() == 1)
        return baseMatrix->at(0)->at(0);

    matrix copiedMatrixWithoutJthColumn, matrixWithoutJthColumn;
    qreal determinant = 0, addend;
    int columnIndex = 0;

    fillCopiedMatrix(baseMatrix, &matrixWithoutJthColumn);

    // Remove first value of each row (first column)
    removeMatrixColumn(&matrixWithoutJthColumn, columnIndex);

    for(int rowIndex = 0; rowIndex < baseMatrix->size(); ++rowIndex)
    {
        fillCopiedMatrix(&matrixWithoutJthColumn, &copiedMatrixWithoutJthColumn);

        // Remove i-th row
        removeMatrixRow(&copiedMatrixWithoutJthColumn, rowIndex);

        addend = qPow(-1.0, rowIndex + columnIndex);
        addend *= baseMatrix->at(rowIndex)->at(columnIndex);
        addend *= countMatrixDeterminantRecursively(&copiedMatrixWithoutJthColumn);

        determinant += addend;
    }

    return determinant;
}

void fillInverseMatrix(matrixPtr baseMatrix, matrixPtr inverseMatrix)
{

}

void fillTransposedMatrix(matrixPtr baseMatrix, matrixPtr transposedMatrix)
{
    transposedMatrix->clear();

    // Insert row for each column
    foreach(auto column, *(baseMatrix->at(0)))
        transposedMatrix->append(new QVector<qreal>());

    for(int rowIndex = 0; rowIndex < baseMatrix->size(); ++rowIndex)
    {
        for(int columnIndex = 0; columnIndex < baseMatrix->at(rowIndex)->size(); ++columnIndex)
        {
            transposedMatrix->at(columnIndex)->append(baseMatrix->at(rowIndex)->at(columnIndex));
        }
    }
}

void fillCofactorMatrix(matrixPtr baseMatrix, matrixPtr cofactorMatrix)
{

}

void fillCopiedMatrix(matrixPtr baseMatrix, matrixPtr copy)
{
    copy->clear();

    foreach(QVector<qreal>* row, *baseMatrix)
    {
        copy->append(new QVector<qreal>());

        foreach (qreal position, *row)
            copy->last()->append(position);
    }
}

void removeMatrixColumn(matrixPtr baseMatrix, int columnIndex)
{
    // Remove first value of each row (first column)
    foreach (QVector<qreal>* row, *baseMatrix)
        row->removeAt(columnIndex);
}

void removeMatrixRow(matrixPtr baseMatrix, int rowIndex)
{
    baseMatrix->removeAt(rowIndex);
}
