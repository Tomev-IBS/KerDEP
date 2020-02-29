#include "matrixoperationslibrary.h"
#include <cmath>

void fillCovarianceMatrix(double correlationCoefficient, vector<double> *stDevs, matrix *covarianceMatrix)
{
    for(size_t i = 0; i < stDevs->size(); ++i)
    {
        covarianceMatrix->push_back(std::make_shared<vector<double>>());

        for(size_t j = 0; j < stDevs->size(); ++j)
        {
            if(i==j)
                covarianceMatrix->at(i)->push_back(stDevs->at(i) * stDevs->at(i));
            else
                covarianceMatrix->at(i)->push_back(stDevs->at(i) * stDevs->at(j) * correlationCoefficient);
        }
    }
}

void fillCholeskyDecompositionMatrix(matrixPtr baseMatrix, matrixPtr decomposedMatrix)
{
    // Decomposing covariance matrix using Cholesky decomposition
    // Note that it was used only for normal distributions covariance matrices.
    // Warning: May note work correctly.

    double value;
    decomposedMatrix->clear();

    for(size_t rowNum = 0; rowNum < baseMatrix->size(); ++rowNum)
    {
        decomposedMatrix->push_back(std::make_shared<vector<double>>());

        for(size_t columnNum = 0; columnNum < baseMatrix->size(); ++columnNum)
        {
            if(columnNum > rowNum)
            {
                value = 0;
            }
            else if(rowNum == columnNum)
            {
                value = baseMatrix->at(columnNum)->at(rowNum);

                for(size_t k = 0; k < rowNum -1 && rowNum != 0; ++k)
                    value -= pow(decomposedMatrix->at(k)->at(rowNum), 2);

                value = sqrt(value);
            }
            else
            {
                value = baseMatrix->at(rowNum)->at(columnNum);

                for(size_t k = 0; k < rowNum; ++k)
                    value -= decomposedMatrix->at(k)->at(rowNum) * decomposedMatrix->at(k)->at(columnNum);

                value /= decomposedMatrix->at(columnNum)->at(columnNum);
            }

            decomposedMatrix->at(rowNum)->push_back(value);
        }
    }
}

double countMatrixDeterminantRecursively(matrixPtr baseMatrix)
{
    // Recursive definition from https://pl.wikipedia.org/wiki/Wyznacznik
    // Column index is always = 0

    // Warning: not checking if matrix is N x N

    if(baseMatrix->size() == 1)
        return baseMatrix->at(0)->at(0);

    matrix copiedMatrixWithoutJthColumn, matrixWithoutJthColumn;
    double determinant = 0, addend;
    int columnIndex = 0;

    fillCopiedMatrix(baseMatrix, &matrixWithoutJthColumn);

    // Remove first value of each row (first column)
    removeMatrixColumn(&matrixWithoutJthColumn, columnIndex);

    for(size_t rowIndex = 0; rowIndex < baseMatrix->size(); ++rowIndex)
    {
        fillCopiedMatrix(&matrixWithoutJthColumn, &copiedMatrixWithoutJthColumn);

        // Remove i-th row
        removeMatrixRow(&copiedMatrixWithoutJthColumn, rowIndex);

        addend = pow(-1.0, rowIndex + columnIndex);
        addend *= baseMatrix->at(rowIndex)->at(columnIndex);
        addend *= countMatrixDeterminantRecursively(&copiedMatrixWithoutJthColumn);

        determinant += addend;
    }

    return determinant;
}

void fillInverseMatrix(matrixPtr baseMatrix, matrixPtr inverseMatrix)
{
    matrix cofactorMatrix, transposedCofactorMatrix;
    double determinant;

    fillCofactorMatrix(baseMatrix, &cofactorMatrix);
    fillTransposedMatrix(&cofactorMatrix, &transposedCofactorMatrix);

    determinant = countMatrixDeterminantRecursively(baseMatrix);

    if(determinant == 0)
      return;

    inverseMatrix->clear();

    for(auto row : transposedCofactorMatrix)
    {
      inverseMatrix->push_back(std::make_shared<vector<double>>());

      for(double value : *row)
        inverseMatrix->back()->push_back(value/determinant);
    }
}

void fillTransposedMatrix(matrixPtr baseMatrix, matrixPtr transposedMatrix)
{
    transposedMatrix->clear();

    // Insert row for each column
    for(size_t i = 0; i < baseMatrix->at(0)->size(); ++i)
        transposedMatrix->push_back(std::make_shared<vector<double>>());

    for(size_t rowIndex = 0; rowIndex < baseMatrix->size(); ++rowIndex)
    {
        for(size_t columnIndex = 0; columnIndex < baseMatrix->at(rowIndex)->size(); ++columnIndex)
        {
            transposedMatrix->at(columnIndex)->push_back(baseMatrix->at(rowIndex)->at(columnIndex));
        }
    }
}

void fillCofactorMatrix(matrixPtr baseMatrix, matrixPtr cofactorMatrix)
{
    // According to https://pl.wikipedia.org/wiki/Macierz_do%C5%82%C4%85czona

    if(baseMatrix->size() != baseMatrix->at(0)->size())
        return;

    if(baseMatrix->size() == 1)
    {
        fillCopiedMatrix(baseMatrix, cofactorMatrix);
        return;
    }

    cofactorMatrix->clear();

    matrix submatrix;

    for(size_t rowIndex = 0; rowIndex < baseMatrix->size(); ++rowIndex)
    {
        cofactorMatrix->push_back(std::make_shared<vector<double>>());

        for(size_t columnIndex = 0; columnIndex < baseMatrix->at(rowIndex)->size(); ++columnIndex)
        {
            fillCopiedMatrix(baseMatrix, &submatrix);
            removeMatrixRow(&submatrix, rowIndex);
            removeMatrixColumn(&submatrix, columnIndex);

            cofactorMatrix->back()->push_back(pow(-1, rowIndex + columnIndex) * countMatrixDeterminantRecursively(&submatrix));
        }
    }
}

void fillCopiedMatrix(matrixPtr baseMatrix, matrixPtr copy)
{
    copy->clear();

    for(std::shared_ptr<vector<double>> row : *baseMatrix)
    {
        copy->push_back(row);
    }
}

void removeMatrixColumn(matrixPtr baseMatrix, int columnIndex)
{
    // Remove first value of each row (first column)
    for(std::shared_ptr<vector<double>> row : *baseMatrix)
        row->erase(row->begin() + columnIndex);
}

void removeMatrixRow(matrixPtr baseMatrix, int rowIndex)
{
    baseMatrix->erase(baseMatrix->begin() + rowIndex);
}
