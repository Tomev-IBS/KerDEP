#ifndef MATRIXOPERATIONSLIBRARY_H
#define MATRIXOPERATIONSLIBRARY_H

// TODO: Write Matrix Class
#include <memory>
#include <vector>

using std::vector;

typedef vector<std::shared_ptr<vector<double>>>* matrixPtr;
typedef vector<std::shared_ptr<vector<double>>> matrix;

void fillCovarianceMatrix(double correlationCoefficient, vector<double>* stDevs, matrix *covarianceMatrix);

double countMatrixDeterminantRecursively(matrixPtr baseMatrix);

void fillCholeskyDecompositionMatrix(matrixPtr baseMatrix, matrixPtr decomposedMatrix);

void fillTransposedMatrix(matrixPtr baseMatrix, matrixPtr transposedMatrix);

void fillCofactorMatrix(matrixPtr baseMatrix, matrixPtr cofactorMatrix);

void fillInverseMatrix(matrixPtr baseMatrix, matrixPtr inverseMatrix);

void fillCopiedMatrix(matrixPtr baseMatrix, matrixPtr copy);

void removeMatrixColumn(matrixPtr baseMatrix, int columnIndex);

void removeMatrixRow(matrixPtr baseMatrix, int rowIndex);


#endif // MATRIXOPERATIONSLIBRARY_H
