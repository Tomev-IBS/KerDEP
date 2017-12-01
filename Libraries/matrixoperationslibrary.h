#ifndef MATRIXOPERATIONSLIBRARY_H
#define MATRIXOPERATIONSLIBRARY_H

// TODO: Write Matrix Class

#include <QObject>
#include <QVector>
#include <memory>

typedef QVector<std::shared_ptr<QVector<qreal>>>* matrixPtr;
typedef QVector<std::shared_ptr<QVector<qreal>>> matrix;

void fillCovarianceMatrix(qreal correlationCoefficient, QVector<qreal>* stDevs, matrix *covarianceMatrix);

qreal countMatrixDeterminantRecursively(matrixPtr baseMatrix);

void fillCholeskyDecompositionMatrix(matrixPtr baseMatrix, matrixPtr decomposedMatrix);

void fillTransposedMatrix(matrixPtr baseMatrix, matrixPtr transposedMatrix);

void fillCofactorMatrix(matrixPtr baseMatrix, matrixPtr cofactorMatrix);

void fillInverseMatrix(matrixPtr baseMatrix, matrixPtr inverseMatrix);

void fillCopiedMatrix(matrixPtr baseMatrix, matrixPtr copy);

void removeMatrixColumn(matrixPtr baseMatrix, int columnIndex);

void removeMatrixRow(matrixPtr baseMatrix, int rowIndex);


#endif // MATRIXOPERATIONSLIBRARY_H
