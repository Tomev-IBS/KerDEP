#ifndef MATRIXOPERATIONSLIBRARY_H
#define MATRIXOPERATIONSLIBRARY_H

// TODO: Write Matrix Class

#include <QObject>
#include <QVector>

typedef QVector<QVector<qreal>*>* matrixPtr;
typedef QVector<QVector<qreal>*> matrix;

void fillCovarianceMatrix(qreal correlationCoefficient, QVector<qreal>* stDevs, matrixPtr covarianceMatrix);

qreal countMatrixDeterminantRecursively(matrixPtr baseMatrix);

void fillCholeskyDecompositionMatrix(matrixPtr baseMatrix, matrixPtr decomposedMatrix);

void fillTransposedMatrix(matrixPtr baseMatrix, matrixPtr transposedMatrix);

void fillCofactorMatrix(matrixPtr baseMatrix, matrixPtr cofactorMatrix);

void fillInverseMatrix(matrixPtr baseMatrix, matrixPtr inverseMatrix);

void fillCopiedMatrix(matrixPtr baseMatrix, matrixPtr copy);


#endif // MATRIXOPERATIONSLIBRARY_H
