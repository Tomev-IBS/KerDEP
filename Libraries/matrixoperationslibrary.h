#ifndef MATRIXOPERATIONSLIBRARY_H
#define MATRIXOPERATIONSLIBRARY_H

// TODO: Write Matrix Class

#include <QObject>
#include <QVector>

typedef QVector<QVector<qreal>*>* matrixPtr;
typedef QVector<QVector<qreal>*> matrix;

void fillCovarianceMatrix(qreal correlationCoefficient, QVector<qreal>* stDevs, matrixPtr covarianceMatrix);

qreal countMatrixDeterminantRecursively(matrixPtr matrix);

void fillCholeskyDecompositionMatrix(matrixPtr matrix, matrixPtr decomposedMatrix);

void fillTransposedMatrix(matrixPtr matrix, matrixPtr transposedMatrix);

void fillInverseMatrix(matrixPtr matrix, matrixPtr inverseMatrix);


#endif // MATRIXOPERATIONSLIBRARY_H
