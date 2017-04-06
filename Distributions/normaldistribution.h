#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H

#include "distribution.h"

#include "random"
#include <QVector>

class normalDistribution : public distribution
{
    public:
        normalDistribution(int seed, QVector<qreal>* means, QVector<qreal>* stDevs);

        void getValue(QVector<qreal>* result);

    private:
        QVector<qreal>* means;
        QVector<std::normal_distribution<qreal>*> distributions;
        qreal correlationCoefficient = 0.5;

        QVector<QVector<qreal>*> A;

        void fillCovarianceMatrix(QVector<QVector<qreal> *> *covarianceMatrix, QVector<qreal>* stDev);
        void fillA(QVector<QVector<qreal> *> *covarianceMatrix);

};

#endif // NORMALDISTRIBUTION_H
