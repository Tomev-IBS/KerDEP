#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H

#include "distribution.h"
#include "../Libraries/matrixoperationslibrary.h"

#include "random"
#include <QVector>

class normalDistribution : public distribution
{
    public:
        normalDistribution(int seed, QVector<qreal>* means, QVector<qreal>* stDevs, double maxMean);

        void getValue(QVector<qreal>* result);
        void increaseMeans(qreal addend);

    private:      
        QVector<qreal>* means;
        QVector<qreal>* stDevs;
        QVector<std::normal_distribution<qreal>*> distributions;

        matrix A;

        double _maxMean = 0;

        void fillA(QVector<QVector<qreal> *> *covarianceMatrix);
};

#endif // NORMALDISTRIBUTION_H
