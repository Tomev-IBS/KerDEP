#ifndef COMPLEXDISTRIBUTION_H
#define COMPLEXDISTRIBUTION_H

#include "distribution.h"
#include <random>
#include <QVector>
#include <QObject>

class complexDistribution : public distribution
{
    public:

        complexDistribution(int seed, QVector<distribution*>* elementalDistributions, QVector<qreal>* contributions);

        void getValue(QVector<qreal>* result);
        void increaseMeans(qreal addend);

    private:
        QVector<distribution*>* elementalDistributions;
        QVector<qreal>* contributions;

        const std::uniform_real_distribution<qreal> uniformDistribution;

        int randomizeDistributionIndex();
};

#endif // COMPLEXDISTRIBUTION_H
