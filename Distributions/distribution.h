#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <QObject>
#include <random>

class distribution
{
    public:
        virtual void getValue(QVector<qreal>* result) = 0;

        virtual void increaseMeans(qreal addend) = 0;

    protected:

        std::default_random_engine generator;
};

#endif // DISTRIBUTION_H
