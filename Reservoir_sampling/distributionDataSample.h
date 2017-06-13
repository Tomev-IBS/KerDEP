#ifndef DISTRIBUTIONDATASAMPLE_H
#define DISTRIBUTIONDATASAMPLE_H

#include <QVector>

#include "sample.h"

class distributionDataSample : public sample
{
    public:
        QVector<qreal> values;
        void print();
};

#endif // DISTRIBUTIONDATASAMPLE_H
