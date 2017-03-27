#ifndef FUNCTION_H
#define FUNCTION_H

#include <QObject>
#include <QVector>

typedef QVector<qreal> point;

class function
{
    public:
        virtual qreal getValue(point* arguments) = 0;
};

#endif // FUNCTION_H
