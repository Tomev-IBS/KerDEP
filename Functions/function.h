#ifndef FUNCTION_H
#define FUNCTION_H

#include <QObject>
#include <QVector>

class function
{
    public:
        virtual qreal getValue(QVector<qreal>* arguments) = 0;
};

#endif // FUNCTION_H
