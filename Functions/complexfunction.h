#ifndef COMPLEXFUNCTION_H
#define COMPLEXFUNCTION_H

#include "function.h"

class complexFunction : public function
{
    public:
        complexFunction(QVector<qreal>* contributions, QVector<function*>* elementalFunctions);

        qreal getValue(point* args);

    private:
        QVector<qreal>* contributions;
        QVector<function*>* elementalFunctions;

};

#endif // COMPLEXFUNCTION_H
