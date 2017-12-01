#ifndef COMPLEXFUNCTION_H
#define COMPLEXFUNCTION_H

#include "function.h"

#include <memory>

class complexFunction : public function
{
    public:
        complexFunction(QVector<qreal>* contributions, QVector<std::shared_ptr<function>>* elementalFunctions);
        ~complexFunction();

        qreal getValue(point* args);

    private:
        QVector<qreal> contributions;
        QVector<std::shared_ptr<function>> elementalFunctions;

};

#endif // COMPLEXFUNCTION_H
