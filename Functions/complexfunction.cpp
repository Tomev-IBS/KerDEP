#include "complexfunction.h"

complexFunction::complexFunction(QVector<qreal> *contributions, QVector<function *> *elementalFunctions)
{
    this->contributions         = QVector<qreal>(*contributions);
    this->elementalFunctions    = QVector<function*>(*elementalFunctions);
}

qreal complexFunction::getValue(point *args)
{
    qreal result = 0;

    for(int functionIndex = 0; functionIndex < elementalFunctions.size(); ++functionIndex)
        result += elementalFunctions.at(functionIndex)->getValue(args) * contributions.at(functionIndex)/100;

    return result;
}
