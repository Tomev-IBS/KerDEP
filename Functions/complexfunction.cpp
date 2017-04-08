#include "complexfunction.h"

complexFunction::complexFunction(QVector<qreal> *contributions, QVector<function *> *elementalFunctions) :
    contributions(contributions), elementalFunctions(elementalFunctions)
{}

qreal complexFunction::getValue(point *args)
{
    qreal result = 0;

    for(int functionIndex = 0; functionIndex < elementalFunctions->size(); ++functionIndex)
        result += elementalFunctions->at(functionIndex)->getValue(args) * contributions->at(functionIndex)/100;

    return result;
}
