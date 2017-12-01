#include "complexfunction.h"

#include "QDebug"

complexFunction::complexFunction(QVector<qreal> *contributions, QVector<std::shared_ptr<function>> *elementalFunctions)
{
    this->contributions         = QVector<qreal>(*contributions);
    this->elementalFunctions    = QVector<std::shared_ptr<function>>(*elementalFunctions);
}

complexFunction::~complexFunction()
{
    contributions.clear();
    elementalFunctions.clear();
}

qreal complexFunction::getValue(point *args)
{
    qreal result = 0;

    for(int functionIndex = 0; functionIndex < elementalFunctions.size(); ++functionIndex)
        result += elementalFunctions.at(functionIndex).get()->getValue(args) * contributions.at(functionIndex)/100;

    return result;
}
