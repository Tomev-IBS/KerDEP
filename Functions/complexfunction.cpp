#include "complexfunction.h"

complexFunction::complexFunction(vector<double> *contributions, vector<std::shared_ptr<function>> *elementalFunctions)
{
    this->contributions         = vector<double>(*contributions);
    this->elementalFunctions    = vector<std::shared_ptr<function>>(*elementalFunctions);
}

complexFunction::~complexFunction()
{
    contributions.clear();
    elementalFunctions.clear();
}

double complexFunction::getValue(point *args)
{
    double result = 0;

    for(size_t functionIndex = 0; functionIndex < elementalFunctions.size(); ++functionIndex)
        result += elementalFunctions.at(functionIndex).get()->getValue(args) * contributions.at(functionIndex)/100.0;

    return result;
}
