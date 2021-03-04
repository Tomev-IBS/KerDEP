#ifndef COMPLEXFUNCTION_H
#define COMPLEXFUNCTION_H

#include "function.h"

#include <memory>

class complexFunction : public function
{
    public:
        complexFunction(vector<double> *contributions, vector<std::shared_ptr<function>> *elementalFunctions);
        ~complexFunction();

        double getValue(point* args);

    private:
        vector<double> contributions;
        vector<std::shared_ptr<function>> elementalFunctions;

};

#endif // COMPLEXFUNCTION_H
