#ifndef FUNCTION_H
#define FUNCTION_H

#include <vector>

using std::vector;
typedef vector<double> point;

class function
{
    public:
        virtual double getValue(point* arguments) = 0;
        virtual ~function(){};
};

#endif // FUNCTION_H
