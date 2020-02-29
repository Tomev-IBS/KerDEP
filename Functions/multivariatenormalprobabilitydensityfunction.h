#ifndef MULTIVARIATENORMALPROBABILITYDENSITYFUNCTION_H
#define MULTIVARIATENORMALPROBABILITYDENSITYFUNCTION_H

#include "function.h"
#include "../Libraries/matrixoperationslibrary.h"

class multivariateNormalProbabilityDensityFunction : public function
{
    public:
        multivariateNormalProbabilityDensityFunction(vector<double>* means, vector<double>* stDevs);

        double getValue(point* arguments);

    private:

        matrix inverseCovarianceMatrix;
        vector<double> means;
        double covarianceMatrixDeterminant;
};

#endif // MULTIVARIATENORMALPROBABILITYDENSITYFUNCTION_H
