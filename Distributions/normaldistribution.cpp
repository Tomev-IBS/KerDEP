#include "normaldistribution.h"
#include "../Libraries/matrixoperationslibrary.h"

#include <QDebug>

normalDistribution::normalDistribution(int seed, vector<double> *means,
                                       vector<double> *stDevs, double maxMean) :
    means(means), stDevs(stDevs), _maxMean(maxMean)
{
    generator = std::default_random_engine(seed);

    //this->means = vector<double>(*means);
    //this->stDevs = vector<double>(*stDevs);

    double correlationCoefficient = 0.5;
    matrix covarianceMatrix;

    fillCovarianceMatrix(correlationCoefficient, stDevs, &covarianceMatrix);

    fillCholeskyDecompositionMatrix(&covarianceMatrix, &A);
}

void normalDistribution::getValue(vector<double> *result)
{
    // Generate vector Z of n values from random distribution

    vector<double> Z;
    std::normal_distribution<double> normalDis(0,1);

    for(size_t i = 0; i < A.size(); ++i)
        Z.push_back(normalDis(generator));

    // Generete result according to X = u + AZ

    double value;

    for(size_t i = 0; i < A.size(); ++i)
    {
        value = 0.0;

        for(size_t j = 0; j < A.size(); ++j)
            value += A.at(i)->at(j) * Z.at(j);

        value += means->at(i);

        result->push_back(value);
    }
}

void normalDistribution::increaseMeans(double addend, int index)
{
  // Update mean at index, if it has been provided.
  if(index > -1 || means->size() > index){
    if(means->at(index) < _maxMean){
      (*means)[index] += addend;
    }
    return;
  }

  // Otherwise update all means
  for(size_t i = 0; i < means->size(); ++i)
  {
    // TODO: FIXED THRESHOLD FOR RESEARCHES
    if(means->at(i) < _maxMean)
    {
      (*means)[i] += addend;
    }
  }
}
