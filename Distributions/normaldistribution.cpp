#include "normaldistribution.h"
#include "../Libraries/matrixoperationslibrary.h"

#include <QDebug>

normalDistribution::normalDistribution(int seed, vector<double> *means,
                                       vector<double> *stDevs) :
    means(means), stDevs(stDevs)
{
    generator = std::default_random_engine(seed);

    //this->means = vector<double>(*means);
    //this->stDevs = vector<double>(*stDevs);

    double correlationCoefficient = 0;
    matrix covarianceMatrix;

    fillCovarianceMatrix(correlationCoefficient, stDevs, &covarianceMatrix);

    fillCholeskyDecompositionMatrix(&covarianceMatrix, &A);

    qDebug() << A.size() << " is the size of A.";

    for(size_t i = 0; i < A.size(); ++i){
      auto s = seed + i;
      generators.push_back(std::default_random_engine(s));
    }
}

void normalDistribution::getValue(vector<double> *result)
{
    // Generate vector Z of n values from random distribution

    vector<double> Z;

    for(size_t i = 0; i < A.size(); ++i) {
      std::normal_distribution<double> normalDis(0,1);
      Z.push_back(normalDis(generators[i]));
    }

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
    (*means)[index] += addend;
    return;
  }

  // Otherwise update all means
  for(size_t i = 0; i < means->size(); ++i)
  {
    (*means)[i] += addend;
  }
}

void normalDistribution::setMeans(double newMean, int index) {
  // Update mean at index, if it has been provided.
  if(index > -1 || means->size() > index){
    (*means)[index] = newMean;
    return;
  }

  // Otherwise update all means
  for(size_t i = 0; i < means->size(); ++i)
  {
    (*means)[i] = newMean;
  }
}
