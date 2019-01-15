#include "normaldistribution.h"
#include "../Libraries/matrixoperationslibrary.h"

#include <QDebug>

normalDistribution::normalDistribution(int seed, QVector<qreal> *means, QVector<qreal> *stDevs) :
    means(means), stDevs(stDevs)
{
    generator = std::default_random_engine(seed);

    //this->means = QVector<qreal>(*means);
    //this->stDevs = QVector<qreal>(*stDevs);

    qreal correlationCoefficient = 0.5;
    matrix covarianceMatrix;

    fillCovarianceMatrix(correlationCoefficient, stDevs, &covarianceMatrix);

    fillCholeskyDecompositionMatrix(&covarianceMatrix, &A);
}

void normalDistribution::getValue(QVector<qreal> *result)
{
    // Generate vector Z of n values from random distribution

    QVector<qreal> Z;
    std::normal_distribution<qreal> normalDis(0,1);

    for(int i = 0; i < A.size(); ++i)
        Z.append(normalDis(generator));

    // Generete result according to X = u + AZ

    qreal value;

    for(int i = 0; i < A.size(); ++i)
    {
        value = 0.0;

        for(int j = 0; j < A.size(); ++j)
            value += A.at(i)->at(j) * Z.at(j);

        value += means->at(i);

        result->append(value);
    }
}

void normalDistribution::increaseMeans(qreal addend)
{
    for(int i = 0; i < means->size(); ++i)
    {
        // TODO: FIXED THRESHOLD FOR RESEARCHES
        if(means->at(i) < 1 && means->at(i) < 2)
        {
            means->push_back(means->at(i) + addend);
            means->pop_front();
        }
    }
}
