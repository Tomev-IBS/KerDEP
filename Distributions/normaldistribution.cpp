#include "normaldistribution.h"
#include "../Libraries/matrixoperationslibrary.h"

#include "QDebug"
#include "QtMath"

normalDistribution::normalDistribution(int seed, QVector<qreal> *means, QVector<qreal> *stDevs)
    :   means(means)
{
    generator = std::default_random_engine(seed);

    for(int distributionNumber = 0; distributionNumber < means->size(); ++distributionNumber)
        distributions.append(new std::normal_distribution<qreal>(means->at(distributionNumber),
                                                                 stDevs->at(distributionNumber)));

    qreal correlationCoefficient = 0.5;
    QVector<QVector<qreal>*> covarianceMatrix;

    fillCovarianceMatrix(correlationCoefficient, stDevs, &covarianceMatrix);

    fillCholeskyDecompositionMatrix(&covarianceMatrix, &A);
}

void normalDistribution::getValue(QVector<qreal> *result)
{
    // Generate vector Z of n values from random distribution

    QVector<qreal> Z;

    for(int i = 0; i < A.size(); ++i)
        Z.append((*distributions.last())(generator));

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

void normalDistribution::fillA(QVector<QVector<qreal>* > *covarianceMatrix)
{

}
