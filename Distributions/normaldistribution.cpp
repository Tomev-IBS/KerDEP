#include "normaldistribution.h"

#include "QDebug"
#include "QtMath"

normalDistribution::normalDistribution(int seed)
{
    generator = std::default_random_engine(seed);
    distributions.append(new std::normal_distribution<qreal>(0.0, 1.0));
}

normalDistribution::normalDistribution(int seed, QVector<qreal> *mean, QVector<qreal> *stDev)
{
    generator = std::default_random_engine(seed);

    for(int distributionNumber = 0; distributionNumber < mean->size(); ++distributionNumber)
        distributions.append(new std::normal_distribution<qreal>(mean->at(distributionNumber),
                                                                 stDev->at(distributionNumber)));


    QVector<QVector<qreal>*> covarianceMatrix;

    fillCovarianceMatrix(&covarianceMatrix, stDev);

    fillA(&covarianceMatrix);

    foreach (QVector<qreal>* a, A)
    {
        qDebug() << *a;
    }

}

void normalDistribution::getValue(QVector<qreal> *result)
{


    result->append((*distributions.last())(generator));
}

void normalDistribution::fillCovarianceMatrix(QVector<QVector<qreal>* > *covarianceMatrix, QVector<qreal> *stDev)
{
    for(int i = 0; i < stDev->size(); ++i)
    {

        covarianceMatrix->append(new QVector<qreal>());

        for(int j = 0; j < stDev->size(); ++j)
        {
            if(i==j)
            {
                covarianceMatrix->at(i)->append(stDev->at(i) * stDev->at(i));
            }
            else
            {
                covarianceMatrix->at(i)->append(stDev->at(i) * stDev->at(j) * correlationCoefficient);
            }
        }
    }
}

void normalDistribution::fillA(QVector<QVector<qreal>* > *covarianceMatrix)
{
    // Decomposing covariance matrix using Cholesky decomposition

    qreal value;
    A.clear();

    for(int rowNum = 0; rowNum < covarianceMatrix->size(); ++rowNum)
    {
        A.append(new QVector<qreal>());

        for(int columnNum = 0; columnNum < covarianceMatrix->size(); ++columnNum)
        {
            if(columnNum > rowNum)
            {
                value = 0;
            }
            else if(rowNum == columnNum)
            {
                value = covarianceMatrix->at(columnNum)->at(rowNum);

                for(int k = 0; k < rowNum -1; ++k)
                    value -= qPow(A.at(k)->at(rowNum), 2);

                value = qSqrt(value);
            }
            else
            {
                value = covarianceMatrix->at(rowNum)->at(columnNum);

                for(int k = 0; k < rowNum; ++k)
                    value -= A.at(k)->at(rowNum) * A.at(k)->at(columnNum);

                value /= A.at(columnNum)->at(columnNum);
            }

            A.at(rowNum)->append(value);
        }
    }
}
