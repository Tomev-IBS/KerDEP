#include "pluginsmoothingparametercounter.h"

#include <QDebug>

pluginSmoothingParameterCounter::pluginSmoothingParameterCounter(QVector<qreal>* samples)
    : samples(samples){}

qreal pluginSmoothingParameterCounter::count2ndRankPluginSmoothingParameter()
{
    // Page 82 of Kernel Estimators in System Analysis

    qreal h2 = -2.0;
    h2 *= countK6thDerivativeInPoint(0);
    h2 /= U;
    h2 /= samples->size();
    h2 /= countSmallC(8);
    h2 = qPow(h2, 1.0/9.0);

    return  countPluginSmoothingParameter(h2);
}

qreal pluginSmoothingParameterCounter::count3rdRankPluginSmoothingParameter()
{
    // Page 82 of Kernel Estimators in System Analysis

    qreal h3 = -2.0;
    h3 *= countK8thDerivativeInPoint(0);
    h3 /= U;
    h3 /= countSmallC(10);
    h3 /= samples->size();

    qDebug() << h3;

    h3 = qPow(qAbs(h3), 1.0/11.0);
    h3 = -h3;

    qDebug() << h3;

    qreal h2 = -2.0;
    h2 *= countK6thDerivativeInPoint(0);
    h2 /= U;
    h2 /= countCapitalC(8.0, h3);
    h2 /= samples->size();

    qDebug() << h2;

    h2 = qPow(h2, 1.0/9.0);

    return countPluginSmoothingParameter(h2);
}

qreal pluginSmoothingParameterCounter::countCapitalC(int xsi, qreal smoothingParameter)
{
    // Page 81 of Kernel Estimators in System Analysis, P. Kulczycki

    qreal (pluginSmoothingParameterCounter::* xsithKDerivative)(qreal);

    switch(xsi)
    {
        case 8:
            xsithKDerivative = &countK8thDerivativeInPoint;
        break;
        case 6:
            xsithKDerivative = &countK6thDerivativeInPoint;
        break;
        case 4:
        default:
            xsithKDerivative = &countK4thDerivativeInPoint;
        break;
    }

    qreal C = 0.0;

    foreach(qreal xi, *samples)
    {
        foreach(qreal xj, *samples)
        {
            C += (this->*xsithKDerivative)((xi - xj)/smoothingParameter);
        }
    }

    C /= qPow(samples->size(), 2);
    C /= qPow(smoothingParameter, xsi+1);

    return C;
}

qreal pluginSmoothingParameterCounter::countSmallC(int xsi)
{
    // Page 81 of Kernel Estimators in System Analysis

    qreal c = 1.0;
    for(int i = 1; i <= xsi; ++i) c *= qreal(i);
    for(int i = 1; i <= (xsi/2); ++i) c /= qreal(i);
    c /= qSqrt(M_PI);
    c /= qPow(2.0 * countStandardDeviationEstimator(), xsi + 1);

    return c;
}

qreal pluginSmoothingParameterCounter::count1stRankPluginSmoothingParameter(qreal h2)
{
    // Page 82 of Kernel Estimators in System Analysis

    qreal h1 = -2.0;
    h1 *= countK4thDerivativeInPoint(0);
    h1 /= U;
    h1 /= countCapitalC(6.0, h2);
    h1 /= samples->size();
    h1 = qPow(h1, 1.0/7.0);

    return h1;
}

qreal pluginSmoothingParameterCounter::countPluginSmoothingParameter(qreal h2)
{
    qreal h1 = count1stRankPluginSmoothingParameter(h2);

    qreal Zf = countCapitalC(4, h1);

    qreal h0 = 0.6;
    h0 /= qPow(0.2, 2.0);
    h0 /= samples->size();
    h0 /= Zf;
    h0 = qPow(h0, 0.2);

    return h0;
}

qreal pluginSmoothingParameterCounter::countStandardDeviationEstimator()
{
    // Page 38, Kernel Estimators in system analysis, P. Kulczycki

    if(samples == NULL)
    {
        qDebug() << "Samples pointer is numm.";
        return -1.0;
    }

    if(samples->size() < 0)
    {
        qDebug() << "Samples size < 0 in smoothingParameterCounter.";
        return -1.0;
    }

    qreal substractor = 0.0, V = 0.0;

    foreach(qreal sample, *samples)
    {
        V += qPow(sample, 2.0);
        substractor += sample;
    }

    V /= (samples->size() - 1);

    substractor *= substractor;
    substractor /= samples->size() * (samples->size() - 1);

    V -= substractor;

    return qSqrt(V);
}

qreal pluginSmoothingParameterCounter::countK4thDerivativeInPoint(qreal point)
{
    // Page 83, Kernel Estimators in system analysis, P. Kulczycki

    qreal result = qPow(point, 4.0);
    result -= 6.0 * qPow(point, 2.0);
    result += 3.0;
    result /= qSqrt(2.0 * M_PI);
    result *= qExp(- 0.5 * qPow(point,2));

    return result;
}

qreal pluginSmoothingParameterCounter::countK6thDerivativeInPoint(qreal point)
{
    // Page 83, Kernel Estimators in system analysis, P. Kulczycki

    qreal result = qPow(point, 6.0);
    result -= 15 * qPow(point, 4.0);
    result += 45.0 * qPow(point, 2.0);
    result -= 15.0;
    result /= qSqrt(2.0 * M_PI);
    result *= qExp(- 0.5 * qPow(point,2));

    return result;
}

qreal pluginSmoothingParameterCounter::countK8thDerivativeInPoint(qreal point)
{
    // Page 83, Kernel Estimators in system analysis, P. Kulczycki

    qreal result = qPow(point, 8);
    result -= 28.0 * qPow(point, 6);
    result += 210.0 * qPow(point, 4);
    result -= 420.0 * qPow(point, 2);
    result += 105.0;
    result /= qSqrt(2.0 * M_PI);
    result *= qExp(- 0.5 * qPow(point,2));

    return result;
}
