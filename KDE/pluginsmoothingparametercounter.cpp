#include "pluginsmoothingparametercounter.h"

#include <QDebug>

pluginSmoothingParameterCounter::pluginSmoothingParameterCounter(){}

pluginSmoothingParameterCounter::pluginSmoothingParameterCounter(
    QVector<qreal> *samples, int rank) : samples(samples), rank(rank){}

double pluginSmoothingParameterCounter::countSmoothingParameterValue()
{
  switch(rank)
  {
    case 4:
      return count4thRankPluginSmoothingParameter();
    case 3:
      return count3rdRankPluginSmoothingParameter();
    case 2:
      return count2ndRankPluginSmoothingParameter();
    case 0:
    default:
      return count0RankPluginSmoothingParameter();
  }

  return -1.0;
}

void pluginSmoothingParameterCounter::setSamples(QVector<qreal>* samples)
{
    this->samples = samples;
}

qreal pluginSmoothingParameterCounter::count4thRankPluginSmoothingParameter()
{
  qreal h4 = countH4();
  qreal h3 = countH3(h4);
  qreal h2 = countH2(h3);
  qreal h1 = countH1(h2);
  return countPluginSmoothingParameter(h1);
}

qreal pluginSmoothingParameterCounter::count3rdRankPluginSmoothingParameter()
{
  qreal h3 = countH3();
  qreal h2 = countH2(h3);
  qreal h1 = countH1(h2);
  return countPluginSmoothingParameter(h1);
}

qreal pluginSmoothingParameterCounter::count2ndRankPluginSmoothingParameter()
{
  qreal h2 = countH2();
  qreal h1 = countH1(h2);
  return countPluginSmoothingParameter(h1);
}

qreal pluginSmoothingParameterCounter::count1stRankPluginSmoothingParameter()
{
  qreal h1 = countH1();
  return countPluginSmoothingParameter(h1);
}

qreal pluginSmoothingParameterCounter::count0RankPluginSmoothingParameter()
{
  // Calculation for normal kernel!
  double h = pow(4 * M_PI / (3 * samples->size()), 0.2);
  h *= countStandardDeviationEstimator();
  return h;
}

inline bool pluginSmoothingParameterCounter::isNearlyEqual(double x, double y)
{
  const double epsilon = 1e-5;
  return std::abs(x - y) <= epsilon * std::abs(x);
}

qreal pluginSmoothingParameterCounter::countCapitalC(int xsi, qreal smoothingParameter)
{
    // Page 81 of Kernel Estimators in System Analysis, P. Kulczycki

    qreal (pluginSmoothingParameterCounter::* xsithKDerivative)(qreal);

    switch(xsi)
    {
      case 10:
        xsithKDerivative = &pluginSmoothingParameterCounter::countK10thDerivativeInPoint;
        break;
      case 8:
        xsithKDerivative = &pluginSmoothingParameterCounter::countK8thDerivativeInPoint;
        break;
      case 6:
        xsithKDerivative = &pluginSmoothingParameterCounter::countK6thDerivativeInPoint;
        break;
      case 4:
      default:
        xsithKDerivative = &pluginSmoothingParameterCounter::countK4thDerivativeInPoint;
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

qreal pluginSmoothingParameterCounter::countPluginSmoothingParameter(qreal h1)
{
    qreal Zf = 0;

    if(! isNearlyEqual(h1, 0)){
      Zf = countCapitalC(4, h1);
    }
    else{
      Zf = countSmallC(4);
    }

    qreal h0 = 0.5 / qSqrt(M_PI);
    h0 /= qPow(1, 2.0);
    h0 /= samples->size();
    h0 /= Zf;

    h0 = qPow(h0, 1.0/5.0);

    return h0;
}

qreal pluginSmoothingParameterCounter::countStandardDeviationEstimator()
{
    // Page 38, Kernel Estimators in system analysis, P. Kulczycki

    if(samples == NULL)
    {
        qDebug() << "Samples pointer is null.";
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
    substractor /= samples->size();
    substractor /= (samples->size() - 1); // Has to be separated for large ints.

    V -= substractor;

    return qSqrt(V);
}

qreal pluginSmoothingParameterCounter::countH4(qreal h5)
{
  qreal h4 = -2 * countK10thDerivativeInPoint(0);
  h4 /= U;
  if(! isNearlyEqual(h5, 0)){
    h4 /= countCapitalC(12.0, h5);
  }
  else{
    h4 /= countSmallC(12);
  }
  h4 /= samples->size();

  if(h4 < 0)
  {
      h4 = qPow(qAbs(h4), 1.0/13.0);
      h4 = -h4;
  }
  else
  {
      h4 = qPow(h4, 1.0/13.0);
  }
  return h4;
}

qreal pluginSmoothingParameterCounter::countH3(qreal h4)
{
  qreal h3 = 2.0; // In the original work there's -2, but it generates negative h.
  h3 *= countK8thDerivativeInPoint(0);
  h3 /= U;
  if(! isNearlyEqual(h4, 0)){
    h4 /= countCapitalC(10.0, h4);
  }
  else{
    h3 /= countSmallC(10);
  }
  h3 /= samples->size();

  // qPow has a problem with getting roots from negative numbers,
  // hence if must be used

  if(h3 < 0)
  {
      h3 = qPow(qAbs(h3), 1.0/11.0);
      h3 = -h3;
  }
  else
  {
      h3 = qPow(h3, 1.0/11.0);
  }
  return h3;
}

qreal pluginSmoothingParameterCounter::countH2(qreal h3)
{
  qreal h2 = -2.0;
  h2 *= countK6thDerivativeInPoint(0);
  h2 /= U;
  if(! isNearlyEqual(h3, 0)){
    h2 /= countCapitalC(8.0, h3);
  }
  else{
    h2 /= countSmallC(8);
  }
  h2 /= samples->size();

  // Same issue with qPow as in h3 case

  if(h2 < 0)
  {
      h2 = qPow(qAbs(h2), 1.0/9.0);
      h2 = -h2;
  }
  else
  {
      h2 = qPow(h2, 1.0/9.0);
  }

  return h2;
}

qreal pluginSmoothingParameterCounter::countH1(qreal h2)
{
  qreal h1 = -2.0;
  h1 *= countK4thDerivativeInPoint(0);
  h1 /= U;
  if(! isNearlyEqual(h2, 0)){
    h1 /= countCapitalC(6.0, h2);
  }
  else{
    h1 /= countSmallC(6.0);
  }
  h1 /= samples->size();

  // Same issue with qPow as in h3 case

  if(h2 < 0)
  {
      h2 = qPow(qAbs(h2), 1.0/9.0);
      h2 = -h2;
  }
  else
  {
      h2 = qPow(h2, 1.0/9.0);
  }

  return h2;
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

qreal pluginSmoothingParameterCounter::countK10thDerivativeInPoint(qreal point)
{
  qreal result = qPow(point, 10);
  result -= 45.0 * qPow(point, 8);
  result += 630.0 * qPow(point, 6);
  result -= 3150.0 * qPow(point, 4);
  result += 4725.0 * qPow(point, 2);
  result -= 945.0;
  result /= qSqrt(2.0 * M_PI);
  result *= qExp(- 0.5 * qPow(point,2));

  return result;
}
