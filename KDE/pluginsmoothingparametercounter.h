#ifndef PLUGINSMOOTHINGPARAMETERCOUNTER_H
#define PLUGINSMOOTHINGPARAMETERCOUNTER_H

#include <QObject>
#include <QtMath>
#include <QVector>

#include "smoothingparametercounter.h"

// Page 80 of Kernel Estimators in System Analysis.

class pluginSmoothingParameterCounter : public smoothingParameterCounter
{
public:
  pluginSmoothingParameterCounter();
  pluginSmoothingParameterCounter(QVector<qreal> *samples, int rank);

  double countSmoothingParameterValue();

  void setSamples(QVector<qreal>* samples);

  qreal count4thRankPluginSmoothingParameter();
  qreal count3rdRankPluginSmoothingParameter();
  qreal count2ndRankPluginSmoothingParameter();
  qreal count1stRankPluginSmoothingParameter();

private:

  QVector<qreal>* samples;

  int rank;

  bool isNearlyEqual(double x, double y);

  qreal countCapitalC(int xsi, qreal smoothingParameter);
  qreal countSmallC(int xsi);
  qreal countPluginSmoothingParameter(qreal h1 = 0);
  qreal countStandardDeviationEstimator();

  qreal countH4(qreal h5 = 0);
  qreal countH3(qreal h4 = 0);
  qreal countH2(qreal h3 = 0);
  qreal countH1(qreal h2 = 0);

  qreal countK4thDerivativeInPoint(qreal point);
  qreal countK6thDerivativeInPoint(qreal point);
  qreal countK8thDerivativeInPoint(qreal point);
  qreal countK10thDerivativeInPoint(qreal point);

  const qreal U = 1.0;
};

#endif // PLUGINSMOOTHINGPARAMETERCOUNTER_H
