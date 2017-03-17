#ifndef PLUGINSMOOTHINGPARAMETERCOUNTER_H
#define PLUGINSMOOTHINGPARAMETERCOUNTER_H

#include <QObject>
#include <QtMath>
#include <QVector>

class pluginSmoothingParameterCounter
{
    public:
        pluginSmoothingParameterCounter(QVector<qreal>* samples);

        qreal count2ndRankPluginSmoothingParameter();
        qreal count3rdRankPluginSmoothingParameter();

    private:

        QVector<qreal>* samples;

        qreal countCapitalC(int xsi, qreal smoothingParameter);
        qreal countSmallC(int xsi);
        qreal count1stRankPluginSmoothingParameter(qreal h2);
        qreal countPluginSmoothingParameter(qreal h2);
        qreal countStandardDeviationEstimator();

        qreal countK4thDerivativeInPoint(qreal point);
        qreal countK6thDerivativeInPoint(qreal point);
        qreal countK8thDerivativeInPoint(qreal point);

        const qreal U = 1.0;
};

#endif // PLUGINSMOOTHINGPARAMETERCOUNTER_H
