#ifndef QWTCONTOURPLOTUI_H
#define QWTCONTOURPLOTUI_H

#include "DESDA.h"

#include <qwt_plot.h>
#include <qwt_plot_textlabel.h>

class QwtContourPlotUI
{
  public:
    QwtContourPlotUI(int *currentStep, const int& imagesPeriod, const int& seed,
                     DESDA *DESDAAlgorithm, double* L1Error, double* L2Error,
                     double* supError, double* modError);
    void updateTexts();
    void attach(QwtPlot* plot);
  private:
    QwtText _leftColumnText;
    QwtText _rightColumnText;
    QwtPlotTextLabel _rightColumnLabel;
    QwtPlotTextLabel _leftColumnLabel;

    int *_currentStep;
    double *_L1Error;
    double *_L2Error;
    double *_supError;
    double *_modError;
    QString _imagesPeriodString;
    QString _seedString;
    QString _mKPSSString;
    QString _mMaxString;
    QString _mMinString;

    DESDA *_DESDAAlgorithm;

    void updateLeftColumnText();
    void updateRightColumnText();
    QString formatNumberForDisplay(const double& number);
};

#endif // QWTCONTOURPLOTUI_H
