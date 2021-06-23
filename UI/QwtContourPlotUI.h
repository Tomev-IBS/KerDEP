#ifndef QWTCONTOURPLOTUI_H
#define QWTCONTOURPLOTUI_H

#include "DESDA.h"
#include <QDateTime>

#include <qwt_plot.h>
#include <qwt_plot_textlabel.h>

class QwtContourPlotUI
{
  public:
    QwtContourPlotUI(int *currentStep, const int& imagesPeriod, const int& seed,
                     DESDA *DESDAAlgorithm, double* L1Error, double* L2Error,
                     double* supError, double* modError, double* actual_l1_error, double* actual_l2_error,
                     double* actual_sup_error, double* actual_mod_error, QDateTime* date_time, double* level_density,
                     const QString &experiment_description);
    void updateTexts();
    void attach(QwtPlot* plot);
    void SetErrorsPrinting(const bool &should_print_errors);
  private:
    QwtText _leftColumnText;
    QwtText _rightColumnText;
    QwtPlotTextLabel _rightColumnLabel;
    QwtPlotTextLabel _leftColumnLabel;

    bool should_print_errors = true;

    int *_currentStep;
    double *_L1Error;
    double *_L2Error;
    double *_supError;
    double *_modError;
    double *actual_l1_;
    double *actual_l2_;
    double *actual_sup_;
    double *actual_mod_;
    double *level_density_;
    QDateTime *date_time_;
    QString _imagesPeriodString;
    QString _seedString;
    QString _levelsString;
    QString _mKPSSString;
    QString _mMaxString;
    QString _mMinString;
    QString experiment_description_;

    DESDA *_DESDAAlgorithm;

    void updateLeftColumnText();
    void updateRightColumnText();
    QString formatNumberForDisplay(const double& number);
};

#endif // QWTCONTOURPLOTUI_H
