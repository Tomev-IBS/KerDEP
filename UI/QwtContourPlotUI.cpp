#include "QwtContourPlotUI.h"

QwtContourPlotUI::QwtContourPlotUI(int *currentStep, const int& imagesPeriod,
                                   const int& seed, DESDA *DESDAAlgorithm,
                                   double* L1Error, double* L2Error,
                                   double* supError, double* modError, double* actual_l1_error, double* actual_l2_error,
                                   double* actual_sup_error, double* actual_mod_error, QDateTime *date_time)
 : _currentStep(currentStep), _L1Error(L1Error), _L2Error(L2Error),
   _supError(supError), _modError(modError), _DESDAAlgorithm(DESDAAlgorithm), actual_l1_(actual_l1_error),
   actual_l2_(actual_l2_error), actual_mod_(actual_mod_error), actual_sup_(actual_sup_error),
   date_time_(date_time)
{
  // Set up the QwtTexts.
  QFont uiFont("Courier New");
  uiFont.setStyleHint(QFont::TypeWriter);
  uiFont.setPointSize(16);

  _leftColumnText.setRenderFlags(Qt::AlignLeft | Qt::AlignTop);
  _rightColumnText.setRenderFlags(Qt::AlignRight | Qt::AlignTop);
  _leftColumnText.setFont(uiFont);
  _rightColumnText.setFont(uiFont);

  // Set the constant strings.
  _imagesPeriodString = "iw   = " + QString::number(imagesPeriod) + "\n";
  _seedString = "seed = " + QString::number(seed) + "\n";
  _mKPSSString = "mKPSS = " + QString::number(_DESDAAlgorithm->_kpssM) + "\n";
  _mMaxString = "m0    = " + QString::number(_DESDAAlgorithm->_maxM) + "\n";
  _mMinString = "mmin  = " + QString::number(_DESDAAlgorithm->_minM) + "\n";
}

void QwtContourPlotUI::attach(QwtPlot *plot)
{
  /**
   * A method attaching both columns to the specified plot.
   */
  _leftColumnLabel.attach(plot);
  _rightColumnLabel.attach(plot);
}

void QwtContourPlotUI::updateTexts()
{
  /**
   * This method updates both columns texts.
   */
  updateLeftColumnText();
  updateRightColumnText();
}

void QwtContourPlotUI::updateLeftColumnText()
{
  /**
    * Method for updating left column texts. Expected values are:
    * i -> current step
    * iw -> image generation interval
    * seed -> current seed
    * / empty line /
    * KPSS -> current KPSS test value
    * sgmKPSS -> current sgmKPSS value
    * / empty line /
    * mKPSS -> m value for KPSS test (size of sample for KPSS)
    * m0 -> m0 value (max reservoir size)
    * mmin -> mMin value (min reservoir size)
    * m -> current effective reservoir size
    * / empty line /
    * beta0 -> current beta0 value
    * / empty line /
    * r ->  r value for rare elements calculation
    * q -> q value for rare elements calculation
    * rare -> rare elements number
    * hx -> smoothing parameter for first dimension
    * hy -> smoothing parameter for second dimension
    */

  QString leftColumnText = "";
  leftColumnText += date_time_->toString() + "\n";
  leftColumnText += "i    = " + QString::number(*_currentStep) + "\n";
  leftColumnText += _imagesPeriodString;
  leftColumnText += _seedString;
  leftColumnText += "\n";
  leftColumnText += "KPSS    =" + formatNumberForDisplay(_DESDAAlgorithm->getStationarityTestValue()) + "\n";
  leftColumnText += "sgmKPSS =" + formatNumberForDisplay(_DESDAAlgorithm->_sgmKPSS) + "\n";
  leftColumnText += "\n";
  leftColumnText += _mKPSSString;
  leftColumnText += _mMaxString;
  leftColumnText += _mMinString;
  leftColumnText += "m     = " + QString::number(_DESDAAlgorithm->_m) + "\n";
  leftColumnText += "\n";
  leftColumnText += "beta0 =" + formatNumberForDisplay(_DESDAAlgorithm->_beta0) + "\n";
  leftColumnText += "\n";
  leftColumnText += "r    = " + formatNumberForDisplay(_DESDAAlgorithm->_r) + "\n";
  leftColumnText += "q    = " + formatNumberForDisplay(_DESDAAlgorithm->_quantileEstimator) + "\n";
  leftColumnText += "rare =  " + QString::number(_DESDAAlgorithm->_rareElementsNumber) + "\n";

  if(_DESDAAlgorithm->_smoothingParametersVector.size() == 2) {
    leftColumnText += "\n";
    leftColumnText += "h1 = " + formatNumberForDisplay(_DESDAAlgorithm->_smoothingParametersVector[0]) + "\n";
    leftColumnText += "h2 = " + formatNumberForDisplay(_DESDAAlgorithm->_smoothingParametersVector[1]);
  }

  _leftColumnText.setText(leftColumnText);
  _leftColumnLabel.setText(_leftColumnText);
}

void QwtContourPlotUI::updateRightColumnText()
{
  /**
    * Method for updating right column texts. Right column holds
    * distances from model plot. We've used L1, L2, sup and mod distances with
    * several lower indexes denoting:
    * _w - estimator build on m0 objects
    * _m - estimator build on m objects
    * _d - estimator build on m weighted objects
    * _p - estimator build on m weighted objects with prediction
    * _n - estimator build on m weighted objects with prediction and rare elements
    *
    * For now we're only using _n estimator.
    */

  QString rightColumnText = "";

  rightColumnText += "L1   = " + formatNumberForDisplay(*_L1Error) + "\n";
  rightColumnText += "L1a  = " + formatNumberForDisplay(*actual_l1_) + "\n";
  rightColumnText += "\n";
  rightColumnText += "L2   = " + formatNumberForDisplay(*_L2Error) + "\n";
  rightColumnText += "L2a  = " + formatNumberForDisplay(*actual_l2_) + "\n";
  rightColumnText += "\n";
  rightColumnText += "sup  = " + formatNumberForDisplay(*_supError) + "\n";
  rightColumnText += "supa = " + formatNumberForDisplay(*actual_sup_) + "\n";
  rightColumnText += "\n";
  rightColumnText += "mod  = " + formatNumberForDisplay(*_modError) + "\n";
  rightColumnText += "moda = " + formatNumberForDisplay(*actual_mod_) + "\n";

  _rightColumnText.setText(rightColumnText);
  _rightColumnLabel.setText(_rightColumnText);
}

QString QwtContourPlotUI::formatNumberForDisplay(const double& number)
{
  // According to PK the number should be displayed as #.######
  QString result = " ";

  if(number < 0) result = "";

  QStringList splitNumber = QString::number(number, 'f', 7).split(".");
  result += splitNumber[0];

  if(splitNumber.size() == 1) return result;

  result += ".";

  for(int i = 0; i < 6 && i < splitNumber[1].size(); ++i)
    result += splitNumber[1][i];

  return result;
}
