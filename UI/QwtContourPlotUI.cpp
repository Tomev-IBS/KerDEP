#include "QwtContourPlotUI.h"

QwtContourPlotUI::QwtContourPlotUI(int *currentStep, const int& imagesPeriod,
                                   const int& seed, DESDA *DESDAAlgorithm,
                                   double* L1Error, double* L2Error,
                                   double* supError, double* modError, double* actual_l1_error, double* actual_l2_error,
                                   double* actual_sup_error, double* actual_mod_error, QDateTime *date_time, double* level_density,
                                   const QString &experiment_description)
 : _currentStep(currentStep), _L1Error(L1Error), _L2Error(L2Error),
   _supError(supError), _modError(modError), _DESDAAlgorithm(DESDAAlgorithm), actual_l1_(actual_l1_error),
   actual_l2_(actual_l2_error), actual_mod_(actual_mod_error), actual_sup_(actual_sup_error),
   date_time_(date_time), level_density_(level_density), experiment_description_(experiment_description)
{
  // Set up the QwtTexts.
  QFont uiFont("Courier New");
  uiFont.setStyleHint(QFont::TypeWriter);
  uiFont.setPointSize(16);

  colored_column_text_.setColor(Qt::red);
  _leftColumnText.setColor(Qt::black);
  right_column_.setColor(Qt::black);

  _leftColumnText.setFont(uiFont);
  colored_column_text_.setFont(uiFont);
  right_column_.setFont(uiFont);

  _leftColumnText.setRenderFlags(Qt::AlignLeft | Qt::AlignTop);
  colored_column_text_.setRenderFlags(Qt::AlignLeft | Qt::AlignTop);
  right_column_.setRenderFlags(Qt::AlignRight | Qt::AlignTop);

  // Set the constant strings.
  _imagesPeriodString = "iw   = " + QString::number(imagesPeriod) + "\n";
  _seedString         = "seed = " + QString::number(seed) + "\n";
  _levelsString       = "lvls = " + QString::number(*level_density_) + "\n";
  _mKPSSString        = "mKPSS = " + QString::number(_DESDAAlgorithm->_kpssM) + "\n";
  _mMaxString         = "m0    = " + QString::number(_DESDAAlgorithm->_maxM) + "\n";
  _mMinString         = "mmin  = " + QString::number(_DESDAAlgorithm->_minM) + "\n";
}

void QwtContourPlotUI::attach(QwtPlot *plot)
{
  /**
   * A method attaching both columns to the specified plot.
   */
  _leftColumnLabel.attach(plot);
  colored_column_label_.attach(plot);
  right_column_label_.attach(plot);
}

void QwtContourPlotUI::updateTexts()
{
  /**
   * This method updates both columns texts.
   */
  updateLeftColumnText();
  updateRightColumnText();
  updateColoredColumnText();
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
  leftColumnText += experiment_description_;
  leftColumnText += "\n\n";

  if(!should_print_errors){
    leftColumnText += QLocale(QLocale::English).toString(*date_time_, "dd MMM yyyy, hh:mm");
    leftColumnText += "\n";
  }

  leftColumnText += "t         = " + QString::number(*_currentStep) + "\n";
  //leftColumnText += _imagesPeriodString;
  //leftColumnText += _levelsString;
  //leftColumnText += _seedString;
  leftColumnText += "\n";
  leftColumnText += "KPSS      = " + formatNumberForDisplay(_DESDAAlgorithm->getStationarityTestValue()) + "\n";
  leftColumnText += "sgmKPSS   = " + formatNumberForDisplay(_DESDAAlgorithm->_sgmKPSS) + "\n";
  leftColumnText += "\n";
  //leftColumnText += _mKPSSString;
  //leftColumnText += _mMaxString;
  //leftColumnText += _mMinString;
  leftColumnText += "m         = " + QString::number(_DESDAAlgorithm->_m) + "\n";
  leftColumnText += "\n";
  //leftColumnText += "beta0 =" + formatNumberForDisplay(_DESDAAlgorithm->_beta0) + "\n";
  //leftColumnText += "\n";
  leftColumnText += "r         = " + formatNumberForDisplay(_DESDAAlgorithm->_r) + "\n";
  leftColumnText += "q         = " + formatNumberForDisplay(_DESDAAlgorithm->_quantileEstimator) + "\n";
  leftColumnText += "#atypical = " + QString::number(_DESDAAlgorithm->_rareElementsNumber) + "\n";

  if(should_print_errors){
    leftColumnText += "\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
    leftColumnText += "";
    leftColumnText += "estimated";
    //leftColumnText += "\n\n";
    //leftColumnText += "L^2      = " + formatNumberForDisplay(*_L2Error);
  }

  /*
  if(_DESDAAlgorithm->_smoothingParametersVector.size() == 2) {
    leftColumnText += "\n";
    leftColumnText += "h1 = " + formatNumberForDisplay(_DESDAAlgorithm->_smoothingParametersVector[0]) + "\n";
    leftColumnText += "h2 = " + formatNumberForDisplay(_DESDAAlgorithm->_smoothingParametersVector[1]);
  }
  */

  _leftColumnText.setText(leftColumnText);
  _leftColumnLabel.setText(_leftColumnText);
}

void QwtContourPlotUI::updateColoredColumnText()
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

  QString colored_column_text = "";

  if(should_print_errors){
    colored_column_text = "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\ntheoretical";
    colored_column_text_.setText(colored_column_text);
    colored_column_label_.setText(colored_column_text_);
    return;
  }

  /*
  colored_column_text += "L1   = " + formatNumberForDisplay(*_L1Error) + "\n";
  colored_column_text += "L1a  = " + formatNumberForDisplay(*actual_l1_) + "\n";
  colored_column_text += "\n";
  colored_column_text += "L2   = " + formatNumberForDisplay(*_L2Error) + "\n";
  colored_column_text += "L2a  = " + formatNumberForDisplay(*actual_l2_) + "\n";
  colored_column_text += "\n";
  colored_column_text += "sup  = " + formatNumberForDisplay(*_supError) + "\n";
  colored_column_text += "supa = " + formatNumberForDisplay(*actual_sup_) + "\n";
  colored_column_text += "\n";
  colored_column_text += "mod  = " + formatNumberForDisplay(*_modError) + "\n";
  colored_column_text += "moda = " + formatNumberForDisplay(*actual_mod_) + "\n";
   //*/

  colored_column_text_.setText(colored_column_text);
  colored_column_label_.setText(colored_column_text_);
}

QString QwtContourPlotUI::formatNumberForDisplay(const double& number)
{
  // According to PK the number should be displayed as #.######
  QString result = "";

  if(number < 0) result = "";

  QStringList splitNumber = QString::number(number, 'f', 7).split(".");
  result += splitNumber[0];

  if(splitNumber.size() == 1) return result;

  result += ".";

  for(int i = 0; i < 3 && i < splitNumber[1].size(); ++i)
    result += splitNumber[1][i];

  return result;
}

void QwtContourPlotUI::SetErrorsPrinting(const bool &should_print_errors) {
  this->should_print_errors = should_print_errors;
}

void QwtContourPlotUI::updateRightColumnText() {
  QString right_column_text = "";

  if(should_print_errors){
    right_column_text = "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nL2 = " + formatNumberForDisplay(*_L2Error);
  }

  right_column_.setText(right_column_text);
  right_column_label_.setText(right_column_);
}
