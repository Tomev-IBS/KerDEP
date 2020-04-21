#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <math.h>
#include <climits>
#include <float.h>
#include <set>
#include <QDebug>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "UI/plotLabel.h"

#include "Functions/multivariatenormalprobabilitydensityfunction.h"
#include "Functions/complexfunction.h"

#include "Distributions/distributions.h"
#include "KDE/pluginsmoothingparametercounter.h"
#include "KDE/weightedSilvermanSmoothingParameterCounter.h"

#include "Reservoir_sampling/biasedReservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/basicReservoirSamplingAlgorithm.h"

#include "Reservoir_sampling/distributiondataparser.h"
#include "Reservoir_sampling/progressivedistributiondatareader.h"

#include "groupingThread/groupingThread.h"

#include "groupingThread/kMedoidsAlgorithm/numericalAttributeData.h"

#include "DESDA.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    kernelTypes << "Normal" << "Triangle" << "Epanecznikow" << "Dull";

    ui->setupUi(this);

    setupValidators();
    setupPlot();
    setupKernelsTable();

    testNewFunctionalities();

    // Set font
    //auto appFont = QFont("Monospace");
    auto appFont = QFont("Courier New");
    appFont.setStyleHint(QFont::TypeWriter);
    setFont(appFont);
}

void MainWindow::testNewFunctionalities()
{
  qDebug() << "Start test.";

  qDebug() << "Nothing to test.";

  qDebug() << "Finish test.";
}

QVector<double> MainWindow::getTargetFunctionValuesOnDomain(QVector<double> *domain)
{
  QVector<double> values = {};
  for(auto val : *domain){
    std::vector<double> x = {val};
    values.push_back(targetFunction->getValue(&x));
  }
  return values;
}

double MainWindow::calculateL1Error(const QVector<double> &model, const QVector<double> estimated, const QVector<double> &domain)
{
  QVector<qreal> errorHolder = {};

  for(int i = 0; i < model.size(); ++i)
    errorHolder.push_back(fabs(model[i] - estimated[i]));

  double avg = 0;

  for(auto val : errorHolder)
    avg += val;

  if(errorHolder.size() > 0) avg /= errorHolder.size();

  double min = domain[0];
  double max = domain[domain.size() - 1];
  double len = max - min;

  return avg * len;
}

double MainWindow::calculateL2Error(const QVector<double> &model, const QVector<double> estimated, const QVector<double> &domain)
{
  QVector<qreal> errorHolder = {};

  for(int i = 0; i < model.size(); ++i)
    errorHolder.push_back(pow(model[i] - estimated[i], 2));

  double avg = 0;

  for(auto val : errorHolder)
    avg += val;

  if(errorHolder.size() > 0) avg /= errorHolder.size();

  double min = domain[0];
  double max = domain[domain.size() - 1];
  double len = max - min;

  return avg * len;
}

double MainWindow::calculateSupError(const QVector<double> &model, const QVector<double> estimated)
{
  double sup = fabs(model[0] - estimated[0]);

  for(int i = 0; i < model.size(); ++i){
    double val = fabs(model[i] - estimated[i]);
    sup = val > sup ? val : sup;
  }

  return sup;
}

double MainWindow::findExtrema(const QVector<double> &values, const QVector<double> domain)
{
  double extrema = domain[0];
  double maxVal = values[0];

  for(int i = 0; i < values.size(); ++i){
    double val = values[i];

    if(val > maxVal){
      maxVal = val;
      extrema = domain[i];
    }
  }

  return extrema;
}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::setupValidators()
{
    QLocale locale = QLocale::English;
    locale.setNumberOptions(QLocale::c().numberOptions());

    std::unique_ptr<QIntValidator> naturalNumbersValidator(new QIntValidator(0, std::numeric_limits<int>::max(), this));
    std::unique_ptr<QIntValidator> positiveNaturalNumbersValidator(new QIntValidator(1, std::numeric_limits<int>::max(), this));

    std::unique_ptr<QDoubleValidator> xAxisValidator(new QDoubleValidator(MIN_X, MAX_X, DECIMAL_NUMBERS, this));
    xAxisValidator->setLocale(locale);
    xAxisValidator->setNotation(QDoubleValidator::StandardNotation);

    std::unique_ptr<QDoubleValidator> yAxisValidator(new QDoubleValidator(MIN_Y, MAX_Y, DECIMAL_NUMBERS, this));
    yAxisValidator->setLocale(locale);
    yAxisValidator->setNotation(QDoubleValidator::StandardNotation);

    ui->lineEdit_sampleSize->setValidator(positiveNaturalNumbersValidator.get());
    ui->lineEdit_seed->setValidator(naturalNumbersValidator.get());
    ui->lineEdit_milisecondsDelay->setValidator(naturalNumbersValidator.get());

    ui->lineEdit_minX->setValidator(xAxisValidator.get());
    ui->lineEdit_maxX->setValidator(xAxisValidator.get());

    ui->lineEdit_minY->setValidator(yAxisValidator.get());
    ui->lineEdit_maxY->setValidator(yAxisValidator.get());
}

void MainWindow::setupPlot()
{
    //ui->widget_plot->xAxis->setLabel("x");
    //ui->widget_plot->yAxis->setLabel("f(x)");

    ui->widget_plot->xAxis->setRange(DEFAULT_MIN_X, DEFAULT_MAX_X);
    ui->widget_plot->yAxis->setRange(DEFAULT_MIN_Y, DEFAULT_MAX_Y);
}

void MainWindow::setupKernelsTable()
{
    ui->tableWidget_dimensionKernels
      ->horizontalHeader()->setStretchLastSection(true);

    refreshKernelsTable();
    refreshTargetFunctionTable();
}

void MainWindow::drawPlots(DESDA *DESDAAlgorithm)
{
    clearPlot();
    resizePlot();

    // Generate plot of model function
    if(ui->checkBox_showEstimatedPlot->isChecked()){
      QVector<qreal> modelDistributionY = getTargetFunctionValuesOnDomain(&_drawableDomain);
      addPlot(&modelDistributionY, _MODEL_PLOT_PEN);
    }

    // Generate less elements KDE plot (navy blue)
    if(ui->checkBox_showEstimationPlot->isChecked()){
      _lessElementsEstimatorY =
          DESDAAlgorithm->getKDEValues(&_drawableDomain);
      addPlot(&_lessElementsEstimatorY, _KDE_PLOT_PEN);
    }

    // Generate weighted estimator plot (light blue)
    if(ui->checkBox_showWeightedEstimationPlot->isChecked()){
        _weightedEstimatorY =
            DESDAAlgorithm->getWeightedKDEValues(&_drawableDomain);
      addPlot(&_weightedEstimatorY, _WEIGHTED_PLOT_PEN);
    }

    // Generate full estimator plot (BLACK)
    if(ui->checbox_showFullEstimator->isChecked()){
        _windowedEstimatorY =
            DESDAAlgorithm->getWindowKDEValues(&_drawableDomain);
      addPlot(&_windowedEstimatorY, _WINDOWED_PLOT_PEN);
    }

    // Generate plot for kernel prognosis derivative
    if(ui->checkBox_kernelPrognosedPlot->isChecked())
      addPlot(&_kernelPrognosisDerivativeValues, _DERIVATIVE_PLOT_PEN);
      //addPlot(&_kernelPrognosisDerivativeValues, _DERIVATIVE_PLOT_PEN);

    // Generate plot for standarized prognosis derivative, assuming that
    // normal derivative was generated first
    if(ui->checkBox_standarizedDerivative->isChecked()){
      QVector<double> standarizedDerivativeY = {};
      for(auto val : _kernelPrognosisDerivativeValues){
        standarizedDerivativeY.push_back(
          0.1 * val / DESDAAlgorithm->_maxAbsDerivativeValueInCurrentStep
        );
      }
      addPlot(&standarizedDerivativeY, _STANDARIZED_DERIVATIVE_PLOT_PEN);
    }

    if(ui->checkBox_sigmoidallyEnhancedKDE->isChecked()){
        _sigmoidallyEnhancedPlotY =
            DESDAAlgorithm->getEnhancedKDEValues(&_drawableDomain);
        addPlot(&_sigmoidallyEnhancedPlotY, _DESDA_KDE_PLOT_PEN);
    }

    if(ui->checkBox_showUnusualClusters->isChecked()){
        _atypicalElementsValuesAndDerivatives =
            DESDAAlgorithm->getAtypicalElementsValuesAndDerivatives();
        _quantileEstimatorValue = DESDAAlgorithm->_quantileEstimator;
      markUncommonClusters();
    }

    if(ui->checkBox_REESEKDE->isChecked()){
        _rareElementsEnhancedPlotY =
            DESDAAlgorithm->getRareElementsEnhancedKDEValues(&_drawableDomain);
        addPlot(&_rareElementsEnhancedPlotY, _DESDA_RARE_ELEMENTS_KDE_PLOT_PEN);
    }
    // Draw plots
    ui->widget_plot->replot();

    oldKerernelY = KDEEstimationY;
}

void MainWindow::addPlot(const QVector<qreal> *Y, const QPen &pen)
{
  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(_drawableDomain, *Y);
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(pen);
}


void MainWindow::resizePlot()
{
    // Resize plot
    qreal   minX = ui->lineEdit_minX->text().toDouble(),
            maxX = ui->lineEdit_maxX->text().toDouble(), // Standard
            //maxX = 3 + ui->lineEdit_distributionProgression->text().toDouble() * 3000, // For progression
            //maxX = 3 + ui->lineEdit_distributionProgression->text().toDouble(), // /* should jump */
            minY = ui->lineEdit_minY->text().toDouble(),
            maxY = ui->lineEdit_maxY->text().toDouble();

    // Check if sizes are entered correctly
    if(minX < maxX)
    {
        // If so change
        ui->widget_plot->xAxis->setRange(minX, maxX);
    }
    else
    {
        // If not log it and correct
        qDebug() << "Minimal x value cannot be lower than it's maximal value.";
        minX = ui->widget_plot->xAxis->range().minRange;
        maxX = ui->widget_plot->xAxis->range().maxRange;
    }

    if(minY < maxY)
    {
        // If so change
        ui->widget_plot->yAxis->setRange(minY, maxY);
    }
    else
    {
        // If not log it and correct
        qDebug() << "Minimal y value cannot be lower than it's maximal value.";
        minY = ui->widget_plot->yAxis->range().minRange;
        maxY = ui->widget_plot->yAxis->range().maxRange;
    }

   ui->widget_plot->xAxis->setTickLabelFont(QFont(font().family(), 14));
   ui->widget_plot->yAxis->setTickLabelFont(QFont(font().family(), 14));

   QVector<double> ticks;
   QVector<QString> labels;
   // Switch to an array

   for(double i = minX; i <= maxX + 1; i += 1)
   {
     ticks << i;
     labels << QString::number(i);
   }

   QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
   textTicker->addTicks(ticks, labels);

   ui->widget_plot->xAxis->setTicker(textTicker);

}

void MainWindow::clearPlot()
{
    while(ui->widget_plot->graphCount() != 0)
        ui->widget_plot->removeGraph(0);

    for(auto a : _linesOnPlot){
      ui->widget_plot->removeItem(a);
    }

    _linesOnPlot.clear();
}

unsigned long long MainWindow::markUncommonClusters()
{
  for(auto x : _atypicalElementsValuesAndDerivatives){
    // Only works for distribution data
    QCPItemLine *verticalLine = new QCPItemLine(ui->widget_plot);
    verticalLine->start->setCoords(x.first, 0);
    verticalLine->end->setCoords(x.first, -_quantileEstimatorValue);
    if(x.second > 0)
      verticalLine->setPen(QPen(Qt::green));
    else
      verticalLine->setPen(QPen(Qt::red));
    _linesOnPlot.push_back(verticalLine);
  }

  return _atypicalElementsValuesAndDerivatives.size();
}

QString MainWindow::formatNumberForDisplay(double number)
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

void MainWindow::fillStandardDeviations(vector<std::shared_ptr<vector<double>>> *stDevs)
{
    int dimensionsNumber                = ui->spinBox_dimensionsNumber->value(),
        targetFunctionElementsNumber    = ui->tableWidget_targetFunctions->rowCount();

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        stDevs->push_back(std::make_shared<vector<double>>());

        for(int dimensionIndex = 0; dimensionIndex < dimensionsNumber; ++dimensionIndex)
        {
            stDevs->back().get()->push_back
            (
                (static_cast<QLineEdit*>(
                    (static_cast<QTableWidget*>(ui->tableWidget_targetFunctions->cellWidget(functionIndex, STDEV_COLUMN_INDEX)))
                    ->cellWidget(dimensionIndex, 0)
                ))
                ->text().toDouble()
            );
        }
    }
}

void MainWindow::fillMeans(vector<std::shared_ptr<vector<double>>> *means)
{
    int dimensionsNumber = ui->spinBox_dimensionsNumber->value(),
        targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        means->push_back(std::make_shared<vector<double>>());

        for(int dimensionIndex = 0; dimensionIndex < dimensionsNumber; ++dimensionIndex)
        {
            means->back()->push_back
            (
                (static_cast<QLineEdit*>(
                    (static_cast<QTableWidget*>(ui->tableWidget_targetFunctions->cellWidget(functionIndex, MEAN_COLUMN_INDEX)))
                    ->cellWidget(dimensionIndex, 0)
                ))
                ->text().toDouble()
            );
        }
    }
}

void MainWindow::fillDomain(QVector<std::shared_ptr<point>>* domain,  std::shared_ptr<point> *prototypePoint) {
    // Check if domain is nullpointer
    if(domain == nullptr) return;

    std::shared_ptr<point> pPoint;

   // Check if prototype is a null pointer
    if(prototypePoint == nullptr)
    {
        // If so make it a point pointer
        pPoint = std::make_shared<point>();
    }
    else
    {
        pPoint = *prototypePoint;
    }

    qreal val = ui->lineEdit_minX->text().toDouble();
    //qreal maxVal = 3 + ui->lineEdit_distributionProgression->text().toDouble() * 3000;  // Traveling case
    //qreal maxVal = 3 + ui->lineEdit_distributionProgression->text().toDouble(); /* should jump */
    qreal maxVal = ui->lineEdit_maxX->text().toDouble();

    //while(val <= ui->lineEdit_maxX->text().toDouble())
    while(val <= maxVal)
    {
        pPoint.get()->push_back(val);

        if(pPoint.get()->size() == (size_t) ui->spinBox_dimensionsNumber->value())
        {
            domain->append(std::make_shared<point>());

            foreach(qreal dimensionVal, *(pPoint.get()))
            {
                domain->back()->push_back(dimensionVal);
            }
        }
        else
        {
            fillDomain(domain, prototypePoint);
        }

        pPoint->erase(pPoint->end() - 1);

        val += ui->lineEdit_domainDensity->text().toDouble();
    }
}

distribution* MainWindow::generateTargetDistribution(
  vector<std::shared_ptr<vector<double>>> *means,
  vector<std::shared_ptr<vector<double>>> *stDevs)
{
    int seed = ui->lineEdit_seed->text().toInt();

    vector<double> contributions;
    vector<std::shared_ptr<distribution>> elementalDistributions;

    int targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

    //double maxMean = ui->lineEdit_distributionProgression->text().toDouble() *  3000;
    double maxMean = ui->lineEdit_maxX->text().toDouble();

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        contributions.push_back
        (
            (static_cast<QLineEdit*>(ui->tableWidget_targetFunctions
              ->cellWidget(functionIndex, CONTRIBUTION_COLUMN_INDEX)))
              ->text().toDouble()
        );

        elementalDistributions.push_back(
              std::shared_ptr<distribution>(
                new normalDistribution(seed,
                                       (*means)[functionIndex].get(),
                                       (*stDevs)[functionIndex].get(),
                                       maxMean)));
    }

    return new complexDistribution(seed, &elementalDistributions, &contributions);
}

reservoirSamplingAlgorithm *MainWindow::generateReservoirSamplingAlgorithm(dataReader *reader,
                                                                           dataParser *parser)
{
    int sampleSize = ui->lineEdit_sampleSize->text().toInt(),
        stepsNumber = ui->lineEdit_iterationsNumber->text().toInt(),
        samplingAlgorithmID = ui->comboBox_samplingAlgorithm->currentIndex();

    switch(samplingAlgorithmID)
    {
        case BIASED_RESERVOIR_SAMPLING_ALGORITHM:
            return new biasedReservoirSamplingAlgorithm(reader, parser, sampleSize, stepsNumber);
        case BASIC_RESERVOIR_SAMPLING_ALGORITHM:
        default:
            return new basicReservoirSamplingAlgorithm(reader, parser, sampleSize, stepsNumber);
    }
}

kernelDensityEstimator* MainWindow::generateKernelDensityEstimator(
    int dimensionsNumber)
{
    vector<int> kernelsIDs;
    vector<double> smoothingParameters;
    vector<std::string> carriersRestrictions;

    for(int rowNumber = 0; rowNumber < dimensionsNumber; ++rowNumber)
    {
        kernelsIDs.push_back((static_cast<QComboBox*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, KERNEL_COLUMN_INDEX)))->currentIndex());
        smoothingParameters.push_back((static_cast<QLineEdit*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->text().toDouble());
        carriersRestrictions.push_back((static_cast<QLineEdit*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX)))->text().toStdString());
    }

    return new kernelDensityEstimator(
                    &samples,
                    &smoothingParameters,
                    &carriersRestrictions,
                    PRODUCT,
                    &kernelsIDs
                );
}

function* MainWindow::generateTargetFunction(
    vector<std::shared_ptr<vector<double>>>* means,
    vector<std::shared_ptr<vector<double>>>* stDevs)

{
  vector<double> contributions;
  vector<std::shared_ptr<function>> elementalFunctions;

  int targetFunctionElementsNumber = ui->tableWidget_targetFunctions
                                       ->rowCount();

  // Check if contributions are set correctly. If they are, then last
  // contribution is >= 0;
  if(static_cast<QLineEdit*>(ui->tableWidget_targetFunctions
       ->cellWidget(targetFunctionElementsNumber -1, CONTRIBUTION_COLUMN_INDEX)
      )->text().toDouble() <= 0
    )
  {
      // If not then uniform distributions and log error
      qDebug() << "Contributions aren't set correctly. Uniforming "
                  "contributions.";
      uniformContributions();
  }

  for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
  {

      contributions.push_back
      (
          (static_cast<QLineEdit*>(ui->tableWidget_targetFunctions->cellWidget(functionIndex, CONTRIBUTION_COLUMN_INDEX)))
          ->text().toDouble()
      );


      elementalFunctions.push_back(std::shared_ptr<function>(new multivariateNormalProbabilityDensityFunction(means->at(functionIndex).get(),
                                                                                                          stDevs->at(functionIndex).get())));
  }

  return new complexFunction(&contributions, &elementalFunctions);
}

int MainWindow::canAnimationBePerformed(int dimensionsNumber)
{
  switch(dimensionsNumber)
  {
    case 1:
      return 1;
    default:
      qDebug() << "Dimensions number is not equal 1. Animation cannot be performed.";
      return -2;
  }
}

void MainWindow::delay(int ms)
{
    QTime dieTime = QTime::currentTime().addMSecs( ms );

    while( QTime::currentTime() < dieTime )
    {
        QCoreApplication::processEvents( QEventLoop::AllEvents, 100 );
    }
}

void MainWindow::on_spinBox_dimensionsNumber_editingFinished()
{
    refreshKernelsTable();
    refreshTargetFunctionTable();
}

void MainWindow::refreshKernelsTable()
{
    // Get new number of rows
    int newNumberOfRows = ui->spinBox_dimensionsNumber->value();

    // If new number of rows is equal to current number of rows do nothing
    if(newNumberOfRows == ui->tableWidget_dimensionKernels->rowCount()) return;

    // Set new row count
    ui->tableWidget_dimensionKernels->setRowCount(newNumberOfRows);

    QLocale locale = QLocale::English;
    locale.setNumberOptions(QLocale::c().numberOptions());

    QDoubleValidator* smoothingParameterValidator = new QDoubleValidator(MIN_SMOOTHING_P, MAX_SMOOTHING_P, DECIMAL_NUMBERS, this);
    smoothingParameterValidator->setLocale(locale);
    smoothingParameterValidator->setNotation(QDoubleValidator::StandardNotation);

    for(int rowNumber = 0; rowNumber < newNumberOfRows; ++rowNumber) addKernelToTable(rowNumber, smoothingParameterValidator);
}

void MainWindow::addKernelToTable(int rowNumber,
                                  QDoubleValidator* smoothingParameterValidator)
{
    // Add combobox with kernels

    ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, KERNEL_COLUMN_INDEX, new QComboBox());

    (static_cast<QComboBox*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, KERNEL_COLUMN_INDEX)))->insertItems(0, kernelTypes);

    // Add input box with validator for smoothing parameters
    ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX, new QLineEdit());

    (static_cast<QLineEdit*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->setText("1.0");
    (static_cast<QLineEdit*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->setValidator(smoothingParameterValidator);

    // Add input box for carrier restriction value
    ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX, new QLineEdit());

    (static_cast<QLineEdit*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX)))->setText("None.");
}

void MainWindow::refreshTargetFunctionTable()
{
    int numberOfRows = ui->tableWidget_targetFunctions->rowCount(),
        dimensionsNumber = ui->spinBox_dimensionsNumber->value();

    // Ensure that rows number is at least 1
    if(numberOfRows == 0)
        numberOfRows = 1;

    // Set row count
    ui->tableWidget_targetFunctions->setRowCount(numberOfRows);

    QLocale locale = QLocale::English;
    locale.setNumberOptions(QLocale::c().numberOptions());

    QDoubleValidator* meanValidator = new QDoubleValidator(-10.0, 10.0, 3, this);
    meanValidator->setLocale(locale);
    meanValidator->setNotation(QDoubleValidator::StandardNotation);

    QDoubleValidator* stDevValidator = new QDoubleValidator(-5.0, 5.0, 3, this);
    stDevValidator->setLocale(locale);
    stDevValidator->setNotation(QDoubleValidator::StandardNotation);

    QDoubleValidator* contributionValidator = new QDoubleValidator(0.0, 100.0, 3, this);
    contributionValidator->setLocale(locale);
    contributionValidator->setNotation(QDoubleValidator::StandardNotation);

    QTableWidget *targetFunctionTablePointer = static_cast<QTableWidget*>(ui->tableWidget_targetFunctions),
                 *meansTablePointer, *stDevsTablePointer;

    for(int rowIndex = 0; rowIndex < numberOfRows; ++rowIndex)
    {
        // TODO TR: Ensure that this doesn't result in memory leaks
        targetFunctionTablePointer->setCellWidget(rowIndex, MEAN_COLUMN_INDEX, new QTableWidget());

        meansTablePointer = static_cast<QTableWidget*>(ui->tableWidget_targetFunctions->cellWidget(rowIndex, MEAN_COLUMN_INDEX));
        meansTablePointer->setRowCount(dimensionsNumber);
        meansTablePointer->setColumnCount(1);
        meansTablePointer->horizontalHeader()->hide();

        // TODO TR: Ensure that this doesn't result in memory leaks
        targetFunctionTablePointer->setCellWidget(rowIndex, STDEV_COLUMN_INDEX, new QTableWidget());

        stDevsTablePointer = static_cast<QTableWidget*>(ui->tableWidget_targetFunctions->cellWidget(rowIndex, STDEV_COLUMN_INDEX));
        stDevsTablePointer->setRowCount(dimensionsNumber);
        stDevsTablePointer->setColumnCount(1);
        stDevsTablePointer->horizontalHeader()->hide();

        for(int dimensionNumber = 0; dimensionNumber < dimensionsNumber; ++dimensionNumber)
        {
           meansTablePointer->setCellWidget(dimensionNumber, 0, new QLineEdit());
           (static_cast<QLineEdit*>(meansTablePointer->cellWidget(dimensionNumber, 0)))->setText("0.0");
           (static_cast<QLineEdit*>(meansTablePointer->cellWidget(dimensionNumber, 0)))->setValidator(meanValidator);

           stDevsTablePointer->setCellWidget(dimensionNumber, 0, new QLineEdit());
           (static_cast<QLineEdit*>(stDevsTablePointer->cellWidget(dimensionNumber, 0)))->setText("1.0");
           (static_cast<QLineEdit*>(stDevsTablePointer->cellWidget(dimensionNumber, 0)))->setValidator(stDevValidator);
        }

        // TODO TR: Ensure that this doesn't result in memory leaks
        targetFunctionTablePointer->setCellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX, new QLineEdit());
        (static_cast<QLineEdit*>(targetFunctionTablePointer->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->setMaxLength(6);
        (static_cast<QLineEdit*>(targetFunctionTablePointer->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->setValidator(contributionValidator);
        QObject::connect((static_cast<QLineEdit*>(targetFunctionTablePointer->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX))), SIGNAL(textEdited(QString)), this, SLOT(updateLastContribution()));
    }

    // Disable last contribution cell, as it's filled automatically
    (static_cast<QLineEdit*>(targetFunctionTablePointer->cellWidget(numberOfRows -1, CONTRIBUTION_COLUMN_INDEX)))->setEnabled(false);

    uniformContributions();
}

void MainWindow::uniformContributions()
{
    int numberOfRows = ui->tableWidget_targetFunctions->rowCount(),
        lastRowIndex = numberOfRows - 1;

    for(int rowIndex = 0; rowIndex < lastRowIndex; ++rowIndex)
    {
      static_cast<QLineEdit*>(
        ui->tableWidget_targetFunctions
          ->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)
      )->setText(QString::number(100.0/numberOfRows));
    }

    static_cast<QLineEdit*>(
      ui->tableWidget_targetFunctions
        ->cellWidget(lastRowIndex, CONTRIBUTION_COLUMN_INDEX)
    )->setText(QString::number(countLastContribution()));
}

qreal MainWindow::countLastContribution()
{
    qreal result = 100.0;

    int lastRowIndex = ui->tableWidget_targetFunctions->rowCount()-1;

    for(int rowIndex = 0; rowIndex < lastRowIndex; ++rowIndex)
      result -= (static_cast<QLineEdit*>
        (ui->tableWidget_targetFunctions
           ->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX))
      )->text().toDouble();

    return result;
}

void MainWindow::updateLastContribution()
{
    int lastRowIndex = ui->tableWidget_targetFunctions->rowCount()-1;
    qreal lastContributionValue = countLastContribution();

    (static_cast<QLineEdit*>(ui
      ->tableWidget_targetFunctions
      ->cellWidget(lastRowIndex, CONTRIBUTION_COLUMN_INDEX))
    )->setText(QString::number(lastContributionValue));
}

void MainWindow::on_pushButton_addTargetFunction_clicked()
{
    int newRowsNumber = ui->tableWidget_targetFunctions->rowCount() +1;

    // Set row count
    ui->tableWidget_targetFunctions->setRowCount(newRowsNumber);

    refreshTargetFunctionTable();
}

void MainWindow::on_pushButton_removeTargetFunction_clicked()
{
    int newRowsNumber = ui->tableWidget_targetFunctions->rowCount() -1;

    // Set row count
    ui->tableWidget_targetFunctions->setRowCount(newRowsNumber);

    refreshTargetFunctionTable();
}

double MainWindow::numericIntegral(const QVector<qreal> *Y)
{
  double  integral = 0.0,
          domainDensity = ui->lineEdit_domainDensity->text().toDouble();

  for(auto y : *Y)
    integral += y * domainDensity;

  return integral;
}

void MainWindow::on_pushButton_start_clicked()
{
  int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

  if(!canAnimationBePerformed(dimensionsNumber)) return;

  QString seedString = ui->lineEdit_seed->text();

  // Log that application started generating KDE
  // Standard seed was 5625.
  qDebug() << "KDE animation started.";
  qDebug() << "Seed: " + seedString +
              ", Sample size: " + ui->lineEdit_sampleSize->text();

  srand(static_cast<unsigned int>(ui->lineEdit_seed->text().toInt()));

  fillMeans(&means);
  fillStandardDeviations(&stDevs);

  targetFunction.reset(generateTargetFunction(&means, &stDevs));

  std::shared_ptr<kernelDensityEstimator>
    estimator(generateKernelDensityEstimator(dimensionsNumber));

  estimator->_shouldConsiderWeights = false;

  kernelPrognoser.reset(generateKernelDensityEstimator(dimensionsNumber));
  _enchancedKDE.reset(generateKernelDensityEstimator(dimensionsNumber));

  std::shared_ptr<distribution>
    targetDistribution(generateTargetDistribution(&means, &stDevs));

  parser.reset(new distributionDataParser(&attributesData));

  qreal progressionSize =
      ui->lineEdit_distributionProgression->text().toDouble();

  reader.reset(
    new progressiveDistributionDataReader(targetDistribution.get(),
                                          progressionSize,
                                          0 /* delay */,
                                          false /* should jump */)
  );

  reader->gatherAttributesData(&attributesData);
  parser->setAttributesOrder(reader->getAttributesOrder());

  positionalSecondGradeEstimatorCountingMethod =
    ui->comboBox_rareElementsMethod->currentIndex();

  reservoirSamplingAlgorithm* algorithm =
    generateReservoirSamplingAlgorithm(reader.get(), parser.get());

  objects.clear();

  int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();
  int medoidsNumber = 50;

  groupingThread gt(&storedMedoids, parser);

  gt.setAttributesData(&attributesData);

  qDebug() << "Attributes data set.";

  int sampleSize = ui->lineEdit_sampleSize->text().toInt();
  gt.initialize(medoidsNumber, sampleSize);

  _longestStepExecutionInSecs = 0;

  double newWeightB = 0.5;
  int mE = ui->lineEdit_sampleSize->text().toInt() / 2;

  clusters = &storedMedoids;

  weightedSilvermanSmoothingParameterCounter smoothingParamCounter(clusters, 0);

  kernelPrognoser->_shouldConsiderWeights = false;

  DESDA DESDAAlgorithm(
    estimator,
    kernelPrognoser,
    _enchancedKDE,
    ui->lineEdit_weightModifier->text().toDouble(),
    algorithm,
    clusters,
    &storedMedoids,
    ui->lineEdit_rarity->text().toDouble(),
    &gt, newWeightB, mE
  );


  QString expNum = "700";
  this->setWindowTitle("Experiment #" + expNum);
  QString expDesc = "reservoir, v=0.1, beta0=0.2/3, mMin=" +
                      QString::number(DESDAAlgorithm._minM) +", sz001";
  screenGenerationFrequency = 10;

  //QString driveDir = "D:\\Dysk Google\\"; // Home
  QString driveDir = "\\\\beabourg\\private\\"; // WIT PCs

  QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " ("
                    + expDesc + ")\\";

  clearPlot();
  resizePlot();

  plotLabel expNumLabel(ui->widget_plot, 0.02, 0.25, "Exp." + expNum);
  expNumLabel.setFont(QFont("Courier New", 250));

  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  QString imageName = dirPath + QString::number(0) + ".png";

  qDebug() << "Image saved: " << ui->widget_plot->savePng(imageName,
                                                           0, 0, 1, -1);
  expNumLabel.setText("");

  double horizontalOffset = 0.01, verticalOffset = 0.01, verticalStep = 0.03;

  plotLabel iTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "i     = 0");
  verticalOffset += verticalStep;

  plotLabel iwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "iw    = " + QString::number(screenGenerationFrequency));
  verticalOffset += verticalStep;

  plotLabel seedTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "seed  = " + ui->lineEdit_seed->text());
  verticalOffset += verticalStep;

  plotLabel betaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
              "beta0 = " + QString::number(DESDAAlgorithm._beta0));
  verticalOffset += verticalStep;

  plotLabel m0TextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "m0    = " + ui->lineEdit_sampleSize->text());
  verticalOffset += verticalStep;

  plotLabel mMinTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "mmin  = " + QString::number(DESDAAlgorithm._minM));
  verticalOffset += verticalStep;

  plotLabel mTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "m     = " + QString::number(sampleSize));
  verticalOffset += verticalStep;

  plotLabel vTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "v     =" + formatNumberForDisplay(DESDAAlgorithm._v));
  verticalOffset += verticalStep;

  plotLabel rTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "r     =" + formatNumberForDisplay(DESDAAlgorithm._r));
  verticalOffset += verticalStep;

  plotLabel qTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "q     =" + formatNumberForDisplay(DESDAAlgorithm._quantileEstimator));

  verticalOffset += verticalStep;
  plotLabel rareElementsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "rare  = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel KPSSTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "KPSS     = 0");
  verticalOffset += verticalStep;

  plotLabel sgmKPSSTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "sgmKPSS  = 0");
  verticalOffset += verticalStep;

  plotLabel sgmKPSS2TextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "sgmKPSS2 = 0"); // 2 sgmKPSS - 1
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel maxAbsATextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "max(|g|)     = 0");
  verticalOffset += verticalStep;

  std::vector<QString> gTextLabelsLabels = {"g(int(0.2m)) = ", "g(int(0.5m)) = ", "g(int(0.8m)) = "};
  std::vector<plotLabel> gTextLabels = {};

  for(int i = 0; i < gTextLabelsLabels.size(); ++i){
    gTextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset, gTextLabelsLabels[i] + "0"));
    verticalOffset += verticalStep;
  }

  verticalOffset += verticalStep;

  std::vector<plotLabel> wTextLabels = {};
  std::vector<plotLabel> wStarTextLabels = {};
  std::vector<plotLabel> wStar2TextLabels = {};
  std::vector<plotLabel> wStar3TextLabels = {};

  std::vector<QString> wTextLabelsLabels =      {"w(int(0.2m))    = ", "w(int(0.5m))    = ", "w(int(0.8m))    = "};
  std::vector<QString> wStarTextLabelsLabels =  {"w*(int(0.2m))   = ", "w*(int(0.5m))   = ", "w*(int(0.8m))   = "};
  std::vector<QString> wStar2TextLabelsLabels = {"w**(int(0.2m))  = ", "w**(int(0.5m))  = ", "w**(int(0.8m))  = "};
  std::vector<QString> wStar3TextLabelsLabels = {"w***(int(0.2m)) = ", "w***(int(0.5m)) = ", "w***(int(0.8m)) = "};

  for(int i = 0; i < wTextLabelsLabels.size(); ++i){
    wStarTextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset, wStarTextLabelsLabels[i] + "0"));
    wStar2TextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset + 3 * verticalStep, wStar2TextLabelsLabels[i] + "0"));
    wStar3TextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset + 6 * verticalStep, wStar3TextLabelsLabels[i] + "0"));
    wTextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset + 9 * verticalStep, wTextLabelsLabels[i] + "0"));
    verticalOffset += verticalStep;
  }


  //====================  SECOND COLUMN =================//

  horizontalOffset = 0.20;
  verticalOffset = 0.01 + 9 * verticalStep;

  //==================== SUMMARIC ERRORS=================//

  horizontalOffset = 0.85;
  verticalOffset = 0.01;

  /*
  plotLabel expNumTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "Exp" + expNum);
   verticalOffset += verticalStep;
   verticalOffset += verticalStep;

  verticalOffset += verticalStep;

  plotLabel hTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "h  = " + formatNumberForDisplay(DESDAAlgorithm._h));
  verticalOffset += verticalStep;
  plotLabel hwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "hw = " + formatNumberForDisplay(DESDAAlgorithm._hWindowed));
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;
  */

  plotLabel L1WTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_w = 0");
  verticalOffset += verticalStep;

  plotLabel L1MTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_m = 0");
  verticalOffset += verticalStep;

  plotLabel L1DTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_d = 0");
  verticalOffset += verticalStep;

  plotLabel L1PTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_p = 0");
  verticalOffset += verticalStep;

  plotLabel L1NTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_n = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel L2WTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_w = 0");
  verticalOffset += verticalStep;

  plotLabel L2MTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_m = 0");
  verticalOffset += verticalStep;

  plotLabel L2DTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_d = 0");
  verticalOffset += verticalStep;

  plotLabel L2PTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_p = 0");
  verticalOffset += verticalStep;

  plotLabel L2NTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_n = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel supWTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_w = 0");
  verticalOffset += verticalStep;

  plotLabel supMTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_m = 0");
  verticalOffset += verticalStep;

  plotLabel supDTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_d = 0");
  verticalOffset += verticalStep;

  plotLabel supPTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_p = 0");
  verticalOffset += verticalStep;

  plotLabel supNTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_n = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel modWTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_w = 0");
  verticalOffset += verticalStep;

  plotLabel modMTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_m = 0");
  verticalOffset += verticalStep;

  plotLabel modDTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_d = 0");
  verticalOffset += verticalStep;

  plotLabel modPTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_p = 0");
  verticalOffset += verticalStep;

  plotLabel modNTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_n = 0");
  verticalOffset += verticalStep;

  fillDomain(&_domain, nullptr);
  for(auto pt : _domain) _drawableDomain.push_back(pt->at(0));

  ui->widget_plot->replot();
  qApp->processEvents();

  int numberOfErrorCalculations = 1;

  QVector<int> additionalScreensSteps = {};

  /*
  for(int i = 990; i < 1011; ++i){
      additionalScreensSteps.append(i);
  }
  */

  for(stepNumber = 1; stepNumber < stepsNumber; ++stepNumber)
  {
    clock_t executionStartTime = clock();

    DESDAAlgorithm.performStep();

    targetFunction.reset(generateTargetFunction(&means, &stDevs));

    if(stepNumber % screenGenerationFrequency == 0 || stepNumber < 10
       || additionalScreensSteps.contains(stepNumber))
    {
      qDebug() << "Drawing in step number " << stepNumber << ".";

      _kernelPrognosisDerivativeValues =
          DESDAAlgorithm.getKernelPrognosisDerivativeValues(&_drawableDomain);

      // Error calculations
      if(stepNumber >= 1000)
      {
        _windowedErrorDomain = DESDAAlgorithm.getWindowedErrorDomain();
        _errorDomain = DESDAAlgorithm.getErrorDomain();

        _windowedModelPlotY = getTargetFunctionValuesOnDomain(&_windowedErrorDomain);
        _windowedEstimatorErrorY = DESDAAlgorithm.getWindowKDEValues(&_windowedErrorDomain);
        _modelPlotErrorY = getTargetFunctionValuesOnDomain(&_errorDomain);
        _lessElementsEstimatorErrorY = DESDAAlgorithm.getKDEValues(&_errorDomain);
        _weightedEstimatorErrorY = DESDAAlgorithm.getWeightedKDEValues(&_errorDomain);
        _sigmoidallyEnhancedErrorPlotY = DESDAAlgorithm.getEnhancedKDEValues(&_errorDomain);
        _rareElementsEnhancedErrorPlotY = DESDAAlgorithm.getRareElementsEnhancedKDEValues(&_errorDomain);

        _L1_w += calculateL1Error(_windowedModelPlotY, _windowedEstimatorErrorY, _windowedErrorDomain);
        _L1_m += calculateL1Error(_modelPlotErrorY, _lessElementsEstimatorErrorY, _errorDomain);
        _L1_d += calculateL1Error(_modelPlotErrorY, _weightedEstimatorErrorY, _errorDomain);
        _L1_p += calculateL1Error(_modelPlotErrorY, _sigmoidallyEnhancedErrorPlotY, _errorDomain);
        _L1_n += calculateL1Error(_modelPlotErrorY, _rareElementsEnhancedErrorPlotY, _errorDomain);

        _L2_w += calculateL2Error(_windowedModelPlotY, _windowedEstimatorErrorY, _windowedErrorDomain);
        _L2_m += calculateL2Error(_modelPlotErrorY, _lessElementsEstimatorErrorY, _errorDomain);
        _L2_d += calculateL2Error(_modelPlotErrorY, _weightedEstimatorErrorY, _errorDomain);
        _L2_p += calculateL2Error(_modelPlotErrorY, _sigmoidallyEnhancedErrorPlotY, _errorDomain);
        _L2_n += calculateL2Error(_modelPlotErrorY, _rareElementsEnhancedErrorPlotY, _errorDomain);

        _sup_w += calculateSupError(_windowedModelPlotY, _windowedEstimatorErrorY);
        _sup_m += calculateSupError(_modelPlotErrorY, _lessElementsEstimatorErrorY);
        _sup_d += calculateSupError(_modelPlotErrorY, _weightedEstimatorErrorY);
        _sup_p += calculateSupError(_modelPlotErrorY, _sigmoidallyEnhancedErrorPlotY);
        _sup_n += calculateSupError(_modelPlotErrorY, _rareElementsEnhancedErrorPlotY);

        _mod_w += fabs(findExtrema(_windowedModelPlotY, _windowedErrorDomain) - findExtrema(_windowedEstimatorErrorY, _windowedErrorDomain));
        _mod_m += fabs(findExtrema(_modelPlotErrorY, _errorDomain) - findExtrema(_lessElementsEstimatorErrorY, _errorDomain));
        _mod_d += fabs(findExtrema(_modelPlotErrorY, _errorDomain) - findExtrema(_weightedEstimatorErrorY, _errorDomain));
        _mod_p += fabs(findExtrema(_modelPlotErrorY, _errorDomain) - findExtrema(_sigmoidallyEnhancedErrorPlotY, _errorDomain));
        _mod_n += fabs(findExtrema(_modelPlotErrorY, _errorDomain) - findExtrema(_rareElementsEnhancedErrorPlotY, _errorDomain));
        _mod_w += fabs(findExtrema(_windowedModelPlotY, _windowedErrorDomain) - findExtrema(_windowedEstimatorErrorY, _windowedErrorDomain));
        _mod_m += fabs(findExtrema(_modelPlotErrorY, _errorDomain) - findExtrema(_lessElementsEstimatorErrorY, _errorDomain));
        _mod_d += fabs(findExtrema(_modelPlotErrorY, _errorDomain) - findExtrema(_weightedEstimatorErrorY, _errorDomain));
        _mod_p += fabs(findExtrema(_modelPlotErrorY, _errorDomain) - findExtrema(_sigmoidallyEnhancedErrorPlotY, _errorDomain));
        _mod_n += fabs(findExtrema(_modelPlotErrorY, _errorDomain) - findExtrema(_rareElementsEnhancedErrorPlotY, _errorDomain));

        ++numberOfErrorCalculations;
      }

      // ============ SUMS =========== //

      L1WTextLabel
          .setText("L1_w  =" + formatNumberForDisplay(_L1_w / numberOfErrorCalculations));
      L1MTextLabel
          .setText("L1_m  =" + formatNumberForDisplay(_L1_m / numberOfErrorCalculations));
      L1DTextLabel
          .setText("L1_d  =" + formatNumberForDisplay(_L1_d / numberOfErrorCalculations));
      L1PTextLabel
          .setText("L1_p  =" + formatNumberForDisplay(_L1_p / numberOfErrorCalculations));
      L1NTextLabel
          .setText("L1_n  =" + formatNumberForDisplay(_L1_n / numberOfErrorCalculations));
      L2WTextLabel
          .setText("L2_w  =" + formatNumberForDisplay(_L2_w / numberOfErrorCalculations));
      L2MTextLabel
          .setText("L2_m  =" + formatNumberForDisplay(_L2_m / numberOfErrorCalculations));
      L2DTextLabel
          .setText("L2_d  =" + formatNumberForDisplay(_L2_d / numberOfErrorCalculations));
      L2PTextLabel
          .setText("L2_p  =" + formatNumberForDisplay(_L2_p / numberOfErrorCalculations));
      L2NTextLabel
          .setText("L2_n  =" + formatNumberForDisplay(_L2_n / numberOfErrorCalculations));
      supWTextLabel
          .setText("sup_w =" + formatNumberForDisplay(_sup_w / numberOfErrorCalculations));
      supMTextLabel
          .setText("sup_m =" + formatNumberForDisplay(_sup_m / numberOfErrorCalculations));
      supDTextLabel
          .setText("sup_d =" + formatNumberForDisplay(_sup_d / numberOfErrorCalculations));
      supPTextLabel
          .setText("sup_p =" + formatNumberForDisplay(_sup_p / numberOfErrorCalculations));
      supNTextLabel
          .setText("sup_n =" + formatNumberForDisplay(_sup_n / numberOfErrorCalculations));
      modWTextLabel
          .setText("mod_w =" + formatNumberForDisplay(_mod_w / numberOfErrorCalculations));
      modMTextLabel
          .setText("mod_m =" + formatNumberForDisplay(_mod_m / numberOfErrorCalculations));
      modDTextLabel
          .setText("mod_d =" + formatNumberForDisplay(_mod_d / numberOfErrorCalculations));
      modPTextLabel
          .setText("mod_p =" + formatNumberForDisplay(_mod_p / numberOfErrorCalculations));
      modNTextLabel
          .setText("mod_n =" + formatNumberForDisplay(_mod_n / numberOfErrorCalculations));


      // ============= LEFT SIDE UPDATE ================ //
      iTextLabel
          .setText(      "i     = " + QString::number(stepNumber));

      rareElementsTextLabel
          .setText(      "rare  = " + QString::number(_atypicalElementsValuesAndDerivatives.size()));

      mTextLabel.setText("m     = " + QString::number(DESDAAlgorithm._m));

      qTextLabel.setText("q     =" + formatNumberForDisplay(DESDAAlgorithm._quantileEstimator));

      KPSSTextLabel.setText(   "KPSS     = " + formatNumberForDisplay(
                                DESDAAlgorithm.getStationarityTestValue()));

      sgmKPSSTextLabel.setText( "sgmKPSS  = " + formatNumberForDisplay(
                                    DESDAAlgorithm._sgmKPSS));

      sgmKPSS2TextLabel.setText("sgmKPSS2 = " + formatNumberForDisplay(
                                    2 * DESDAAlgorithm._sgmKPSS - 1));

      rTextLabel.setText("r     =" + formatNumberForDisplay(DESDAAlgorithm._r));

      drawPlots(&DESDAAlgorithm);

      for(int i = 0; i < DESDAAlgorithm._examinedClustersW.size(); ++i){
        wTextLabels[i].setText(wTextLabelsLabels[i] + formatNumberForDisplay(
                                   DESDAAlgorithm._examinedClustersW[i]));
        wStarTextLabels[i].setText(wStarTextLabelsLabels[i] + formatNumberForDisplay(
                                   DESDAAlgorithm._examinedClustersWStar[i]));
        wStar2TextLabels[i].setText(wStar2TextLabelsLabels[i] + formatNumberForDisplay(
                                   DESDAAlgorithm._examinedClustersWStar2[i]));
        wStar3TextLabels[i].setText(wStar3TextLabelsLabels[i] + formatNumberForDisplay(
                                   DESDAAlgorithm._examinedClustersWStar3[i]));
      }

      for(int i = 0; i < DESDAAlgorithm._examinedClustersDerivatives.size(); ++i){
          gTextLabels[i].setText(gTextLabelsLabels[i] + formatNumberForDisplay(
                                  DESDAAlgorithm._examinedClustersDerivatives[i]));
      }

      maxAbsATextLabel.setText("max(|g|)     = " + formatNumberForDisplay(
                                   DESDAAlgorithm._maxAbsDerivativeValueInCurrentStep));


      vTextLabel.setText("v     =" + formatNumberForDisplay(DESDAAlgorithm._v));
      betaTextLabel.setText("beta0 = " + QString::number(DESDAAlgorithm._beta0));

       /*
      // DEBUG
      hTextLabel.setText( "h  = " + formatNumberForDisplay(DESDAAlgorithm._h));
      hwTextLabel.setText("hw = " + formatNumberForDisplay(DESDAAlgorithm._hWindowed));
      */

      ui->widget_plot->replot();
      qApp->processEvents();

      if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

      QString imageName = dirPath + QString::number(stepNumber) + ".png";
      qDebug() << "Image saved: " << ui->widget_plot->savePng(imageName, 0, 0, 1, -1);
    }
  }

  qDebug() << "Animation finished.";
}
