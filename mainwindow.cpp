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

void MainWindow::drawPlots(kernelDensityEstimator* estimator, function* targetFunction)
{
    clearPlot();
    resizePlot();

    QVector<qreal> X;
    QVector<qreal> modelDistributionY;

    // Fill domain with points
    // To keep things simple let's consider only these domains wherein
    // each dimension has equal size.

    ModelValues.clear();
    KDEValues.clear();

    double modelMax = 0.0;

    foreach(auto x, _domain)
    {
      double val = targetFunction->getValue(x.get());

      if(val > modelMax)
      {
        modelMax = val;
        modelExtrema = x->at(0);
      }

      modelDistributionY.append(val);
      ModelValues.push_back(val);
      X.append(x->at(0));
    }

    // Generate plot of model function
    if(ui->checkBox_showEstimatedPlot->isChecked())
      addPlot(&modelDistributionY, _MODEL_PLOT_PEN);

    // Generate a vector of values from selected KDE
    KDEEstimationY.clear();

    double val;
    _maxEstimatorValueOnDomain = 0;

    estimator->_shouldConsiderWeights = false;

    // TODO: Place counting in another thread
    for(int i = 0; i < _domain.size(); ++i)
    {
      auto x = _domain[i];

      val = estimator->getValue(x.get());

      if (val > _maxEstimatorValueOnDomain)
      {
        _maxEstimatorValueOnDomain = val;
        KDEExtrema = X[i];
      }

      oldKerernelY.append(val);

      KDEEstimationY.append(val);

      KDEValues.push_back(val);
    }

    estimator->_shouldConsiderWeights = true;

    // Generate a plot of KDE
    if(ui->checkBox_showEstimationPlot->isChecked())
      addPlot(&KDEEstimationY, _KDE_PLOT_PEN);

    double maxWKDE = 0;

    if(ui->checkBox_showWeightedEstimationPlot->isChecked())
    {
      WKDEValues.clear();

      for(int i = 0; i < _domain.size(); ++i)
      {
        auto x = _domain[i];

        val = estimator->getValue(x.get());

        WKDEValues.push_back(val);

        if(val > maxWKDE)
        {
          maxWKDE = val;
          WKDEExtrema = X[i];
        }
      }

      estimator->_shouldConsiderWeights = false;

      //addWeightedEstimatorPlot(&X, &WKDEValues);
      addPlot(&WKDEValues, _WEIGHTED_PLOT_PEN);
    }

    // Generate full estomator plot (BLACK)

    if(ui->checbox_showFullEstimator->isChecked())
      addPlot(&_windowedEstimatorY, _WINDOWED_PLOT_PEN);

    // Generate plot for kernel prognosis derivative
    if(ui->checkBox_kernelPrognosedPlot->isChecked())
      addPlot(&_kernelPrognosisDerivativeValues, _DERIVATIVE_PLOT_PEN);

    if(ui->checkBox_sigmoidallyEnhancedKDE->isChecked())
      addPlot(&_sigmoidallyEnhancedPlotY, _DESDA_KDE_PLOT_PEN);

    if(ui->checkBox_showUnusualClusters->isChecked())
      markUncommonClusters();

    if(ui->checkBox_REESEKDE->isChecked())
      addPlot(&_rareElementsEnhancedPlotY, _DESDA_RARE_ELEMENTS_KDE_PLOT_PEN);
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
            //maxX = ui->lineEdit_maxX->text().toDouble(), // Standard
            maxX = 3 + ui->lineEdit_distributionProgression->text().toDouble() * 3000, // For progression
            //maxX = 3 + ui->lineEdit_distributionProgression->text().toDouble(), // For jump
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
  for(auto x : _atypicalElementsValues){
    // Only works for distribution data
    QCPItemLine *verticalLine = new QCPItemLine(ui->widget_plot);
    verticalLine->start->setCoords(x, 0);
    verticalLine->end->setCoords(x, -_quantileEstimatorValue);
    verticalLine->setPen(QPen(Qt::blue));
    _linesOnPlot.push_back(verticalLine);
  }

  return _atypicalElementsValues.size();
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
    qreal maxVal = 3 + ui->lineEdit_distributionProgression->text().toDouble() * 3000;  // Traveling case
    //qreal maxVal = 3 + ui->lineEdit_distributionProgression->text().toDouble(); // Jump case

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

    double maxMean = ui->lineEdit_distributionProgression->text().toDouble() *  3000;

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        contributions.push_back
        (
            (static_cast<QLineEdit*>(ui->tableWidget_targetFunctions->cellWidget(functionIndex, CONTRIBUTION_COLUMN_INDEX)))
            ->text().toDouble()
        );

        elementalDistributions.push_back(std::shared_ptr<distribution>(new normalDistribution(seed, (*means)[functionIndex].get(), (*stDevs)[functionIndex].get(), maxMean)));
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

  // Log that application started generating KDE
  qDebug() << "KDE animation started.";
  qDebug() << "Seed: " + ui->lineEdit_seed->text() +
              ", Sample size: " + ui->lineEdit_sampleSize->text();

  srand(static_cast<unsigned int>(ui->lineEdit_seed->text().toInt()));

  fillMeans(&means);
  fillStandardDeviations(&stDevs);

  std::shared_ptr<function>
    targetFunction(generateTargetFunction(&means, &stDevs));

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
                                          4000 /* delay */,
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
  int l = 8;

  groupingThread gt(&storedMedoids, parser);

  gt.setAttributesData(&attributesData);

  qDebug() << "Attributes data set.";

  int sampleSize = ui->lineEdit_sampleSize->text().toInt();
  gt.initialize(medoidsNumber, sampleSize);

  _a = MIN_A;

  _longestStepExecutionInSecs = 0;

  double newWeightB = 0.5;
  int mE = ui->lineEdit_sampleSize->text().toInt() / 2;

  clusters = &storedMedoids;

  weightedSilvermanSmoothingParameterCounter smoothingParamCounter(clusters, 0);

  kernelPrognoser->_shouldConsiderWeights = false;

  int lambda = 500;

  DESDA DESDAAlgorithm(
    estimator,
    kernelPrognoser,
    _enchancedKDE,
    ui->lineEdit_weightModifier->text().toDouble(),
    &smoothingParamCounter,
    algorithm,
    clusters,
    &storedMedoids,
    ui->lineEdit_rarity->text().toDouble(),
    &gt,
    ui->lineEdit_distributionProgression->text().toDouble(),
    newWeightB, mE, l, lambda
  );


  double horizontalOffset = 0.01, verticalOffset = 0.01, verticalStep = 0.03;

  plotLabel iTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "i   = ");
  verticalOffset += verticalStep;

  plotLabel mTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "m   = " + QString::number(sampleSize));
   verticalOffset += verticalStep;

  plotLabel EmETextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "EmE = ");
  verticalOffset += verticalStep;

  plotLabel ESmETextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "a_EmExK = ");
  verticalOffset += verticalStep;

  plotLabel maxATextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "avg(max(a...)) = ");
  verticalOffset += verticalStep;

  plotLabel stationarityTestTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "eta = ");
  verticalOffset += verticalStep;

  plotLabel stationarityPKTestTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "eta_PK = ");
  verticalOffset += verticalStep;

  plotLabel uTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "u = ");
  verticalOffset += verticalStep;

  plotLabel bTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "b = " + QString::number(newWeightB));
  verticalOffset += verticalStep;

  plotLabel rTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "r = " + QString::number(DESDAAlgorithm._r));
  verticalOffset += verticalStep;

  plotLabel qTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "q = " + formatNumberForDisplay(DESDAAlgorithm._quantileEstimator));

  //==================== SUMMARIC ERRORS=================//

  horizontalOffset = 0.72;
  verticalOffset = 0.01;
  verticalStep = 0.03;

  plotLabel m0minTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "m0   = " + ui->lineEdit_sampleSize->text() + ", m_min = " + QString::number(DESDAAlgorithm._minM));
  verticalOffset += verticalStep;

  plotLabel vTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "v = " + ui->lineEdit_distributionProgression->text());
  verticalOffset += verticalStep;

  plotLabel mETextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mE  = " + QString::number(mE));
  verticalOffset += verticalStep;

  plotLabel wTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "w_a = 0.99, w_EmE = " + QString::number(DESDAAlgorithm.w_E));
  verticalOffset += verticalStep;

  plotLabel alfaBetaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "alpha = " + QString::number(DESDAAlgorithm._alpha) + ", beta = " + QString::number(DESDAAlgorithm._beta));
  verticalOffset += verticalStep;

  plotLabel deltaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "delta = " + QString::number(DESDAAlgorithm.delta));
  verticalOffset += verticalStep;

  plotLabel lambdaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "lambda = " + QString::number(lambda));
  verticalOffset += verticalStep;

  plotLabel rareElementsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "#Rare elements = 0");
  verticalOffset += verticalStep;

  verticalOffset += verticalStep;

  plotLabel error1SejwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser1_ejw = ");
  verticalOffset += verticalStep;

  plotLabel error1SejTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser1_ej = ");
  verticalOffset += verticalStep;

  plotLabel error1SejsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser1_ejs = ");
  verticalOffset += verticalStep;

  plotLabel error1SejpTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser1_ejp = ");
  verticalOffset += verticalStep;

  plotLabel error1SejnTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser1_ejn = ");
  verticalOffset += verticalStep;

  plotLabel error2SejwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser2_ejw = ");
  verticalOffset += verticalStep;

  plotLabel error2SejTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser2_ej = ");
  verticalOffset += verticalStep;

  plotLabel error2SejsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser2_ejs = ");
  verticalOffset += verticalStep;

  plotLabel error2SejpTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser2_ejp = ");
  verticalOffset += verticalStep;

  plotLabel error2SejnTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "ser2_ejn = ");
  verticalOffset += verticalStep;

  plotLabel errorSupSejwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sersup_ejw = ");
  verticalOffset += verticalStep;

  plotLabel errorSupSejTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sersup_ej = ");
  verticalOffset += verticalStep;

  plotLabel errorSupSejsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sersup_ejs = ");
  verticalOffset += verticalStep;

  plotLabel errorSupSejpTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sersup_ejp = ");
  verticalOffset += verticalStep;

  plotLabel errorSupSejnTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sersup_ejn = ");
  verticalOffset += verticalStep;

  plotLabel errorModSejwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sermod_ejw = ");
  verticalOffset += verticalStep;

  plotLabel errorModSejTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sermod_ej = ");
  verticalOffset += verticalStep;

  plotLabel errorModSejsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sermod_ejs = ");
  verticalOffset += verticalStep;

  plotLabel errorModSejpTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sermod_ejp = ");
  verticalOffset += verticalStep;

  plotLabel errorModSejnTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sermod_ejn = ");
  verticalOffset += verticalStep;

  fillDomain(&_domain, nullptr);
  for(auto pt : _domain) _drawableDomain.push_back(pt->at(0));

  for(stepNumber = 1; stepNumber < stepsNumber; ++stepNumber)
  {
    clock_t executionStartTime = clock();

    DESDAAlgorithm.performStep();

    targetFunction.reset(generateTargetFunction(&means, &stDevs));

    if(stepNumber % 10 == 0)
    {
      qDebug() << "Drawing in step number " << stepNumber << ".";

      _kernelPrognosisDerivativeValues =
          DESDAAlgorithm.getKernelPrognosisDerivativeValues(&_drawableDomain);

      double maxA = 0, meanA = 0, currentMaxA = 0, avgMaxA = 0;

      for(auto val : _kernelPrognosisDerivativeValues){
        meanA += fabs(val);
        currentMaxA = currentMaxA < fabs(val) ? fabs(val) : currentMaxA;
      }

      if(stepNumber > sampleSize)
        maxAs.pop_front();

      maxAs.push_back(currentMaxA);

      for(auto val : maxAs)
      {
        avgMaxA += val;
        maxA = val > maxA ? val : maxA;
      }

      avgMaxA /= maxAs.size();

      // TR TODO: Notice that this should be placed inside DESDA algorithm
      DESDAAlgorithm.gamma = 1.0 / avgMaxA;

      _sigmoidallyEnhancedPlotY =
          DESDAAlgorithm.getEnhancedKDEValues(&_drawableDomain);

      _rareElementsEnhancedPlotY =
          DESDAAlgorithm.getRareElementsEnhancedKDEValues(&_drawableDomain);

      _windowedEstimatorY =
          DESDAAlgorithm.getWindowKDEValues(&_drawableDomain);

      _atypicalElementsValues =
                DESDAAlgorithm.getAtypicalElementsValues();

      _quantileEstimatorValue = DESDAAlgorithm._quantileEstimator;

      double maxKDEP = 0.0;

      for(int i = 0; i < _sigmoidallyEnhancedPlotY.size(); ++i)
      {
        double val = _sigmoidallyEnhancedPlotY[i];

        if(val > maxKDEP)
        {
          maxKDEP = val;
          KDEPExtrema = _drawableDomain[i];
        }
      }

      maxKDEP = 0;

      for(int i = 0; i < _rareElementsEnhancedPlotY.size(); ++i)
      {
        double val = _rareElementsEnhancedPlotY[i];

        if(val > maxKDEP)
        {
          maxKDEP = val;
          REESEExtrema = _drawableDomain[i];
        }
      }

      windowKDEExtrema = 0;

      for(int i = 0; i < _windowedEstimatorY.size(); ++i)
      {
        double val = _rareElementsEnhancedPlotY[i];

        if(val > maxKDEP)
        {
          maxKDEP = val;
          REESEExtrema = _drawableDomain[i];
        }
      }

      drawPlots(estimator.get(), targetFunction.get());

      cluster emE = DESDAAlgorithm.getEmECluster();

      // E1000
      double avg_mE = emE._currentKDEValue;

      // a_E1000xK
      double EmEEst = emE.predictionParameters[1];

      QVector<qreal> errorHolder = {};

      // ser1_ej + sersup_ej
      double sup_ej = 0.0, sup_ejp = 0.0, sup_ejs = 0.0;

      for(int i = 0; i < ModelValues.size(); ++i){
        double val = fabs(ModelValues[i] - KDEValues[i]);
        errorHolder.push_back(val);
        sup_ej = val > sup_ej ? val : sup_ej;
      }

      _errorEJ = numericIntegral(&errorHolder);
      errorHolder.clear();

      // ser1_ejp + sersup_ejp
      for(int i = 0; i < ModelValues.size(); ++i){
        double val = fabs(ModelValues[i] - _sigmoidallyEnhancedPlotY[i]);
        errorHolder.push_back(val);
        sup_ejp = val > sup_ejp ? val : sup_ejp;
      }

      _errorEJP = numericIntegral(&errorHolder);
      errorHolder.clear();

      // ser_ejs + sersup_ejs
      for(int i = 0; i < ModelValues.size(); ++i){
        double val = fabs(ModelValues[i] - WKDEValues[i]);
        errorHolder.push_back(val);
        sup_ejs = val > sup_ejs ? val : sup_ejs;
      }

      double error_ejs = 0;
      error_ejs = numericIntegral(&errorHolder);
      errorHolder.clear();

      // ser_1ejn + sersup_ejn
      double errorEJN = 0, sup_ejn = 0;
      for(int i = 0; i < ModelValues.size(); ++i){
        double val = fabs(ModelValues[i] - _rareElementsEnhancedPlotY[i]);
        errorHolder.push_back(val);
        sup_ejn = val > sup_ejn ? val : sup_ejn;
      }

      errorEJN = numericIntegral(&errorHolder);
      errorHolder.clear();

      // ser_1ejw + sersup_ejw
      double errorEJW = 0, sup_ejw = 0;
      for(int i = 0; i < ModelValues.size(); ++i){
        double val = fabs(ModelValues[i] - _windowedEstimatorY[i]);
        errorHolder.push_back(val);
        sup_ejw = val > sup_ejw ? val : sup_ejw;
      }

      errorEJW = numericIntegral(&errorHolder);
      errorHolder.clear();

      double error2EJ = 0.0, error2EJP = 0.0, error2EJS = 0.0, error2EJN = 0, error2EJW = 0;

      // ser2_ej
      for(int i = 0; i < ModelValues.size(); ++i){
        errorHolder.push_back(pow(ModelValues[i] - KDEValues[i], 2));
      }

      error2EJ = numericIntegral(&errorHolder);
      errorHolder.clear();

      // ser2_ejp
      for(int i = 0; i < ModelValues.size(); ++i){
        errorHolder.push_back(pow(ModelValues[i] - _sigmoidallyEnhancedPlotY[i], 2));
      }

      error2EJP = numericIntegral(&errorHolder);
      errorHolder.clear();

      // ser2_ejp
      for(int i = 0; i < ModelValues.size(); ++i){
        errorHolder.push_back(pow(ModelValues[i] - WKDEValues[i], 2));
      }

      error2EJS = numericIntegral(&errorHolder);
      errorHolder.clear();

      // ser2_ejn
      for(int i = 0; i < ModelValues.size(); ++i){
        errorHolder.push_back(pow(ModelValues[i] - _rareElementsEnhancedPlotY[i], 2));
      }

      error2EJN = numericIntegral(&errorHolder);
      errorHolder.clear();

      // ser2_ejw
      for(int i = 0; i < ModelValues.size(); ++i){
        errorHolder.push_back(pow(ModelValues[i] - _windowedEstimatorY[i], 2));
      }

      error2EJW = numericIntegral(&errorHolder);
      errorHolder.clear();

      double mod_ejw = 0, mod_ej = 0, mod_ejp = 0, mod_ejs = 0, mod_ejn = 0;
      mod_ej = fabs(modelExtrema - KDEExtrema);
      mod_ejp = fabs(modelExtrema - KDEPExtrema);
      mod_ejs = fabs(modelExtrema - WKDEExtrema);
      mod_ejn = fabs(modelExtrema - REESEExtrema);
      mod_ejw = fabs(modelExtrema - windowKDEExtrema);

      if(stepNumber >= 2000)
      {
        _summaricKDEError1    += _errorEJ;
        _summaricKDEPError1   += _errorEJP;
        _summaricKDESError1   += error_ejs;
        _summaricKDENError1   += errorEJN;
        _summaricWindowKDEError1 += errorEJW;
        _summaricKDEError2    += error2EJ;
        _summaricKDEPError2   += error2EJP;
        _summaricKDESError2   += error2EJS;
        _summaricKDENError2   += error2EJN;
        _summaricWindowKDEError2 += error2EJW;
        _summaricKDEErrorSup  += sup_ej;
        _summaricKDEPErrorSup += sup_ejp;
        _summaricKDESErrorSup += sup_ejs;
        _summaricKDENErrorSup += sup_ejn;
        _summaricWindowKDEErrorSup = sup_ejw;
        _summaricKDEErrorMod  += mod_ej;
        _summaricKDEPErrorMod += mod_ejp;
        _summaricKDESErrorMod += mod_ejs;
        _summaricKDENErrorMod += mod_ejn;
        _summaricWindowKDEErrorMod += mod_ejw;

      }

      EmETextLabel
          .setText("EmE            = " + formatNumberForDisplay(avg_mE));

      ESmETextLabel
          .setText("a_EmE xK       = " + formatNumberForDisplay(EmEEst * 1000));

      maxATextLabel
          .setText("EmE(max(|ai|)) = " + formatNumberForDisplay(avgMaxA));

      uTextLabel
          .setText("u     = " + formatNumberForDisplay(DESDAAlgorithm._u_i));

      // ============ SUMS =========== //
      iTextLabel
          .setText("i   = " + QString::number(stepNumber));
      error1SejwTextLabel
          .setText("ser1_ejw    = " + formatNumberForDisplay(_summaricWindowKDEError1));
      error1SejTextLabel
          .setText("ser1_ej     = " + formatNumberForDisplay(_summaricKDEError1));
      error1SejpTextLabel
          .setText("ser1_ejp    = " + formatNumberForDisplay(_summaricKDEPError1));
      error1SejsTextLabel
          .setText("ser1_ejs    = " + formatNumberForDisplay(_summaricKDESError1));
      error1SejnTextLabel
          .setText("ser1_ejn    = " + formatNumberForDisplay(_summaricKDENError1));
      error2SejwTextLabel
          .setText("ser2_ejw    = " + formatNumberForDisplay(_summaricWindowKDEError2));
      error2SejTextLabel
          .setText("ser2_ej     = " + formatNumberForDisplay(_summaricKDEError2));
      error2SejpTextLabel
          .setText("ser2_ejp    = " + formatNumberForDisplay(_summaricKDEPError2));
      error2SejsTextLabel
          .setText("ser2_ejs    = " + formatNumberForDisplay(_summaricKDESError2));
      error2SejnTextLabel
          .setText("ser2_ejn    = " + formatNumberForDisplay(_summaricKDENError2));
      errorSupSejwTextLabel
          .setText("sersup_ejw  = " + formatNumberForDisplay(_summaricWindowKDEErrorSup));
      errorSupSejTextLabel
          .setText("sersup_ej   = " + formatNumberForDisplay(_summaricKDEErrorSup));
      errorSupSejpTextLabel
          .setText("sersup_ejp  = " + formatNumberForDisplay(_summaricKDEPErrorSup));
      errorSupSejsTextLabel
          .setText("sersup_ejs  = " + formatNumberForDisplay(_summaricKDESErrorSup));
      errorSupSejnTextLabel
          .setText("sersup_ejn  = " + formatNumberForDisplay(_summaricKDENErrorSup));
      errorModSejwTextLabel
          .setText("sermod_ejw  = " + formatNumberForDisplay(_summaricWindowKDEErrorMod));
      errorModSejTextLabel
          .setText("sermod_ej   = " + formatNumberForDisplay(_summaricKDEErrorMod));
      errorModSejpTextLabel
          .setText("sermod_ejp  = " + formatNumberForDisplay(_summaricKDEPErrorMod));
      errorModSejsTextLabel
          .setText("sermod_ejs  = " + formatNumberForDisplay(_summaricKDESErrorMod));
      errorModSejnTextLabel
          .setText("sermod_ejn  = " + formatNumberForDisplay(_summaricKDENErrorMod));

      rareElementsTextLabel
          .setText("#Rare elements = " + QString::number(_atypicalElementsValues.size()));

      deltaTextLabel.setText("delta = " + formatNumberForDisplay(DESDAAlgorithm.delta));

      stationarityTestTextLabel
          .setText("eta = " + formatNumberForDisplay(DESDAAlgorithm.getStationarityTestValue()));

      stationarityPKTestTextLabel
          .setText("eta_PK = " + formatNumberForDisplay(DESDAAlgorithm.getPKStationarityTestValue()));

      bTextLabel
          .setText("b = " + QString::number(DESDAAlgorithm._newWeightB));

      mTextLabel.setText("m   = " + QString::number(DESDAAlgorithm._m));

      qTextLabel.setText("q = " + formatNumberForDisplay(DESDAAlgorithm._quantileEstimator));

      mETextLabel.setText("mE  = " + QString::number(DESDAAlgorithm._mE)
                           + ", mEta = " + QString::number(DESDAAlgorithm._kpssM));


      ui->widget_plot->replot();
      qApp->processEvents();

      QString googleDriveDir = "D:\\Dysk Google\\"; // Home

      QString dirPath = googleDriveDir + "TR Badania\\reEksperyment 449.4 ("
                        "v = " + ui->lineEdit_distributionProgression->text() +
                        "r = " + QString::number(DESDAAlgorithm._r) +
                        ")\\";

      if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

      QString imageName = dirPath + QString::number(stepNumber) + ".png";
      qDebug() << "Image saved: " << ui->widget_plot->savePng(imageName, 0, 0, 1, -1);
    }

    clock_t executionFinishTime = clock();
    double stepExecutionTime =
      static_cast<double>(executionFinishTime - executionStartTime);
    stepExecutionTime /= CLOCKS_PER_SEC;

    if(stepExecutionTime > _longestStepExecutionInSecs)
      _longestStepExecutionInSecs = stepExecutionTime;

    //qDebug() << "Longest execution time: " << _longestStepExecutionInSecs;
    //qDebug() << "Current execution time: " << stepExecutionTime;

    /*
    delay(static_cast<int>(
      (_longestStepExecutionInSecs - stepExecutionTime) * 1000
    ));
    */

  }

  qDebug() << "Animation finished.";
  // Value print for exp300+
  qDebug() << formatNumberForDisplay(_summaricKDEPError1);
  qDebug() << formatNumberForDisplay(_summaricKDEPError2);
  qDebug() << formatNumberForDisplay(_summaricKDEPErrorSup);
  qDebug() << formatNumberForDisplay(_summaricKDEPErrorMod);
}
