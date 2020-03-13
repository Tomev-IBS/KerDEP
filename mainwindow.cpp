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

    ModelValues.clear();

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

    // Generate less elements KDE plot (navy blue)
    if(ui->checkBox_showEstimationPlot->isChecked())
        addPlot(&_lessElementsEstimatorY, _KDE_PLOT_PEN);

    // Generate weighted estimator plot (light blue)
    if(ui->checkBox_showWeightedEstimationPlot->isChecked())
        addPlot(&_weightedEstimatorY, _WEIGHTED_PLOT_PEN);

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
      verticalLine->setPen(QPen(Qt::blue));
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
    qreal maxVal = 3 + ui->lineEdit_distributionProgression->text().toDouble() * 3000;  // Traveling case
    //qreal maxVal = 3 + ui->lineEdit_distributionProgression->text().toDouble(); /* should jump */

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
  int l = 4;

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
    &gt, newWeightB, mE, l, lambda
  );


  double horizontalOffset = 0.01, verticalOffset = 0.01, verticalStep = 0.03;

  plotLabel iTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "i    = 0");
  verticalOffset += verticalStep;

  plotLabel seedTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "seed = " + ui->lineEdit_seed->text());
  verticalOffset += verticalStep;

  plotLabel wTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "w    = 0.99");
  verticalOffset += verticalStep;

  plotLabel m0TextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "m0   = " + ui->lineEdit_sampleSize->text());
  verticalOffset += verticalStep;

  plotLabel mMinTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "mMin = " + QString::number(DESDAAlgorithm._minM));
  verticalOffset += verticalStep;

  plotLabel mTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "m    = " + QString::number(sampleSize));
  verticalOffset += verticalStep;

  plotLabel rTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "r    = " + QString::number(DESDAAlgorithm._r));
  verticalOffset += verticalStep;

  plotLabel qTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "q    = " + formatNumberForDisplay(DESDAAlgorithm._quantileEstimator));

  verticalOffset += verticalStep;
  plotLabel rareElementsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "rare = 0");
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
                       "max(|a|)      = 0");
  verticalOffset += verticalStep;

  std::vector<QString> aTextLabelsLabels = {"a(int(0.2m0)) = ", "a(int(0.5m0)) = ", "a(int(0.8m0)) = "};
  std::vector<plotLabel> aTextLabels = {};

  for(int i = 0; i < aTextLabelsLabels.size(); ++i){
    aTextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset, aTextLabelsLabels[i]));
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
    wStarTextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset, wStarTextLabelsLabels[i]));
    wStar2TextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset + 3 * verticalStep, wStar2TextLabelsLabels[i]));
    wStar3TextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset + 6 * verticalStep, wStar3TextLabelsLabels[i]));
    wTextLabels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset + 9 * verticalStep, wTextLabelsLabels[i]));
    verticalOffset += verticalStep;
  }


  //==================== SUMMARIC ERRORS=================//

  horizontalOffset = 0.72;
  verticalOffset = 0.01;

  verticalOffset += verticalStep;

  plotLabel error1SejwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_w = ");
  verticalOffset += verticalStep;

  plotLabel error1SejTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_m = ");
  verticalOffset += verticalStep;

  plotLabel error1SejsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_d = ");
  verticalOffset += verticalStep;

  plotLabel error1SejpTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_p = ");
  verticalOffset += verticalStep;

  plotLabel error1SejnTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_n = ");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel error2SejwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_w = ");
  verticalOffset += verticalStep;

  plotLabel error2SejTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_m = ");
  verticalOffset += verticalStep;

  plotLabel error2SejsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_d = ");
  verticalOffset += verticalStep;

  plotLabel error2SejpTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_p = ");
  verticalOffset += verticalStep;

  plotLabel error2SejnTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_n = ");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel errorSupSejwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_w = ");
  verticalOffset += verticalStep;

  plotLabel errorSupSejTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_m = ");
  verticalOffset += verticalStep;

  plotLabel errorSupSejsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_d = ");
  verticalOffset += verticalStep;

  plotLabel errorSupSejpTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_p = ");
  verticalOffset += verticalStep;

  plotLabel errorSupSejnTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_n = ");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel errorModSejwTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_n = ");
  verticalOffset += verticalStep;

  plotLabel errorModSejTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_m = ");
  verticalOffset += verticalStep;

  plotLabel errorModSejsTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_d = ");
  verticalOffset += verticalStep;

  plotLabel errorModSejpTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_p = ");
  verticalOffset += verticalStep;

  plotLabel errorModSejnTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_n = ");
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

      _windowedEstimatorY =
          DESDAAlgorithm.getWindowKDEValues(&_drawableDomain);

      _lessElementsEstimatorY =
          DESDAAlgorithm.getKDEValues(&_drawableDomain);

      _weightedEstimatorY =
          DESDAAlgorithm.getWeightedKDEValues(&_drawableDomain);

      _sigmoidallyEnhancedPlotY =
          DESDAAlgorithm.getEnhancedKDEValues(&_drawableDomain);

      _rareElementsEnhancedPlotY =
          DESDAAlgorithm.getRareElementsEnhancedKDEValues(&_drawableDomain);      

      _atypicalElementsValuesAndDerivatives =
                DESDAAlgorithm.getAtypicalElementsValuesAndDerivatives();

      _quantileEstimatorValue = DESDAAlgorithm._quantileEstimator;


      double maxKDEValue = 0.0;

      for(int i = 0; i < _lessElementsEstimatorY.size(); ++i)
      {
        double val = _lessElementsEstimatorY[i];

        if(val > maxKDEValue)
        {
          maxKDEValue = val;
          KDEExtrema = _drawableDomain[i];
        }
      }

      maxKDEValue = 0.0;

      for(int i = 0; i < _weightedEstimatorY.size(); ++i)
      {
        double val = _weightedEstimatorY[i];

        if(val > maxKDEValue)
        {
          maxKDEValue = val;
          WKDEExtrema = _drawableDomain[i];
        }
      }

      maxKDEValue = 0.0;

      for(int i = 0; i < _sigmoidallyEnhancedPlotY.size(); ++i)
      {
        double val = _sigmoidallyEnhancedPlotY[i];

        if(val > maxKDEValue)
        {
          maxKDEValue = val;
          KDEPExtrema = _drawableDomain[i];
        }
      }

      maxKDEValue = 0;

      for(int i = 0; i < _rareElementsEnhancedPlotY.size(); ++i)
      {
        double val = _rareElementsEnhancedPlotY[i];

        if(val > maxKDEValue)
        {
          maxKDEValue = val;
          REESEExtrema = _drawableDomain[i];
        }
      }

      maxKDEValue = 0;

      for(int i = 0; i < _windowedEstimatorY.size(); ++i)
      {
        double val = _windowedEstimatorY[i];

        if(val > maxKDEValue)
        {
          maxKDEValue = val;
          windowKDEExtrema = _drawableDomain[i];
        }
      }

      drawPlots(estimator.get(), targetFunction.get());

      QVector<qreal> errorHolder = {};

      // ser1_ej + sersup_ej
      double sup_ej = 0.0, sup_ejp = 0.0, sup_ejs = 0.0;

      for(int i = 0; i < ModelValues.size(); ++i){
        double val = fabs(ModelValues[i] - _lessElementsEstimatorY[i]);
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
        double val = fabs(ModelValues[i] - _weightedEstimatorY[i]);
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
        errorHolder.push_back(pow(ModelValues[i] - _lessElementsEstimatorY[i], 2));
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
        errorHolder.push_back(pow(ModelValues[i] - _weightedEstimatorY[i], 2));
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
        _summaricWindowKDEErrorSup += sup_ejw;
        _summaricKDEErrorMod  += mod_ej;
        _summaricKDEPErrorMod += mod_ejp;
        _summaricKDESErrorMod += mod_ejs;
        _summaricKDENErrorMod += mod_ejn;
        _summaricWindowKDEErrorMod += mod_ejw;
      }

      // ============ SUMS =========== //
      error1SejwTextLabel
          .setText("L1_w    =" + formatNumberForDisplay(_summaricWindowKDEError1));
      error1SejTextLabel
          .setText("L1_m    =" + formatNumberForDisplay(_summaricKDEError1));
      error1SejpTextLabel
          .setText("L1_p    =" + formatNumberForDisplay(_summaricKDEPError1));
      error1SejsTextLabel
          .setText("L1_d    =" + formatNumberForDisplay(_summaricKDESError1));
      error1SejnTextLabel
          .setText("L1_n    =" + formatNumberForDisplay(_summaricKDENError1));
      error2SejwTextLabel
          .setText("L2_w    =" + formatNumberForDisplay(_summaricWindowKDEError2));
      error2SejTextLabel
          .setText("L2_m    =" + formatNumberForDisplay(_summaricKDEError2));
      error2SejpTextLabel
          .setText("L2_p    =" + formatNumberForDisplay(_summaricKDEPError2));
      error2SejsTextLabel
          .setText("L2_d    =" + formatNumberForDisplay(_summaricKDESError2));
      error2SejnTextLabel
          .setText("L2_n    =" + formatNumberForDisplay(_summaricKDENError2));
      errorSupSejwTextLabel
          .setText("sup_w   =" + formatNumberForDisplay(_summaricWindowKDEErrorSup));
      errorSupSejTextLabel
          .setText("sup_m   =" + formatNumberForDisplay(_summaricKDEErrorSup));
      errorSupSejpTextLabel
          .setText("sup_p   =" + formatNumberForDisplay(_summaricKDEPErrorSup));
      errorSupSejsTextLabel
          .setText("sup_d   =" + formatNumberForDisplay(_summaricKDESErrorSup));
      errorSupSejnTextLabel
          .setText("sup_n   =" + formatNumberForDisplay(_summaricKDENErrorSup));
      errorModSejwTextLabel
          .setText("mod_w   =" + formatNumberForDisplay(_summaricWindowKDEErrorMod));
      errorModSejTextLabel
          .setText("mod_m   =" + formatNumberForDisplay(_summaricKDEErrorMod));
      errorModSejpTextLabel
          .setText("mod_p   =" + formatNumberForDisplay(_summaricKDEPErrorMod));
      errorModSejsTextLabel
          .setText("mod_d   =" + formatNumberForDisplay(_summaricKDESErrorMod));
      errorModSejnTextLabel
          .setText("mod_n   =" + formatNumberForDisplay(_summaricKDENErrorMod));


      // ============= LEFT SIDE UPDATE ================ //
      iTextLabel
          .setText(      "i    = " + QString::number(stepNumber));

      rareElementsTextLabel
          .setText(      "rare = " + QString::number(_atypicalElementsValuesAndDerivatives.size()));

      mTextLabel.setText("m    = " + QString::number(DESDAAlgorithm._m));

      qTextLabel.setText("q    =" + formatNumberForDisplay(DESDAAlgorithm._quantileEstimator));

      KPSSTextLabel.setText(   "KPSS     = " + formatNumberForDisplay(
                                DESDAAlgorithm.getStationarityTestValue()));

      sgmKPSSTextLabel.setText( "sgmKPSS  = " + formatNumberForDisplay(
                                    DESDAAlgorithm._sgmKPSS));

      sgmKPSS2TextLabel.setText("sgmKPSS2 = " + formatNumberForDisplay(
                                    2 * DESDAAlgorithm._sgmKPSS - 1));

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

      for(int i = 0; i < DESDAAlgorithm._examinedClustersAs.size(); ++i){
          aTextLabels[i].setText(aTextLabelsLabels[i] + formatNumberForDisplay(
                                  DESDAAlgorithm._examinedClustersAs[i]));
      }

      maxAbsATextLabel.setText("max(|a|)      = " + formatNumberForDisplay(
                                   DESDAAlgorithm.getMaxAbsAOnLastKPSSMSteps()));

      ui->widget_plot->replot();
      qApp->processEvents();

      QString googleDriveDir = "D:\\Dysk Google\\"; // Home

      QString dirPath = googleDriveDir + "TR Badania\\Eksperyment 462 ("
                        "v = " + ui->lineEdit_distributionProgression->text() +
                        ", r = " + QString::number(DESDAAlgorithm._r) +
                        ", lepsza wersja)\\";

      if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

      QString imageName = dirPath + QString::number(stepNumber) + ".png";
      qDebug() << "Image saved: " << ui->widget_plot->savePng(imageName, 0, 0, 1, -1);
    }
  }

  qDebug() << "Animation finished.";
}
