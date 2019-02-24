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

void MainWindow::keyPressEvent(QKeyEvent *event)
{
  unsigned int keyCode = static_cast<unsigned int>(event->key());

  //qDebug() << keyCode;

  // If key is a number
  if(keyCode >= 48 && keyCode <= 57)
  {
    unsigned int objectsNumber = keyCode - 48;

    if(objectsNumber == 0) objectsNumber += 10;

    insertObjectsBetweenIntervals(objectsNumber);
  }

  // If m is pressed
  if(keyCode == 77) insertMassiveData();
}

unsigned long long MainWindow::insertObjectsBetweenIntervals(
    unsigned int objectsNumber)
{
  std::vector<std::shared_ptr<sample>> interIntervalObjects;

  generateInterIntervalObjects(&interIntervalObjects, objectsNumber);
  selectDesiredNumberOfInterIntervalObjects(&interIntervalObjects);
  insertClustersFromInterIntervalObjects(&interIntervalObjects);

  objects.insert(objects.end(), interIntervalObjects.begin(), interIntervalObjects.end());

  qDebug() << "Inter interval objects inserterd: " << interIntervalObjects.size();

  return interIntervalObjects.size();
}

unsigned long long MainWindow::generateInterIntervalObjects(
    std::vector<std::shared_ptr<sample> > *interIntervalObjects,
    unsigned int objectsNumber)
{
  //TR TODO: Check if reader and parser are initialized

  for(unsigned int i = 0; i < objectsNumber; ++i)
  {
    reader->getNextRawDatum(parser->buffer);

    parser->addDatumToContainer(interIntervalObjects);

    parser->writeDatumOnPosition(
      interIntervalObjects,
      static_cast<int>(interIntervalObjects->size()) - 1
    );
  }

  return interIntervalObjects->size();
}

unsigned long long MainWindow::selectDesiredNumberOfInterIntervalObjects(
    std::vector<std::shared_ptr<sample> > *interIntervalObjects)
{
  int desiredNumberOfClusters =
      ui->lineEdit_interIntervalClusters->text().toInt();

  while(interIntervalObjects->size()
        > static_cast<unsigned int>(desiredNumberOfClusters))
    interIntervalObjects->erase(interIntervalObjects->begin() + static_cast<int>(static_cast<unsigned int>(rand()) % interIntervalObjects->size()));

  return interIntervalObjects->size();
}

unsigned long long MainWindow::insertClustersFromInterIntervalObjects(
    std::vector<std::shared_ptr<sample> > *interIntervalObjects)
{
  std::vector<std::shared_ptr<cluster>> newClusters;

  for(unsigned int i = 0; i < interIntervalObjects->size(); ++i)
  {
    newClusters.push_back(
      std::shared_ptr<cluster>(
        new cluster(static_cast<long>(clusters->size() + i),
                    (*interIntervalObjects)[i])
      )
    );
    newClusters.back()->setTimestamp(stepNumber);
  }

  setInterIntervalClustersWeights(&newClusters);

  clusters->insert(clusters->end(), newClusters.begin(), newClusters.end());

  return newClusters.size();
}

double MainWindow::setInterIntervalClustersWeights(
    std::vector<std::shared_ptr<cluster> > *newClusters)
{
  double weight = countInterIntervalClustersWeight();

  for(unsigned int i = 0; i < newClusters->size(); ++i)
    (*newClusters)[i]->setWeight(weight);

  return weight;
}

double MainWindow::countInterIntervalClustersWeight()
{
  // All values in milliseconds.
  double intervalValue = ui->lineEdit_milisecondsDelay->text().toInt();

  if(intervalValue < 1e-10) return 1.0;

  long long end = std::chrono::duration_cast< std::chrono::milliseconds >(
        std::chrono::system_clock::now().time_since_epoch()).count();

  long long difference = end - start;

  double weightModifier = ui->lineEdit_weightModifier->text().toDouble();

  double power = 1.0 -(difference/intervalValue);


  // TR TODO: This is neccessary for proper functionality when intervals are low.

  if(power < 0)
  {
    qDebug() << "Power is below 0.";
    return 0.0;
  }


  return pow(weightModifier, power);
}

unsigned long long MainWindow::insertMassiveData()
{
  std::vector<std::shared_ptr<sample>> massiveData;

  qDebug() << "Generating data.";

  generateMassiveData(&massiveData);

  qDebug() << "Massive data generated.";

  clusterMassiveData(&massiveData, &storedMedoids);

  qDebug() << "Massive data clustered.";

  return massiveData.size();
}

unsigned long long MainWindow::generateMassiveData(
  std::vector<std::shared_ptr<sample>> *dataContainer)
{
  long dataSize = 10000;

  dataContainer->clear();

  for(int i = 0; i < dataSize; ++i)
  {
    reader->getNextRawDatum(parser->buffer);

    parser->addDatumToContainer(dataContainer);

    parser->writeDatumOnPosition(dataContainer,
                                 static_cast<int>(dataContainer->size())-1);
  }

  return dataContainer->size();
}

void MainWindow::clusterMassiveData(
  std::vector<std::shared_ptr<sample>> *objects,
  std::vector<std::vector<std::shared_ptr<cluster>>> *storage)
{
  // Select medoids
  std::set<int> medoidsIndexes;

  // TODO TR: Add ui control
  unsigned int medoidsNumber = 10;

  do
  {
    medoidsIndexes.insert(
      static_cast<int>(static_cast<unsigned int>(rand()) % objects->size())
    );
  } while(medoidsIndexes.size() < medoidsNumber);

  // Create clusters from medoids
  std::set<int>::iterator it = medoidsIndexes.begin();

  if(storage->size() == 0) storage->push_back(std::vector<std::shared_ptr<cluster>>());

  for(unsigned int i = 0; i < medoidsNumber; ++i)
  {
    storage->at(0).push_back(std::make_shared<cluster>(cluster(static_cast<long>(i), objects->at(static_cast<unsigned int>(*it)))));
    storage->at(0).back()->setTimestamp(stepNumber);
    storage->at(0).back()->setWeight(objects->size() / medoidsNumber);
    std::advance(it, 1);
  }
}

std::vector<std::shared_ptr<cluster>> MainWindow::getClustersForEstimator()
{
  std::vector<std::shared_ptr<cluster>> consideredClusters;

  for(std::vector<std::shared_ptr<cluster>> level : storedMedoids)
  {
    for(std::shared_ptr<cluster> c : level)
    {
      if(c->getWeight() >= positionalSecondGradeEstimator)
        consideredClusters.push_back(c);
    }
  }

  return consideredClusters;
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
    ui->widget_plot->xAxis->setLabel("x");
    ui->widget_plot->yAxis->setLabel("f(x)");

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

void MainWindow::drawPlots(kernelDensityEstimator* estimator,
                           function* targetFunction)
{
    clearPlot();

    resizePlot();

    QVector<qreal> X;
    QVector<qreal> modelDistributionY;

    // Fill domain with points
    // To keep things simple let's consider only these domains wherein
    // each dimension has equal size.

    QVector<std::shared_ptr<point>> domain;
    fillDomain(&domain, nullptr);
    ModelValues.clear();
    KDEValues.clear();

    double modelMax = 0.0;

    foreach(auto x, domain)
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
      addModelPlot(&X, &modelDistributionY);

    // Generate a vector of values from selected KDE
    KDEEstimationY.clear();

    double val;
    _maxEstimatorValueOnDomain = 0;

    // TODO: Place counting in another thread
    //foreach(std::shared_ptr<point> x, domain)
    for(int i = 0; i < domain.size(); ++i)
    {
      auto x = domain[i];

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

    qDebug() << "Max est val on domain: " << _maxEstimatorValueOnDomain;

    // Generate a plot of KDE
    if(ui->checkBox_showEstimatedPlot->isChecked())
      addEstimatedPlot(&X, &KDEEstimationY);

    // Generate a plot of temporal derivative
    KDETemporalDerivativeY.clear();
    double visibilityEnchantCoefficient = 3;
    double derivativeYOffset = 0.0;

    if(oldKerernelY.size() != 0)
    {
      for(int i = 0; i < KDEEstimationY.size(); ++i)
        KDETemporalDerivativeY.push_back(visibilityEnchantCoefficient*
                                         (KDEEstimationY[i] - oldKerernelY[i]) - derivativeYOffset);
    }

    if(ui->checkBox_showTimeDerivativePlot->isChecked())
      addTemporalDerivativePlot(&X, &KDETemporalDerivativeY);

    // Generate plot for kernel prognosis derivative
    if(ui->checkBox_kernelPrognosedPlot->isChecked())
    {
      //countKernelPrognosisDerivativeY(&X);
      addKernelPrognosisDerivativePlot(&X);
    }

    if(ui->checkBox_sigmoidallyEnhancedKDE->isChecked())
      addSigmoidallyEnhancedEstimationPlot(&X, estimator);

    if(ui->checkBox_showUnusualClusters->isChecked())
      markUncommonClusters();

    if(ui->checkBox_showNewTrends->isChecked())
      markNewTrends();

    if(ui->checkBox_negativeC2Clusters->isChecked())
      markClustersWithNegativeDerivative();

    // Draw plots
    ui->widget_plot->replot();

    oldKerernelY = KDEEstimationY;
}

void MainWindow::resizePlot()
{
    // Resize plot
    qreal   minX = ui->lineEdit_minX->text().toDouble(),
            maxX = ui->lineEdit_maxX->text().toDouble(),
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
}

void MainWindow::clearPlot()
{
    while(ui->widget_plot->graphCount() != 0)
        ui->widget_plot->removeGraph(0);
}

void MainWindow::addModelPlot(const QVector<qreal> *X, const QVector<qreal> *Y)
{
    //qDebug() << numericIntegral(Y);

    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, *Y);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)
        ->setPen(QPen(Qt::red));
}

void MainWindow::addEstimatedPlot(const QVector<qreal> *X, const QVector<qreal> *Y)
{
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, *Y);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)
        ->setPen(QPen(Qt::blue));

    // Derivative + KDE
    /*
    QVector<qreal> additionY = {};

    for(int i = 0; i < Y->size(); ++i)
      additionY.append((*Y)[i] + _kernelPrognosisDerivativeValues[i]);

    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, additionY);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)
        ->setPen(QPen(Qt::green));
    */
}

double MainWindow::countNewtonianDerivative(int i, const QVector<qreal> *Y)
{
  if(i + 1 < Y->size())
  {
    double result = 0;

    result += Y->at(i+1) - Y->at(i);

    result /= ui->lineEdit_domainDensity->text().toDouble();

    return result;
  }
  else return 0;
}

void MainWindow::addKernelPrognosisDerivativePlot(const QVector<qreal> *X)
{
  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)
      ->setData(*X, _kernelPrognosisDerivativeValues);
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)
      ->setPen(QPen(Qt::cyan));
}

void MainWindow::countKernelPrognosisDerivativeY(const QVector<qreal> *X)
{
  std::vector<std::shared_ptr<cluster>> currentClusters
      = getClustersForEstimator();

  std::vector<double> prognosisCoefficients;

  _kernelPrognosisDerivativeValues.clear();

  prognosisCoefficients.clear();

  for(auto c : currentClusters)
    prognosisCoefficients.push_back(c->predictionParameters[1]);

  if(prognosisCoefficients.size() == currentClusters.size())
  {
    kernelPrognoser->setAdditionalMultipliers(prognosisCoefficients);
    kernelPrognoser->setClusters(currentClusters);

    for(qreal x: *X)
    {
      QVector<qreal> pt;
      pt.push_back(x);
      _kernelPrognosisDerivativeValues.push_back(
            kernelPrognoser->getValue(&pt)
      );
    }
  }
}

void MainWindow::addSigmoidallyEnhancedEstimationPlot(
    const QVector<qreal> *X, kernelDensityEstimator *estimator)
{
  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, _sigmoidallyEnhancedPlotY);
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(Qt::darkRed));
}

unsigned long long MainWindow::markUncommonClusters()
{
  double x;

  // Clear previously added markers
  ui->widget_plot->clearItems();

  // For each uncommon cluster add a red vertical line to the plot
  for(std::shared_ptr<cluster> c : uncommonClusters)
  {
    // Only works for distribution data samples as programmed
    x = std::stod(c->getRepresentative()->attributesValues["Val0"]);

    QCPItemLine *verticalLine = new QCPItemLine(ui->widget_plot);
    verticalLine->start->setCoords(x, 0.05);
    verticalLine->end->setCoords(x, -0.05);
    verticalLine->setPen(QPen(Qt::green));
  }

  return uncommonClusters.size();
}

void MainWindow::markNewTrends()
{
  double x = 0.0;

  // For each uncommon cluster add a red vertical line to the plot
  for(std::shared_ptr<cluster> c : uncommonClusters)
  {
    if(c->predictionParameters[1] > 0 && c->_KDEDerivativeValue > 0)
    {
      // Only works for distribution data samples as programmed
      x = std::stod(c->getRepresentative()->attributesValues["Val0"]);

      QCPItemLine *verticalLine = new QCPItemLine(ui->widget_plot);
      verticalLine->start->setCoords(x, 0.1);
      verticalLine->end->setCoords(x, -0.1);
      verticalLine->setPen(QPen(Qt::blue));
    }
  }
}

void MainWindow::markClustersWithNegativeDerivative()
{
  std::vector<std::shared_ptr<cluster>> consideredClusters = getClustersForEstimator();

  double x = 0;

  for(auto c : consideredClusters)
  {
    if(c->predictionParameters[1] < 0) // Mark with
    {
      // Only works for distribution data samples as programmed
      x = std::stod(c->getRepresentative()->attributesValues["Val0"]);

      QCPItemLine *verticalLine = new QCPItemLine(ui->widget_plot);
      verticalLine->start->setCoords(x, 0.01);
      verticalLine->end->setCoords(x, -0.01);
      verticalLine->setPen(QPen(Qt::red));
    }
  }
}

unsigned long long MainWindow::findUncommonClusters()
{
  uncommonClusters.clear();

  std::vector<std::shared_ptr<cluster>> consideredClusters = getClustersForEstimator();

  for(std::shared_ptr<cluster> c : consideredClusters)
  {
    if(c->_currentKDEValue < _maxEstimatorValueOnDomain * _a)
      uncommonClusters.push_back(c);
  }

  qDebug() << "Considered clusters number: " << consideredClusters.size();
  qDebug() << "Uncommon clusters number: " << uncommonClusters.size();
  qDebug() << "Positional estimator value: " << positionalSecondGradeEstimator;
  qDebug() << "Dynamic comparitor value: " << _a * _maxEstimatorValueOnDomain;

  return uncommonClusters.size();
}

QString MainWindow::formatNumberForDisplay(double number)
{
  // According to PK the number should be displayed as #.######
  QString result = "";

  QStringList splitNumber = QString::number(number, 'f', 7).split(".");
  result += splitNumber[0];

  if(splitNumber.size() == 1) return result;

  result += ".";

  for(int i = 0; i < 6 && i < splitNumber[1].size(); ++i)
    result += splitNumber[1][i];

  return result;
}

void MainWindow::countKDEValuesOnClusters(
  std::shared_ptr<kernelDensityEstimator> estimator)
{
  QVector<qreal> x;

  std::vector<std::shared_ptr<cluster>> consideredClusters = getClustersForEstimator();
  estimator->setClusters(consideredClusters);

  for(std::shared_ptr<cluster> c : consideredClusters)
  {
    x.clear();
    x.push_back(std::stod(c->getRepresentative()->attributesValues["Val0"]));
    double estimatorValueOnCluster = estimator->getValue(&x);
    c->_currentKDEValue = estimatorValueOnCluster;
  }
}

void MainWindow::addTemporalDerivativePlot(
  const QVector<qreal> *X, const QVector<qreal> *Y)
{
  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, *Y);
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(Qt::magenta));
}

void MainWindow::fillStandardDeviations(
  QVector<std::shared_ptr<QVector<qreal>>> *stDevs)
{
    int dimensionsNumber                = ui->spinBox_dimensionsNumber->value(),
        targetFunctionElementsNumber    = ui->tableWidget_targetFunctions->rowCount();

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        stDevs->append(std::make_shared<QVector<qreal>>());

        for(int dimensionIndex = 0; dimensionIndex < dimensionsNumber; ++dimensionIndex)
        {
            stDevs->last().get()->append
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

void MainWindow::fillMeans(QVector<std::shared_ptr<QVector<qreal>>> *means)
{
    int dimensionsNumber = ui->spinBox_dimensionsNumber->value(),
        targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        means->append(std::make_shared<QVector<qreal>>());

        for(int dimensionIndex = 0; dimensionIndex < dimensionsNumber; ++dimensionIndex)
        {
            means->last().get()->append
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

void MainWindow::fillDomain(
  QVector<std::shared_ptr<point>>* domain,
  std::shared_ptr<point> *prototypePoint)
{
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

    while(val <= ui->lineEdit_maxX->text().toDouble())
    {
        pPoint.get()->append(val);

        if(pPoint.get()->size() == ui->spinBox_dimensionsNumber->value())
        {
            domain->append(std::make_shared<point>());

            foreach(qreal dimensionVal, *(pPoint.get()))
            {
                domain->last()->append(dimensionVal);
            }
        }
        else
        {
            fillDomain(domain, prototypePoint);
        }

        pPoint.get()->removeLast();

        val += ui->lineEdit_domainDensity->text().toDouble();
    }
}

distribution* MainWindow::generateTargetDistribution(
  QVector<std::shared_ptr<QVector<qreal>>> *means,
  QVector<std::shared_ptr<QVector<qreal>>> *stDevs)
{
    int seed = ui->lineEdit_seed->text().toInt();

    QVector<qreal> contributions;
    QVector<std::shared_ptr<distribution>> elementalDistributions;

    int targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        contributions.append
        (
            (static_cast<QLineEdit*>(ui->tableWidget_targetFunctions->cellWidget(functionIndex, CONTRIBUTION_COLUMN_INDEX)))
            ->text().toDouble()
        );

        elementalDistributions.append(std::shared_ptr<distribution>(new normalDistribution(seed, means->at(functionIndex).get(), stDevs->at(functionIndex).get())));
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
    QVector<int> kernelsIDs;
    QVector<qreal> smoothingParameters;
    QVector<QString> carriersRestrictions;

    for(int rowNumber = 0; rowNumber < dimensionsNumber; ++rowNumber)
    {
        kernelsIDs.append((static_cast<QComboBox*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, KERNEL_COLUMN_INDEX)))->currentIndex());
        smoothingParameters.append((static_cast<QLineEdit*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->text().toDouble());
        carriersRestrictions.append((static_cast<QLineEdit*>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX)))->text());
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
    QVector<std::shared_ptr<QVector<qreal>>>* means,
    QVector<std::shared_ptr<QVector<qreal>>>* stDevs)

{
  QVector<qreal> contributions;
  QVector<std::shared_ptr<function>> elementalFunctions;

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

      contributions.append
      (
          (static_cast<QLineEdit*>(ui->tableWidget_targetFunctions->cellWidget(functionIndex, CONTRIBUTION_COLUMN_INDEX)))
          ->text().toDouble()
      );


      elementalFunctions.append(std::shared_ptr<function>(new multivariateNormalProbabilityDensityFunction(means->at(functionIndex).get(),
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
                                          2000 /* delay */)
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

  _a = MIN_A;

  _longestStepExecutionInSecs = 0;

  double newWeightA = 0;//0.5 * 1.0 / pow(sampleSize - 1, 2);
  double newWeightB = 0;//0.5 * 1.0 / (sampleSize - 1);

  storedMedoids.push_back(std::vector<std::shared_ptr<cluster>>());
  clusters = &(storedMedoids[0]);

  weightedSilvermanSmoothingParameterCounter smoothingParamCounter(clusters, 0);

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
    newWeightA, newWeightB
  );

  double horizontalOffset = 0.5, verticalOffset = 0.03, verticalStep = 0.04;

  // add text labels at the top:

  std::shared_ptr<QCPItemText> E1000TextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  E1000TextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  E1000TextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  E1000TextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  E1000TextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  E1000TextLabel->setText("E1000 = ");

  qDebug() << "Font family : " << font().family();

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> ES1000TextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  ES1000TextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  ES1000TextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  ES1000TextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  ES1000TextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  ES1000TextLabel->setText("a_E1000xK = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> StDevES1000TextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  StDevES1000TextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  StDevES1000TextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  StDevES1000TextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  StDevES1000TextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  StDevES1000TextLabel->setText("StDev(a...) = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> uTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  uTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  uTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  uTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  uTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  uTextLabel->setText("u = ");

  std::vector<std::shared_ptr<QCPItemText>> vs = {};
  std::vector<QString> vsLabels = {"v300 = ", "v500 = ", "v700 = "};

  for(unsigned int i = 0; i < 3; ++i){
    verticalOffset += verticalStep;

    vs.push_back(std::make_shared<QCPItemText>(ui->widget_plot));
    vs.back()->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
    vs.back()->position->setType(QCPItemPosition::ptAxisRectRatio);
    vs.back()->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
    vs.back()->setFont(QFont(font().family(), 28)); // make font a bit larger
    vs.back()->setText(vsLabels[i]);
  }

  verticalOffset += verticalStep;

  /*
  std::shared_ptr<QCPItemText> errorEjTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorEjTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorEjTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorEjTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorEjTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorEjTextLabel->setText("");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorEjpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorEjpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorEjpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorEjpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorEjpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorEjpTextLabel->setText("");

  verticalOffset += verticalStep;
  */

  std::shared_ptr<QCPItemText> error1SejTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  error1SejTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  error1SejTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  error1SejTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  error1SejTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  error1SejTextLabel->setText("ser1_ej = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> error1SejpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  error1SejpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  error1SejpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  error1SejpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  error1SejpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  error1SejpTextLabel->setText("ser1_ejp = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> error2SejTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  error2SejTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  error2SejTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  error2SejTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  error2SejTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  error2SejTextLabel->setText("ser2_ej = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> error2SejpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  error2SejpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  error2SejpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  error2SejpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  error2SejpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  error2SejpTextLabel->setText("ser2_ejp = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorSupSejTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorSupSejTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorSupSejTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorSupSejTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorSupSejTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorSupSejTextLabel->setText("sersup_ej = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorSupSejpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorSupSejpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorSupSejpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorSupSejpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorSupSejpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorSupSejpTextLabel->setText("sersup_ejp = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorModSejTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorModSejTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorModSejTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorModSejTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorModSejTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorModSejTextLabel->setText("sermod_ej = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorModSejpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorModSejpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorModSejpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorModSejpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorModSejpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorModSejpTextLabel->setText("sermod_ejp = ");

  verticalOffset += verticalStep;

  /*
  std::shared_ptr<QCPItemText> errorREjTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorREjTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorREjTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorREjTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorREjTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorREjTextLabel->setText("");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorREjpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorREjpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorREjpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorREjpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorREjpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorREjpTextLabel->setText("");

  verticalOffset += verticalStep;
  */

  std::shared_ptr<QCPItemText> error1SRejTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  error1SRejTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  error1SRejTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  error1SRejTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  error1SRejTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  error1SRejTextLabel->setText("srej1_ej = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> error1SRejpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  error1SRejpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  error1SRejpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  error1SRejpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  error1SRejpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  error1SRejpTextLabel->setText("srej1_ejp = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> error2SRejTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  error2SRejTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  error2SRejTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  error2SRejTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  error2SRejTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  error2SRejTextLabel->setText("srej2_ej = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> error2SRejpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  error2SRejpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  error2SRejpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  error2SRejpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  error2SRejpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  error2SRejpTextLabel->setText("srej2_ejp = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorSupSRejTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorSupSRejTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorSupSRejTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorSupSRejTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorSupSRejTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorSupSRejTextLabel->setText("srejsup_ej = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorSupSRejpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorSupSRejpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorSupSRejpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorSupSRejpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorSupSRejpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorSupSRejpTextLabel->setText("srejsup_ej = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorModSRejTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorModSRejTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorModSRejTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorModSRejTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorModSRejTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorModSRejTextLabel->setText("srejmod_ej = ");

  verticalOffset += verticalStep;

  std::shared_ptr<QCPItemText> errorModSRejpTextLabel =
      std::make_shared<QCPItemText>(ui->widget_plot);
  errorModSRejpTextLabel->setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  errorModSRejpTextLabel->position->setType(QCPItemPosition::ptAxisRectRatio);
  errorModSRejpTextLabel->position->setCoords(horizontalOffset, verticalOffset); // place position at center/top of axis rect
  errorModSRejpTextLabel->setFont(QFont(font().family(), 28)); // make font a bit larger
  errorModSRejpTextLabel->setText("srejmod_ej = ");


  // Save data
  /*
  std::ofstream experimentDataFile;

  experimentDataFile.open("d:\\Dysk Google\\Badania\\experimentData.csv");
  experimentDataFile << "x,y,y_est,e,~z,~~z,z,w\n";

  experimentDataFile.close();
  */

  for(stepNumber = 0; stepNumber < stepsNumber; ++stepNumber)
  {
    clock_t executionStartTime = clock();

    DESDAAlgorithm.performStep();

    ui->label_h_parameter_value->setText(
      QString::number(smoothingParamCounter.getSmoothingParameterValue())
    );

    ui->label_sigma_value->setText(
      QString::number(smoothingParamCounter._stDev)
    );

    targetFunction.reset(generateTargetFunction(&means, &stDevs));

    if(stepNumber > 0 && stepNumber % 10 == 0)
    {
      qDebug() << "Drawing in step number " << stepNumber << ".";
      qDebug() << "h_i = " << smoothingParamCounter.getSmoothingParameterValue();
      qDebug() << "sigma_i = " << smoothingParamCounter._stDev;

      QVector<qreal> X;
      QVector<std::shared_ptr<point>> domain;
      fillDomain(&domain, nullptr);

      for(auto x : domain) X.append(x->at(0));

      _kernelPrognosisDerivativeValues =
          DESDAAlgorithm.getKernelPrognosisDerivativeValues(&X);
      _sigmoidallyEnhancedPlotY =
          DESDAAlgorithm.getEnhancedKDEValues(&X);

      double maxKDEP = 0.0;

      //for(auto val : _sigmoidallyEnhancedPlotY)
      for(int i = 0; i < _sigmoidallyEnhancedPlotY.size(); ++i)
      {
        double val = _sigmoidallyEnhancedPlotY[i];

        if(val > maxKDEP)
        {
          maxKDEP = val;
          KDEPExtrema = X[i];
        }
      }

      drawPlots(estimator.get(), targetFunction.get());


      cluster e1000 = DESDAAlgorithm.getE1000Cluster();

      // E1000
      double avg1000 = e1000._currentKDEValue;

      // StDev(a...)
      double stDevE1000 = 0;

      for(auto val : _kernelPrognosisDerivativeValues)
        stDevE1000 += val;

      double denominator = _kernelPrognosisDerivativeValues.size() - 1;

      stDevE1000 *= -stDevE1000;
      stDevE1000 /= _kernelPrognosisDerivativeValues.size() * denominator;

      for(auto val : _kernelPrognosisDerivativeValues){
        stDevE1000 += val * val / denominator;
      }

      stDevE1000 = sqrt(stDevE1000);

      // a_E1000xK
      double avg1000Est = e1000.predictionParameters[1];

      QVector<qreal> errorHolder = {};

      // ser1_ej + sersup_ej

      double sup_ej = 0.0, sup_ejp = 0.0;

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

      double error2EJ = 0.0, error2EJP = 0.0;

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

      if(stepNumber >= 1000)
      {
        _summaricKDEError1    += _errorEJ;
        _summaricKDEPError1   += _errorEJP;
        _summaricKDEError2    += error2EJ;
        _summaricKDEPError2   += error2EJP;
        _summaricKDEErrorSup  += sup_ej;
        _summaricKDEPErrorSup += sup_ejp;
        _summaricKDEErrorMod  += fabs(modelExtrema - KDEExtrema);
        _summaricKDEPErrorMod += fabs(modelExtrema - KDEPExtrema);
      }

      double rejej = 0.0, rejejp = 0.0, rejsupej = 0.0, rejsupejp = 0.0;

      while(_lastSigmoidallyEnhancedPlotY.size() < _sigmoidallyEnhancedPlotY.size())
      {
        _lastSigmoidallyEnhancedPlotY.push_back(0);
        lastKDEValues.push_back(0);
      }

      errorHolder.clear();

      // Rejej + rejsupej
      for(unsigned int i = 0; i < KDEValues.size(); ++i){
        double val = fabs(KDEValues[i] - lastKDEValues[i]);
        errorHolder.push_back(val);
        rejsupej = val > rejsupej ? val : rejsupej;
      }

      rejej = numericIntegral(&errorHolder);
      errorHolder.clear();

      // Rejejp + rejsupejp
      for(unsigned int i = 0; i < _sigmoidallyEnhancedPlotY.size(); ++i){
        double val = fabs(_sigmoidallyEnhancedPlotY[i] - _lastSigmoidallyEnhancedPlotY[i]);
        errorHolder.push_back(val);
        rejsupejp = val > rejsupejp ? val : rejsupejp;
      }

      rejejp = numericIntegral(&errorHolder);
      errorHolder.clear();

      double rejej2 = 0.0, rejejp2 = 0.0;

      // Rejej
      for(unsigned int i = 0; i < KDEValues.size(); ++i){
        errorHolder.push_back(pow(KDEValues[i] - lastKDEValues[i], 2));
      }

      rejej2 = numericIntegral(&errorHolder);
      errorHolder.clear();

      // Rejejp
      for(unsigned int i = 0; i < _sigmoidallyEnhancedPlotY.size(); ++i){
        errorHolder.push_back(pow(_sigmoidallyEnhancedPlotY[i] - _lastSigmoidallyEnhancedPlotY[i], 2));
      }

      rejejp2 = numericIntegral(&errorHolder);
      errorHolder.clear();

      if(stepNumber >= 1000)
      {
        _summaricKDEsError1    += rejej;
        _summaricKDEPsError1   += rejejp;
        _summaricKDEsError2    += rejej2;
        _summaricKDEPsError2   += rejejp2;
        _summaricKDEsErrorSup  += rejsupej;
        _summaricKDEPsErrorSup += rejsupejp;
        _summaricKDEsErrorMod  += fabs(KDEExtrema - lastKDEExtrema);
        _summaricKDEPsErrorMod += fabs(KDEPExtrema - lastKDEPExtrema);
      }


      E1000TextLabel
          ->setText("E1000         = " + formatNumberForDisplay(avg1000));
      ES1000TextLabel
          ->setText("a_E1000 x K   = " + formatNumberForDisplay(avg1000Est / progressionSize));
      StDevES1000TextLabel
          ->setText("stDev(a ... ) = " + formatNumberForDisplay(stDevE1000));
      uTextLabel
          ->setText("u    = " + formatNumberForDisplay(DESDAAlgorithm._u_i));

      for(unsigned int i = 0; i < DESDAAlgorithm._selectedVValues.size(); ++i){
        vs[i]->setText(
          vsLabels[i] + formatNumberForDisplay(DESDAAlgorithm._selectedVValues[i])
        );
      }

      /*
      errorEjTextLabel
          ->setText("er_ej     = " + formatNumberForDisplay(_errorEJ));
      errorEjpTextLabel
          ->setText("er_ejp   = " + formatNumberForDisplay(_errorEJP));
      */
      error1SejTextLabel
          ->setText("ser1_ej     = " + formatNumberForDisplay(_summaricKDEError1));
      error1SejpTextLabel
          ->setText("ser1_ejp    = " + formatNumberForDisplay(_summaricKDEPError1));
      error2SejTextLabel
          ->setText("ser2_ej     = " + formatNumberForDisplay(_summaricKDEError2));
      error2SejpTextLabel
          ->setText("ser2_ejp    = " + formatNumberForDisplay(_summaricKDEPError2));
      errorSupSejTextLabel
          ->setText("sersup_ej   = " + formatNumberForDisplay(_summaricKDEErrorSup));
      errorSupSejpTextLabel
          ->setText("sersup_ejp  = " + formatNumberForDisplay(_summaricKDEPErrorSup));
      errorModSejTextLabel
          ->setText("sermod_ej   = " + formatNumberForDisplay(_summaricKDEErrorMod));
      errorModSejpTextLabel
          ->setText("sermod_ejp  = " + formatNumberForDisplay(_summaricKDEPErrorMod));
      /*
      errorREjTextLabel
          ->setText("rej_ej    = " + formatNumberForDisplay(rejej));
      errorREjpTextLabel
          ->setText("rej_ejp  = " + formatNumberForDisplay(rejejp));
      */
      error1SRejTextLabel
          ->setText("srej1_ej    = " + formatNumberForDisplay(_summaricKDEsError1));
      error1SRejpTextLabel
          ->setText("srej1_ejp   = " + formatNumberForDisplay(_summaricKDEPsError1));
      error2SRejTextLabel
          ->setText("srej2_ej    = " + formatNumberForDisplay(_summaricKDEsError2));
      error2SRejpTextLabel
          ->setText("srej2_ejp   = " + formatNumberForDisplay(_summaricKDEPsError2));
      errorSupSRejTextLabel
          ->setText("srejsup_ej  = " + formatNumberForDisplay(_summaricKDEsErrorSup));
      errorSupSRejpTextLabel
          ->setText("srejsup_ejp = " + formatNumberForDisplay(_summaricKDEPsErrorSup));
      errorModSRejTextLabel
          ->setText("srejmod_ej  = " + formatNumberForDisplay(_summaricKDEsErrorMod));
      errorModSRejpTextLabel
          ->setText("srejmod_ejp = " + formatNumberForDisplay(_summaricKDEPsErrorMod));

      _lastSigmoidallyEnhancedPlotY.clear();
      lastKDEValues.clear();

      for(auto val : _sigmoidallyEnhancedPlotY)
        _lastSigmoidallyEnhancedPlotY.push_back(val);

      for(auto val : KDEValues)
        lastKDEValues.push_back(val);

      lastModelExtrema = modelExtrema;
      lastKDEExtrema = KDEExtrema;
      lastKDEPExtrema = KDEPExtrema;

      qApp->processEvents();

      QString dirPath = "D:\\Dysk Google\\TR Badania\\Eksperyment 63\\";
      //QString dirPath = "D:\\Dysk Google\\TR Badania\\test\\";

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
}
