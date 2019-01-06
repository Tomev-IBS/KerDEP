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
}

void MainWindow::testNewFunctionalities()
{
  qDebug() << "Start test.";
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

  // Clusters are now pointer to 0 level of storedMedoids
  /*
  for(std::shared_ptr<cluster> c : *clusters)
  {
    if(c->getWeight() >= positionalSecondGradeEstimator)
      consideredClusters.push_back(c);
  }
  */

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

void MainWindow::updateA()
{
  std::vector<std::shared_ptr<cluster>> allConsideredClusters = getClustersForEstimator();
  double summaricClustersWeight = 0.0;

  for(std::shared_ptr<cluster> c : allConsideredClusters)
  {
    summaricClustersWeight += c->getWeight();
  }

  double currentUncommonClusterWeight = 0.0;

  for(std::shared_ptr<cluster> uc : uncommonClusters)
  {
    currentUncommonClusterWeight += uc->getWeight();
  }

  currentUncommonClusterWeight /= summaricClustersWeight;

  qDebug() << "Previous uncommon clusters part: " << _previousUncommonClustersWeight;
  qDebug() << "Current uncommon clusters part: " << currentUncommonClusterWeight;
  qDebug() << "Desired uncommon clusters part: " << ui->lineEdit_rarity->text().toDouble();
  qDebug() << "Previous a value: " << _a;

  _a += (ui->lineEdit_rarity->text().toDouble() - currentUncommonClusterWeight)
            * _maxEstimatorValueOnDomain;

  if(_a > MAX_A) _a = MAX_A;
  if(_a < MIN_A) _a = MIN_A;

  qDebug() << "Current a value: " << _a;

  _previousUncommonClustersWeight = currentUncommonClusterWeight;
}

void MainWindow::drawPlots(kernelDensityEstimator* estimator,
                           function* targetFunction)
{
    clearPlot();

    resizePlot();

    QVector<std::shared_ptr<point>> domain;

    QVector<qreal> X;
    QVector<qreal> normalDistributionY;

    // Fill domain with points
    // To keep things simple let's consider only these domains wherein
    // each dimension has equal size.

    fillDomain(&domain, nullptr);

    foreach(auto x, domain)
    {
        normalDistributionY.append(targetFunction->getValue(x.get()));
        X.append(x->at(0));
    }

    // Generate plot of model function
    if(ui->checkBox_showEstimatedPlot->isChecked())
    addModelPlot(&X, &normalDistributionY);

    // Generate a vector of values from selected KDE
    KDEEstimationY.clear();

    double val;
    _maxEstimatorValueOnDomain = 0;

    // TODO: Place counting in another thread
    foreach(std::shared_ptr<point> x, domain)
    {
      val = estimator->getValue(x.get());

      if (val > _maxEstimatorValueOnDomain) _maxEstimatorValueOnDomain = val;

      oldKerernelY.append(val);
      KDEEstimationY.append(val);
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
      countKernelPrognosisDerivativeY(&X);
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

void MainWindow::addPlot(const QVector<qreal> *X, const QVector<qreal> *Y)
{
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, *Y);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)
        ->setPen(QPen(getRandomColor()));
}

void MainWindow::addModelPlot(const QVector<qreal> *X, const QVector<qreal> *Y)
{
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
}

double MainWindow::countNewtonianDerivative(int i, const QVector<qreal> *Y)
{
  if(i + 1 < Y->size())
  {
    double result = 0;

    result += Y->at(i+1) - Y->at(i);

    result /= this->ui->lineEdit_domainDensity->text().toDouble();

    return result;
  }
  else return 0;
}

void MainWindow::addKernelPrognosisDerivativePlot(const QVector<qreal> *X)
{
  std::vector<std::shared_ptr<cluster>> currentClusters
      = getClustersForEstimator();

  QVector<double> offsetKernelPrognosisPlotValues;

  for(int i = 0; i < _kernelPrognosisDerivativeValues.size(); ++i)
  {
    offsetKernelPrognosisPlotValues.push_back(
      _kernelPrognosisDerivativeValues[i]
    );
  }

  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)
      ->setData(*X, offsetKernelPrognosisPlotValues);
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
  auto currentClusters = getClustersForEstimator();

  QVector<qreal> sigmoidallyEnhancedPlotY;
  QVector<double> starredVs;
  std::vector<double> vs;
  double starredVsSum = 0.0;

  for(auto c : currentClusters)
  {
    starredVs.push_back(2.0 / (1.0 + exp(- c->predictionParameters[1] * adaptivePredictionPowerParameter )));
    starredVsSum += starredVs.last();
  }

  for(auto starredV : starredVs)
  {
    vs.push_back(starredVs.length() * starredV / starredVsSum);
  }

  estimator->setAdditionalMultipliers(vs);

  for(auto x : *X)
  {
    QVector<qreal> pt = {x};
    sigmoidallyEnhancedPlotY.push_back(estimator->getValue(&pt));
  }

  estimator->setAdditionalMultipliers({});

  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, sigmoidallyEnhancedPlotY);
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(Qt::darkYellow));
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

void MainWindow::generateSamples(
  QVector<std::shared_ptr<QVector<qreal>> > *means,
  QVector<std::shared_ptr<QVector<qreal>> > *stDevs)
{
    samples.clear();
    objects.clear();

    int sampleSize = ui->lineEdit_sampleSize->text().toInt();

    if(sampleSize < 1)
    {
        qDebug() << "Sample size < 1.";
        return;
    }

    std::shared_ptr<distribution> targetDistribution(generateTargetDistribution(means, stDevs));

    dataParser *parser = new distributionDataParser(&attributesData);

    qreal progressionSize = ui->lineEdit_distributionProgression->text().toDouble();

    dataReader *reader = new progressiveDistributionDataReader(targetDistribution.get(), progressionSize, 1000);

    reader->gatherAttributesData(&attributesData);
    parser->setAttributesOrder(reader->getAttributesOrder());

    reservoirSamplingAlgorithm *samplingAlgorithm = generateReservoirSamplingAlgorithm(reader, parser);

    samplingAlgorithm->fillReservoir(&objects);

    foreach(auto object, objects)
    {
      samples.push_back(std::make_shared<QVector<qreal>>());

      for(auto nameValue : static_cast<distributionDataSample*>(object.get())->attributesValues)
      {
        samples.last()->push_back(stod(nameValue.second));
      }
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

QColor MainWindow::getRandomColor()
{
    return QColor(rand()%110 + 50, rand()%110 + 50, rand()%110 + 50);
}

void MainWindow::on_pushButton_animate_clicked()
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

    std::shared_ptr<distribution>
      targetDistribution(generateTargetDistribution(&means, &stDevs));

    parser.reset(new distributionDataParser(&attributesData));

    qreal progressionSize =
        ui->lineEdit_distributionProgression->text().toDouble();

    reader.reset(
      new progressiveDistributionDataReader(targetDistribution.get(),
                                            progressionSize,
                                            1000 /* delay */)
    );

    reader->gatherAttributesData(&attributesData);
    parser->setAttributesOrder(reader->getAttributesOrder());

    positionalSecondGradeEstimatorCountingMethod =
      ui->comboBox_rareElementsMethod->currentIndex();

    reservoirSamplingAlgorithm* algorithm =
      generateReservoirSamplingAlgorithm(reader.get(), parser.get());

    objects.clear();

    int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();

    int numberOfClustersForGrouping = 100;
    int medoidsNumber = 50;

    groupingThread gt(&storedMedoids, parser);

    gt.setAttributesData(&attributesData);

    qDebug() << "Attributes data set.";

    int sampleSize = ui->lineEdit_sampleSize->text().toInt();
    gt.initialize(medoidsNumber, sampleSize);

    _a = MIN_A;

    _longestStepExecutionInSecs = 0;

    storedMedoids.push_back(std::vector<std::shared_ptr<cluster>>());
    clusters = &(storedMedoids[0]);

    weightedSilvermanSmoothingParameterCounter smoothingParamCounter(clusters, 0);

    for(stepNumber = 0; stepNumber < stepsNumber; ++stepNumber)
    {
      clock_t executionStartTime = clock();

      updateWeights();

      algorithm->performSingleStep(&objects, stepNumber);

      // TR TODO: It's not working for biased algorithm
      std::shared_ptr<cluster> newCluster =
          std::shared_ptr<cluster>(new cluster(stepNumber, objects.back()));
      newCluster->setTimestamp(stepNumber);

      clusters->push_back(newCluster);

      //qDebug() << "Counting of smoothing parameters.";
      smoothingParamCounter.updateSmoothingParameterValue(
        ui->lineEdit_weightModifier->text().toDouble(),
        std::stod(clusters->back()->getObject()->attributesValues["Val0"])
      );

      double smoothingParameterMultiplier = 1;
      std::vector<double> smoothingParameters =
      {
        smoothingParamCounter.getSmoothingParameterValue() * smoothingParameterMultiplier
      };

      estimator->setSmoothingParameters(smoothingParameters);

      ui->label_h_parameter_value->setText(
        QString::number(smoothingParameters[0])
      );

      ui->label_sigma_value->setText(
        QString::number(smoothingParamCounter._stDev)
      );

      //qDebug() << "Smoothing params counted.";

      // Write h to file for data analysis
      /*
      std::ofstream myfile;
      //myfile.open("h_params.csv", std::ios_base::app);
      myfile << smoothingParameters[0] << ",";
      myfile.close();
      */

      auto currentClusters = getClustersForEstimator();

      qDebug() << "Reservoir size in step "
                 << stepNumber << " is: " << currentClusters.size();
      //qDebug() << "Objects size: " << objects.size();

      //qDebug() << "Counting KDE values on clusters.";
      countKDEValuesOnClusters(estimator);
      //qDebug() << "Counted. Finding uncommon clusters.";
      findUncommonClusters();

      if(currentClusters.size() >= algorithm->getReservoidMaxSize())
      {
        qDebug() << "============ Clustering function started ============";

        std::vector<std::shared_ptr<cluster>> clustersForGrouping;

        for(int cNum = 0; cNum < numberOfClustersForGrouping; ++cNum)
        {
          clustersForGrouping.push_back((*clusters)[0]);
          clusters->erase(clusters->begin(), clusters->begin()+1);
        }

        objects.erase(objects.begin(), objects.begin() + (numberOfClustersForGrouping - medoidsNumber));
        gt.getClustersForGrouping(clustersForGrouping);
        gt.run();
        currentClusters = getClustersForEstimator();
        clusters = &(storedMedoids[0]); // Reassigment needed, as (probably) grouping algorithm resets it
      }

      //updateA();

      estimator->setClusters(currentClusters);

      //qDebug() << "Updating prognosis.";
      updatePrognosisParameters();
      //qDebug() << "Updated. Counting derivative values for clusters.";
      countKDEDerivativeValuesOnClusters();

      targetFunction.reset(generateTargetFunction(&means, &stDevs));

      //if(stepNumber < 10000)
      //if( stepNumber == 1001 || stepNumber == 1101 || stepNumber == 2101
      //    || stepNumber == 3101 || stepNumber == 4101 || stepNumber == 5101)
      if(stepNumber > 1000 && (stepNumber - 1) % 25 == 0)
      {
        qDebug() << "Drawing in step number " << stepNumber << ".";
        qDebug() << "h_i = " << smoothingParameters[0];
        qDebug() << "sigma_i = " << smoothingParamCounter._stDev;

        drawPlots(estimator.get(), targetFunction.get());

        qApp->processEvents();

        QString imageName = "D:\\Dysk Google\\Badania\\v=0.01, b=0.999\\"
            + QString::number(stepNumber - 1) + ".png";
        ui->widget_plot->savePng(imageName, 0, 0, 1, -1);
      }
      else break;

      targetFunction.reset(generateTargetFunction(&means, &stDevs));

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

void MainWindow::on_pushButton_clear_clicked()
{
    clearPlot();

    ui->widget_plot->replot();
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

void MainWindow::on_pushButton_countSmoothingParameters_clicked()
{
  QVector<std::shared_ptr<QVector<qreal>>> means, stDevs;

  fillMeans(&means);
  fillStandardDeviations(&stDevs);

  generateSamples(&means, &stDevs);

  // Count smoothing parameter for each dimension
  int numberOfRows = ui->tableWidget_dimensionKernels->rowCount();

  QVector<qreal> samplesColumn;
  QVector<int> weights;

  smoothingParameterCounter* counter
      = generateSmoothingParameterCounter(&samplesColumn);

  qreal value;

  for(int rowNumber = 0; rowNumber < numberOfRows; ++rowNumber)
  {
    // Create vector that consists of variables inside this dimension

    samplesColumn.clear();

    for(std::shared_ptr<QVector<qreal>> sample : samples)
        samplesColumn.append(sample.get()->at(rowNumber));

    value = counter->countSmoothingParameterValue();

    (static_cast<QLineEdit*>((ui
                  ->tableWidget_dimensionKernels
                  ->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))
      ->setText(QString::number(value)));
  }
}

smoothingParameterCounter *MainWindow::generateSmoothingParameterCounter(
    QVector<qreal> *samplesColumn)
{
  int smoothingParameterCounterID = ui
                                    ->comboBox_smoothingParameterCountingMethod
                                    ->currentIndex();

  switch (smoothingParameterCounterID)
  {
    case WEIGHTED_SILVERMAN:
      return new weightedSilvermanSmoothingParameterCounter(samplesColumn,
                                                            nullptr);
    case RANK_3_PLUG_IN:
      return new pluginSmoothingParameterCounter(samplesColumn, 3);
    case RANK_2_PLUG_IN:
    default:
      return new pluginSmoothingParameterCounter(samplesColumn, 2);
  }
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

void MainWindow::updateWeights()
{
  double weightModifier = ui->lineEdit_weightModifier->text().toDouble(),
         weightDeletionThreshold = ui->lineEdit_deletionThreshold
                                     ->text().toDouble();

  for(unsigned int level = 0; level < storedMedoids.size(); ++level)
  {
    for(unsigned int medoidNumber = 0; medoidNumber < storedMedoids[level].size(); ++medoidNumber)
    {
      storedMedoids[level][medoidNumber]->setWeight(
        weightModifier * storedMedoids[level][medoidNumber]->getWeight()
      );
    }
  }

  // Clusters are now pointer to 0 level of storedMedoids
  /*
  for(unsigned int i = 0; i < clusters.size(); ++i)
  {
    clusters[i]->setWeight(weightModifier * clusters[i]->getWeight());
    if(clusters[i]->getWeight() < weightDeletionThreshold)
      clusters.erase(clusters.begin() + static_cast<int>(i));
  }
  */
}

void MainWindow::updatePrognosisParameters()
{
  std::vector<std::shared_ptr<cluster>> currentClusters
      = getClustersForEstimator();

  if(currentClusters.size() == 0) return;

  for(std::shared_ptr<cluster> c : currentClusters)
  {
    QVector<qreal> pt;

    // Assuming it has object (as it should) and have only numerical values.
    for(auto kv : c->getObject()->attributesValues)
      pt.append(std::stod(kv.second));

    if(c->predictionParameters.size() > 0)
    {
      c->updateDeactualizationParameter(c->_currentKDEValue);
      c->updatePredictionParameters(c->_currentKDEValue);
    }
    else
    {
      c->initializePredictionParameters(c->_currentKDEValue);
    }

    c->_lastKDEValue = c->_currentKDEValue;
  }

}

void MainWindow::countKDEDerivativeValuesOnClusters()
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

    QVector<qreal> x;

    for(auto c : currentClusters)
    {
      QVector<qreal> pt;
      pt.push_back(std::stod(c->getRepresentative()->attributesValues["Val0"]));

      c->_KDEDerivativeValue = kernelPrognoser->getValue(&pt);
    }
  }
}
