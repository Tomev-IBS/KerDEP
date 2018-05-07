#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <math.h>
#include <climits>
#include <set>
#include <QDebug>
#include <algorithm>

#include "Functions/multivariatenormalprobabilitydensityfunction.h"
#include "Functions/complexfunction.h"

#include "Distributions/distributions.h"
#include "KDE/pluginsmoothingparametercounter.h"
#include "KDE/weightedSilvermanSmoothingParameterCounter.h"

#include "Reservoir_sampling/biasedReservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/basicReservoirSamplingAlgorithm.h"

#include "Reservoir_sampling/distributiondataparser.h"
#include "Reservoir_sampling/progressivedistributiondatareader.h"

#include "Detectors/rareElementsDetector.h"

#include "VDE/velocityDensityEstimator.h"
#include "VDE/VDEThread.h"

#include "groupingThread/groupingThread.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    kernelTypes << "Normal" << "Triangle" << "Epanecznikow" << "Dull";

    ui->setupUi(this);

    setupValidators();
    setupPlot();
    setupKernelsTable();
}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
  int keyCode = event->key();

  //qDebug() << keyCode;

  // If key is a number
  if(keyCode >= 48 && keyCode <= 57)
  {
    int objectsNumber = keyCode - 48;

    if(objectsNumber == 0) objectsNumber += 10;

    insertObjectsBetweenIntervals(objectsNumber);
  }

  // If m is pressed
  if(keyCode == 77) insertMassiveData();

  // If left arrow pressed
  if(keyCode == 16777234)
  {
    qDebug() << "Left arrow pressed.";
    moveTargetFunctionLeft();
  }

  // If right arrow pressed
  if(keyCode == 16777236) qDebug() << "Right arrow pressed.";

  // If down arrow pressed
  if(keyCode == 16777237) qDebug() << "Down arrow pressed.";
}

int MainWindow::insertObjectsBetweenIntervals(unsigned int objectsNumber)
{
  std::vector<std::shared_ptr<sample>> interIntervalObjects;

  generateInterIntervalObjects(&interIntervalObjects, objectsNumber);
  selectDesiredNumberOfInterIntervalObjects(&interIntervalObjects);
  insertClustersFromInterIntervalObjects(&interIntervalObjects);

  objects.insert(objects.end(), interIntervalObjects.begin(), interIntervalObjects.end());

  qDebug() << "Inter interval objects inserterd: " << interIntervalObjects.size();

  return interIntervalObjects.size();
}

int MainWindow::generateInterIntervalObjects(std::vector<std::shared_ptr<sample> > *interIntervalObjects, unsigned int objectsNumber)
{
  //TR TODO: Check if reader and parser are initialized

  for(unsigned int i = 0; i < objectsNumber; ++i)
  {
    reader->getNextRawDatum(parser->buffer);

    parser->addDatumToContainer(interIntervalObjects);

    parser->writeDatumOnPosition(interIntervalObjects, interIntervalObjects->size()-1);
  }

  return interIntervalObjects->size();
}

int MainWindow::selectDesiredNumberOfInterIntervalObjects(std::vector<std::shared_ptr<sample> > *interIntervalObjects)
{
  unsigned int desiredNumberOfClusters =
      ui->lineEdit_interIntervalClusters->text().toInt();

  while(interIntervalObjects->size() > desiredNumberOfClusters)
    interIntervalObjects->erase(interIntervalObjects->begin() + (rand() % interIntervalObjects->size()));

  return interIntervalObjects->size();
}

int MainWindow::insertClustersFromInterIntervalObjects(std::vector<std::shared_ptr<sample> > *interIntervalObjects)
{
  std::vector<std::shared_ptr<cluster>> newClusters;

  for(unsigned int i = 0; i < interIntervalObjects->size(); ++i)
  {
    newClusters.push_back(std::shared_ptr<cluster>(new cluster(clusters.size()+i, (*interIntervalObjects)[i])));
    newClusters.back()->setTimestamp(stepNumber);
  }

  setInterIntervalClustersWeights(&newClusters);

  clusters.insert(clusters.end(), newClusters.begin(), newClusters.end());
  clustersForVDE.insert(clustersForVDE.end(), newClusters.begin(), newClusters.end());

  return newClusters.size();
}

double MainWindow::setInterIntervalClustersWeights(std::vector<std::shared_ptr<cluster> > *newClusters)
{
  double weight = countInterIntervalClustersWeight();

  for(unsigned int i = 0; i < newClusters->size(); ++i)
    (*newClusters)[i]->setWeight(weight);

  return weight;
}

double MainWindow::countInterIntervalClustersWeight()
{
  // All values in milliseconds.

  long intervalValue = ui->lineEdit_milisecondsDelay->text().toInt();

  if(intervalValue == 0) return 1.0;

  long end = std::chrono::duration_cast< std::chrono::milliseconds >(
        std::chrono::system_clock::now().time_since_epoch()).count();

  long difference = end - start;

  double weightModifier = ui->lineEdit_weightModifier->text().toDouble();

  double power = 1 -((double)difference)/((double)intervalValue);


  // TR TODO: This is neccessary for proper functionality when intervals are low.

  if(power < 0)
  {
    qDebug() << "Power is below 0.";
    return 0.0;
  }


  return pow(weightModifier, power);
}

int MainWindow::insertMassiveData()
{
  std::vector<std::shared_ptr<sample>> massiveData;

  qDebug() << "Generating data.";

  generateMassiveData(&massiveData);

  qDebug() << "Massive data generated.";

  clusterMassiveData(&massiveData, &storedMedoids);

  qDebug() << "Massive data clustered.";

  return massiveData.size();
}

int MainWindow::generateMassiveData(std::vector<std::shared_ptr<sample>> *dataContainer)
{
  long dataSize = 10000;

  dataContainer->clear();

  for(int i = 0; i < dataSize; ++i)
  {
    reader->getNextRawDatum(parser->buffer);

    parser->addDatumToContainer(dataContainer);

    parser->writeDatumOnPosition(dataContainer, dataContainer->size()-1);
  }

  return dataContainer->size();
}

void MainWindow::clusterMassiveData(std::vector<std::shared_ptr<sample>> *objects,
                                    std::vector<std::vector<std::shared_ptr<cluster>>> *storage)
{
  // Select medoids
  std::set<int> medoidsIndexes;

  // TODO TR: Add ui control
  unsigned int medoidsNumber = 10;

  do
  {
    medoidsIndexes.insert(rand() % objects->size());
  } while(medoidsIndexes.size() < medoidsNumber);

  // Create clusters from medoids
  std::set<int>::iterator it = medoidsIndexes.begin();

  if(storage->size() == 0) storage->push_back(std::vector<std::shared_ptr<cluster>>());

  for(unsigned int i = 0; i < medoidsNumber; ++i)
  {
    storage->at(0).push_back(std::make_shared<cluster>(cluster(i, objects->at(*it))));
    storage->at(0).back()->setTimestamp(stepNumber);
    clustersForVDE.push_back(storage->at(0).back());
    storage->at(0).back()->setWeight(objects->size() / medoidsNumber);
    std::advance(it, 1);
  }
}

std::vector<std::shared_ptr<cluster>> MainWindow::getClustersForEstimator()
{
  std::vector<std::shared_ptr<cluster>> consideredClusters;

  for(std::shared_ptr<cluster> c : clusters)
  {
    if(c->getWeight() >= positionalSecondGradeEstimator)
      consideredClusters.push_back(c);
  }

  for(std::vector<std::shared_ptr<cluster>> level : storedMedoids)
  {
    for(std::shared_ptr<cluster> c : level)
    {
      if(c->getWeight() > positionalSecondGradeEstimator)
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
    ui->tableWidget_dimensionKernels->horizontalHeader()->setStretchLastSection(true);

    refreshKernelsTable();
    refreshTargetFunctionTable();
}

void MainWindow::drawPlots(kernelDensityEstimator* estimator, function* targetFunction)
{
    // Check if prior plots should be saved
    if(!ui->checkBox_keepPriorPlots->isChecked())
    {
        // If not clear plot
        clearPlot();
    }

    resizePlot();

    QVector<std::shared_ptr<point>> domain;

    QVector<qreal> X;
    QVector<qreal> normalDistributionY;

    // Fill domain with points
    // To keep things simple let's consider only these domains wherein
    // each dimension has equal size.

    fillDomain(&domain, NULL);

    foreach(auto x, domain)
    {
        normalDistributionY.append(targetFunction->getValue(x.get()));
        X.append(x->at(0));
    }

    // Generate plot of model function
    if(ui->checkBox_showEstimatedPlot->isChecked())
    addModelPlot(&X, &normalDistributionY);

    // Generate a vector of values from selected KDE
    QVector<qreal> KDEEstimationY;

    double val;

    // TODO: Place counting in another thread
    foreach(std::shared_ptr<point> x, domain)
    {
      val = estimator->getValue(x.get());

      oldKerernelY.append(val);
      KDEEstimationY.append(val);
    }

    // Generate a plot of KDE
    if(ui->checkBox_showEstimatedPlot->isChecked())
      addEstimatedPlot(&X, &KDEEstimationY);

    // Generate a plot of temporal derivative

    KDETemporalDerivativeY.clear();
    double visibilityEnchantCoefficient = 1;

    if(oldKerernelY.size() != 0)
    {
      for(int i = 0; i < KDEEstimationY.size(); ++i)
        KDETemporalDerivativeY.push_back(visibilityEnchantCoefficient*
                                         (KDEEstimationY[i] - oldKerernelY[i]));
    }

    if(ui->checkBox_showTimeDerivativePlot->isChecked())
      addTemporalDerivativePlot(&X, &KDETemporalDerivativeY);

    // Generate plot for estimated KDE in the next iteration
    if(ui->checkBox_showPrognosedPlot->isChecked())
      addPrognosedEstimationPlots(&X, &KDEEstimationY);

    if(ui->checkBox_kernelPrognosedPlot->isChecked())
      addKernelPrognosedEstimationPlot();

    if(ui->checkBox_showUnusualClusters->isChecked())
      markUncommonClusters(estimator);

    if(ui->checkBox_showNewTrends->isChecked())
      markNewTrends();

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
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(getRandomColor()));
}

void MainWindow::addModelPlot(const QVector<qreal> *X, const QVector<qreal> *Y)
{
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, *Y);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(Qt::red));
}

void MainWindow::addEstimatedPlot(const QVector<qreal> *X, const QVector<qreal> *Y)
{
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, *Y);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(Qt::blue));
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

void MainWindow::addPrognosedEstimationPlots(const QVector<qreal> *X, const QVector<qreal> *KDEY)
{
  predictKDEValues(X, KDEY);

  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, predictedKDEValues);
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(Qt::green));
}

int MainWindow::predictKDEValues(const QVector<qreal> *X, const QVector<qreal> *KDEY)
{

  //qDebug() << "Updating predictions.";

  updatePointsPredictionParameters(KDEY, &predictedKDEValues, &pointsPredictionParameters);

  QVector<qreal> newPredictions;

  //qDebug() << "Counting new predictions.";

  for(unsigned int i = 0; i < pointsPredictionParameters.size(); ++i)
    newPredictions.push_back(pointsPredictionParameters[i][0] + stepNumber * pointsPredictionParameters[i][1]);

  //qDebug() << "Swapping old predictions with new.";

  predictedKDEValues = newPredictions;

  return predictedKDEValues.size();
}

int MainWindow::updatePointsPredictionParameters(const QVector<qreal> *KDEY,
                                                 QVector<double>* predictedValues,
                                                 std::vector<std::vector<double>>* target)
{
  if(target->size() == 0)
  {
    countInitialPredictionParameters(KDEY, &pointsPredictionParameters);
    return pointsPredictionParameters.size();
  }

  std::vector<std::vector<double>> updatedParameters;

  double upperValue, lowerValue;

  for(unsigned int i = 0; i < pointsPredictionParameters.size(); ++i)
  {
    upperValue = pointsPredictionParameters[i][0] + pointsPredictionParameters[i][1];
    lowerValue = pointsPredictionParameters[i][1];
    upperValue += (1 - pow(deactualizationParameter, 2)) * (KDEY->at(i)- predictedKDEValues[i]);
    lowerValue += pow(1 - deactualizationParameter, 2) * (KDEY->at(i) - predictedKDEValues[i]);
    updatedParameters.push_back(std::vector<double>({
      upperValue, lowerValue
    }));
  }

  pointsPredictionParameters = updatedParameters;

  return pointsPredictionParameters.size();
}

int MainWindow::countInitialPredictionParameters(const QVector<qreal> *KDEY,
                                                 std::vector<std::vector<double>>* target)
{
  std::vector<std::vector<double>> reversedD =
    { {1-pow(deactualizationParameter,2), pow((1- deactualizationParameter), 2)},
      {pow((1- deactualizationParameter), 2), pow((1- deactualizationParameter), 3)/ deactualizationParameter}};

  for(unsigned int i = 0; i < KDEY->size(); ++i)
  {
    target->push_back(std::vector<double>({
      reversedD[0][0] * KDEY->at(i),
      reversedD[1][0] * KDEY->at(i)
    }));
  }

  return target->size();
}

void MainWindow::addKernelPrognosedEstimationPlot()
{
  qDebug() << "Adding Kernel estimation plot.";

  std::vector<std::shared_ptr<cluster>> currentClusters
      = getClustersForEstimator();
  std::vector<std::shared_ptr<cluster>> prognosedClusters;

  for(std::shared_ptr<cluster> c : currentClusters)
  {
    prognosedClusters.push_back(
      std::shared_ptr<cluster>(new cluster(c->getObject()))
    );
  }

  std::vector<double> prognosisCoefficients;

  //ui->widget_plot->addGraph();
  //ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, predictedKDEValues);
  //ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(Qt::yellow));
}

int MainWindow::markUncommonClusters(kernelDensityEstimator* estimator)
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
    verticalLine->setPen(QPen(Qt::red));
  }

  return uncommonClusters.size();
}

int MainWindow::markNewTrends()
{
  double x;
  int trendStepsRequired = ui->spinBox_trendStep->value();

  // For each uncommon cluster add a red vertical line to the plot
  for(std::shared_ptr<cluster> c : uncommonClusters)
  {
    if(c->positiveTemporalDerivativeTimesInARow >= trendStepsRequired)
    {
      // Only works for distribution data samples as programmed
      x = std::stod(c->getRepresentative()->attributesValues["Val0"]);

      QCPItemStraightLine *verticalLine = new QCPItemStraightLine(ui->widget_plot);
      verticalLine->point1->setCoords(x, MIN_Y);
      verticalLine->point2->setCoords(x, MAX_Y);
      verticalLine->setPen(QPen(Qt::green));
    }
  }

}

int MainWindow::findUncommonClusters(kernelDensityEstimator* estimator)
{
  uncommonClusters.clear();

  QVector<qreal> x;

  std::vector<std::shared_ptr<cluster>> consideredClusters = getClustersForEstimator();
  estimator->setClusters(consideredClusters);

  for(std::shared_ptr<cluster> c : consideredClusters)
  {
    x.clear();
    x.push_back(std::stod(c->getRepresentative()->attributesValues["Val0"]));
    if(estimator->getValue(&x) < positionalSecondGradeEstimator)
      uncommonClusters.push_back(c);
  }

  qDebug() << "Considered clusters number: " << consideredClusters.size();
  qDebug() << "Uncommon clusters number: " << uncommonClusters.size();
  qDebug() << "Positional estimator value: " << positionalSecondGradeEstimator;

  return uncommonClusters.size();
}

std::vector<double> MainWindow::countUnsortedReducedEstimatorValuesOnEstimatorClusters(kernelDensityEstimator *estimator)
{
  std::vector<double> unsortedReducedEstimatorValues;
  std::vector<std::shared_ptr<cluster>> consideredClusters =
      getClustersForEstimator();
  std::shared_ptr<cluster> reducedCluster;
  QVector<qreal> x;

  if(consideredClusters.size() == 1)
  {
    unsortedReducedEstimatorValues.push_back(1.0);
    return unsortedReducedEstimatorValues;
  }

  for(unsigned int i = 0; i < consideredClusters.size(); ++i)
  {
    reducedCluster = consideredClusters[0];
    consideredClusters.erase(consideredClusters.begin());
    estimator->setClusters(consideredClusters);
    x.clear();
    x.push_back(std::stod(reducedCluster->getRepresentative()->attributesValues["Val0"]));
    unsortedReducedEstimatorValues.push_back(estimator->getValue(&x));
    consideredClusters.push_back(reducedCluster);
  }

  return unsortedReducedEstimatorValues;
}

double MainWindow::countPositionalSecondGradeEstimator(
    std::vector<double> *unsortedReducedEstimatorValuesOnClusters)
{
  double uncommonnessThreshold = ui->lineEdit_rarity->text().toDouble();
  double mr = uncommonnessThreshold;

  double estimator = 0;

  double j;

  switch (positionalSecondGradeEstimatorCountingMethod)
  {
    case WEIGHTED:
      mr *= getSummaricClustersWeight(getClustersForEstimator());
      j = mr;
    break;
    case STANDARD:
    default:
      mr *= unsortedReducedEstimatorValuesOnClusters->size();
      j = mr + 0.5;
    break;
  }

  std::vector<double> jSortedReducedEstimatorValues =
      sortJReducedEstimatorValues(unsortedReducedEstimatorValuesOnClusters, j+1);

  if(mr < 0.5)
    estimator = jSortedReducedEstimatorValues[0];
  else
  {
    switch (positionalSecondGradeEstimatorCountingMethod)
    {
      case WEIGHTED:
        estimator = jSortedReducedEstimatorValues[jSortedReducedEstimatorValues.size() - 1];
      break;
      case STANDARD:
      default:
        estimator = (0.5 + j - mr) * jSortedReducedEstimatorValues[j-1];
        estimator += (0.5 - j + mr) * jSortedReducedEstimatorValues[j];
      break;
    }
  }

  return estimator;
}

double MainWindow::getSummaricClustersWeight(std::vector<std::shared_ptr<cluster> > clusters)
{
  double result = 0.0;

  for(std::shared_ptr<cluster> c : clusters) result += c->getWeight();

  return result;
}

std::vector<double> MainWindow::sortJReducedEstimatorValues(
    std::vector<double> *unsortedReducedEstimatorValuesOnClusters, double j)
{
  std::vector<double> jSortedReducedEstimatorValues;
  std::vector<std::shared_ptr<cluster>> consideredClusters = getClustersForEstimator();

  if(unsortedReducedEstimatorValuesOnClusters->size() <= j)
    return *unsortedReducedEstimatorValuesOnClusters;

  unsigned int smallestEstimatorValueIndex;

  double summaricWeight = 0.0;

  switch (positionalSecondGradeEstimatorCountingMethod)
  {
    case WEIGHTED:

      while(summaricWeight < j)
      {
        smallestEstimatorValueIndex =
          findSmallestEstimatorValueIndex(unsortedReducedEstimatorValuesOnClusters);
        summaricWeight += consideredClusters[smallestEstimatorValueIndex]->getWeight();
        jSortedReducedEstimatorValues.push_back(
          unsortedReducedEstimatorValuesOnClusters->at(smallestEstimatorValueIndex)
        );
        unsortedReducedEstimatorValuesOnClusters->erase(
          unsortedReducedEstimatorValuesOnClusters->begin() + smallestEstimatorValueIndex
        );
        consideredClusters.erase(
          consideredClusters.begin() + smallestEstimatorValueIndex
        );
      }

    break;
    case STANDARD:
    default:
      while(jSortedReducedEstimatorValues.size() < floor(j))
      {
        smallestEstimatorValueIndex =
          findSmallestEstimatorValueIndex(unsortedReducedEstimatorValuesOnClusters);
        jSortedReducedEstimatorValues.push_back(
          unsortedReducedEstimatorValuesOnClusters->at(smallestEstimatorValueIndex)
        );
        unsortedReducedEstimatorValuesOnClusters->erase(
          unsortedReducedEstimatorValuesOnClusters->begin() + smallestEstimatorValueIndex
        );
      }
    break;
  }

  while(jSortedReducedEstimatorValues.size() < j)
  {
    smallestEstimatorValueIndex =
      findSmallestEstimatorValueIndex(unsortedReducedEstimatorValuesOnClusters);
    jSortedReducedEstimatorValues.push_back(
      unsortedReducedEstimatorValuesOnClusters->at(smallestEstimatorValueIndex)
    );
    unsortedReducedEstimatorValuesOnClusters->erase(
      unsortedReducedEstimatorValuesOnClusters->begin() + smallestEstimatorValueIndex
    );
  }

  return jSortedReducedEstimatorValues;
}

unsigned int MainWindow::findSmallestEstimatorValueIndex(std::vector<double> *unsortedReducedEstimatorValuesOnClusters)
{
  double result = 2;

  for(double value : *unsortedReducedEstimatorValuesOnClusters)
    if(value < result) result = value;

  return result;
}

int MainWindow::updateClustersTemporalDerivativeTimesInARow()
{
  for(std::shared_ptr<cluster> c : clusters)
  {
    if(hasPositiveTemporalDerivative(c))
      c->positiveTemporalDerivativeTimesInARow += 1;
    else
      c->positiveTemporalDerivativeTimesInARow = 0;
  }

  for(std::vector<std::shared_ptr<cluster>> level: storedMedoids)
  {
    for(std::shared_ptr<cluster> c : level)
    {
      if(hasPositiveTemporalDerivative(c))
        c->positiveTemporalDerivativeTimesInARow += 1;
      else
        c->positiveTemporalDerivativeTimesInARow = 0;
    }
  }

  return 1;
}

bool MainWindow::hasPositiveTemporalDerivative(std::shared_ptr<cluster> c)
{
  if(KDETemporalDerivativeY.size() == 0) return true;

  unsigned int clusterPositionIndex = findClusterPositionIndex(c);

  if(clusterPositionIndex == 0)
    return KDETemporalDerivativeY[0] > 0;
  if(clusterPositionIndex == KDETemporalDerivativeY.size()-1)
    return KDETemporalDerivativeY[KDETemporalDerivativeY.size() - 1] > 0;

  return (KDETemporalDerivativeY[clusterPositionIndex] + KDETemporalDerivativeY[clusterPositionIndex + 1]) > 0;
}

unsigned int MainWindow::findClusterPositionIndex(std::shared_ptr<cluster> c)
{
  double clusterPosition =
      std::stod(c->getRepresentative()->attributesValues["Val0"]);

  unsigned int result = 0;

  while((double)(result * ui->lineEdit_domainDensity->text().toDouble()) < clusterPosition ) ++ result;

  return result;
}

int MainWindow::removeUnpromissingClusters()
{
  for(int index = clusters.size() - 1; index >= 0; --index)
  {
    // Remove cluster if its temporal derivative is not equal to 0 and
    // it's weight is insignificant
    if(clusters[index]->getWeight() < positionalSecondGradeEstimator &&
       clusters[index]->positiveTemporalDerivativeTimesInARow == 0)
    {
      clusters.erase(clusters.begin() + index);
    }
  }

  for(std::vector<std::shared_ptr<cluster>> level: storedMedoids)
  {
    for(int index = level.size() - 1; index >= 0; --index)
    {
      // Remove cluster if its temporal derivative is not equal to 0 and
      // it's weight is insignificant
      if(level[index]->getWeight() < positionalSecondGradeEstimator &&
         level[index]->positiveTemporalDerivativeTimesInARow == 0)
      {
        level.erase(level.begin() + index);
      }
    }
  }
}

void MainWindow::addLatestTemporalVelocityDensityProfilePlot()
{

  std::vector<long> keys;

  for(auto kv : temporalVelocityDensityProfile)
    keys.push_back(kv.first);

  std::vector<long>::iterator lastCountedTimestamp =
      std::max_element(keys.begin(), keys.end());

  QVector<qreal> X, Y;

  if(!temporalVelocityDensityProfile.empty())
  {
    for(auto kv : temporalVelocityDensityProfile[*lastCountedTimestamp])
    {
      X.push_back(kv.first.back());
      Y.push_back(kv.second * 1e18);
    }

    temporalVelocityDensityProfile.end();
  }

  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(X, Y);
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(Qt::green));
}

void MainWindow::addTemporalDerivativePlot(const QVector<qreal> *X, const QVector<qreal> *Y)
{
  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, *Y);
  ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(Qt::cyan));
}

void MainWindow::fillStandardDeviations(QVector<std::shared_ptr<QVector<qreal>>> *stDevs)
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
                ((QLineEdit*)(
                    ((QTableWidget*)(ui->tableWidget_targetFunctions->cellWidget(functionIndex, STDEV_COLUMN_INDEX)))
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
                ((QLineEdit*)(
                    ((QTableWidget*)(ui->tableWidget_targetFunctions->cellWidget(functionIndex, MEAN_COLUMN_INDEX)))
                    ->cellWidget(dimensionIndex, 0)
                ))
                ->text().toDouble()
            );
        }
    }
}

void MainWindow::on_pushButton_generate_clicked()
{
    // Log that application started generating KDE
    qDebug() << "KDE generation started.";
    qDebug() << "Seed: " + ui->lineEdit_seed->text() +
                ", Sample size: " + ui->lineEdit_sampleSize->text();

    srand(ui->lineEdit_seed->text().toDouble());

    int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

    QVector<std::shared_ptr<QVector<qreal>>> means, stDevs;

    fillMeans(&means);
    fillStandardDeviations(&stDevs);

    qDebug() << "Means and stDevs filled.";

    // Generate samples
    generateSamples(&means, &stDevs);

    qDebug() << "Samples generated.";

    std::shared_ptr<function> targetFunction(generateTargetFunction(&means, &stDevs));

    kernelDensityEstimator* estimator
        = generateKernelDensityEstimator(dimensionsNumber);

    qDebug() << "Testing...";

    // Test estimator
    testKDE(estimator, targetFunction.get());

    qDebug() << samples.size();

    // Run plot related tasks if dimension number is equal to 1
    if(dimensionsNumber == 1) drawPlots(estimator, targetFunction.get());
}

void MainWindow::fillDomain(QVector<std::shared_ptr<point>>* domain, std::shared_ptr<point> *prototypePoint)
{
    // Check if domain is nullpointer
    if(domain == NULL) return;

    std::shared_ptr<point> pPoint;

   // Check if prototype is a null pointer
    if(prototypePoint == NULL)
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

void MainWindow::generateSamples(QVector<std::shared_ptr<QVector<qreal>> > *means,
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

    dataReader *reader = new progressiveDistributionDataReader(targetDistribution.get(), progressionSize);

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

distribution* MainWindow::generateTargetDistribution(QVector<std::shared_ptr<QVector<qreal>>> *means,
                                                     QVector<std::shared_ptr<QVector<qreal>>> *stDevs)
{
    qreal seed = ui->lineEdit_seed->text().toDouble();

    QVector<qreal> contributions;
    QVector<std::shared_ptr<distribution>> elementalDistributions;

    int targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        contributions.append
        (
            ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(functionIndex, CONTRIBUTION_COLUMN_INDEX)))
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
        stepsNumber = ui->lineEdit_iterationsNumber->text().toDouble(),
        samplingAlgorithmID = ui->comboBox_samplingAlgorithm->currentIndex();

    switch(samplingAlgorithmID)
    {
        case BIASED_RESERVOIR_SAMPLING_ALGORITHM:
            return new biasedReservoirSamplingAlgorithm(reader, parser, sampleSize, stepsNumber);
        break;
        case BASIC_RESERVOIR_SAMPLING_ALGORITHM:
        default:
            return new basicReservoirSamplingAlgorithm(reader, parser, sampleSize, stepsNumber);
    }
}

kernelDensityEstimator* MainWindow::generateKernelDensityEstimator(int dimensionsNumber)
{
    QVector<int> kernelsIDs;
    QVector<qreal> smoothingParameters;
    QVector<QString> carriersRestrictions;

    for(int rowNumber = 0; rowNumber < dimensionsNumber; ++rowNumber)
    {
        kernelsIDs.append(((QComboBox*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, KERNEL_COLUMN_INDEX)))->currentIndex());
        smoothingParameters.append(((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->text().toDouble());
        carriersRestrictions.append(((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX)))->text());
    }

    return new kernelDensityEstimator(
                    &samples,
                    &smoothingParameters,
                    &carriersRestrictions,
                    PRODUCT,
                    &kernelsIDs
                );
}

function* MainWindow::generateTargetFunction(QVector<std::shared_ptr<QVector<qreal>>>* means,
                                             QVector<std::shared_ptr<QVector<qreal>>>* stDevs)

{
  QVector<qreal> contributions;
  QVector<std::shared_ptr<function>> elementalFunctions;

  int targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

  // Check if contributions are set correctly. If they are, then last contribution is >= 0;
  if(((QLineEdit*)(ui
                   ->tableWidget_targetFunctions
                   ->cellWidget(targetFunctionElementsNumber -1, CONTRIBUTION_COLUMN_INDEX))
                  )
          ->text().toDouble() <= 0)
  {
      // If not then uniform distributions and log error
      qDebug() << "Contributions aren't set correctly. Uniforming contributions.";
      uniformContributions();
  }



  for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
  {

      contributions.append
      (
          ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(functionIndex, CONTRIBUTION_COLUMN_INDEX)))
          ->text().toDouble()
      );


      elementalFunctions.append(std::shared_ptr<function>(new multivariateNormalProbabilityDensityFunction(means->at(functionIndex).get(),
                                                                                                          stDevs->at(functionIndex).get())));
  }

  return new complexFunction(&contributions, &elementalFunctions);
  return NULL;
}

QColor MainWindow::getRandomColor()
{
    return QColor(rand()%110 + 50, rand()%110 + 50, rand()%110 + 50);
}

void MainWindow::testKDE(kernelDensityEstimator *KDE, function *targetFunction)
{
  testKDEError(KDE, targetFunction);
  testRareElementsDetector(KDE);
}

int MainWindow::testKDEError(kernelDensityEstimator *KDE, function *targetFunction)
{
  QVector<point *> testDomain;
  fillTestDomain(&testDomain, NULL);

  qreal error = 0;

  foreach(auto arg, testDomain)
  {
    qDebug()  << "Point: " << *arg
              << "Target: " << targetFunction->getValue(arg)
              << "Estimated: " << KDE->getValue(arg)
              << "Difference: " << qAbs(targetFunction->getValue(arg) - KDE->getValue(arg));
    error += qAbs(targetFunction->getValue(arg) - KDE->getValue(arg));
  }

  qDebug() << "Error: " << error;
  qDebug() << "Average error: " << error/testDomain.size();

  return 0;
}

int MainWindow::testRareElementsDetector(kernelDensityEstimator *KDE)
{
  QVector<point *> testDomain;
  fillTestDomain(&testDomain, NULL);

  qreal r = ui->lineEdit_rarity->text().toDouble();

  QVector<int> atypicalElementsIndexes;

  rareElementsDetector* detector = new rareElementsDetector(KDE, r);

  detector->findAtypicalElementsInDomain(&testDomain, &atypicalElementsIndexes);

  qDebug() << "Atypical elements indexes:";
  qDebug() << atypicalElementsIndexes;

  return 0;
}

void MainWindow::fillTestDomain(QVector<point *> *domain, point *prototypePoint)
{
    // Check if domain is nullpointer
    if(domain == NULL) return;

   // Check if prototype is a null pointer
    if(prototypePoint == NULL)
    {
        // If so make it a point pointer
        prototypePoint = new point();
    }

    for(int i = -1; i <= 1; ++i)
    {
        prototypePoint->append(i);

        if(prototypePoint->size() == ui->spinBox_dimensionsNumber->value())
        {
            domain->append(new point());

            foreach(qreal dimensionVal, *prototypePoint) domain->last()->append(dimensionVal);
        }
        else fillTestDomain(domain, prototypePoint);

        prototypePoint->removeLast();
    }
}

void MainWindow::on_pushButton_animate_clicked()
{
    int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

    if(!canAnimationBePerformed(dimensionsNumber)) return;

    // Log that application started generating KDE
    qDebug() << "KDE animation started.";
    qDebug() << "Seed: " + ui->lineEdit_seed->text() +
                ", Sample size: " + ui->lineEdit_sampleSize->text();

    srand(ui->lineEdit_seed->text().toDouble());

    fillMeans(&means);
    fillStandardDeviations(&stDevs);

    std::shared_ptr<function> targetFunction(generateTargetFunction(&means, &stDevs));

    std::shared_ptr<kernelDensityEstimator> estimator(generateKernelDensityEstimator(dimensionsNumber));
    kernelPrognoser.reset(generateKernelDensityEstimator(dimensionsNumber));

    std::shared_ptr<distribution> targetDistribution(generateTargetDistribution(&means, &stDevs));

    qDebug() << "Adding alternative distribution.";

    //*** Alternative distribution created solely for presentation *************
    QVector<std::shared_ptr<QVector<qreal>>> alternativeMeans, alternativeStDevs;

    std::shared_ptr<QVector<qreal>> alternativeMean;
    alternativeMean.reset(new QVector<qreal>());
    alternativeMean->push_back(10);

    alternativeMeans.push_back(alternativeMean);
    alternativeStDevs = stDevs;

    std::shared_ptr<distribution>
        alternativeDistribution(generateTargetDistribution(&alternativeMeans,
                                                           &alternativeStDevs));
    //**************************************************************************

    parser.reset(new distributionDataParser(&attributesData));

    qreal progressionSize = ui->lineEdit_distributionProgression->text().toDouble();

    reader.reset(new progressiveDistributionDataReader(targetDistribution.get(), progressionSize));

    reader->gatherAttributesData(&attributesData);
    parser->setAttributesOrder(reader->getAttributesOrder());

    positionalSecondGradeEstimatorCountingMethod = ui->comboBox_rareElementsMethod->currentIndex();

    reservoirSamplingAlgorithm* algorithm = generateReservoirSamplingAlgorithm(reader.get(), parser.get());

    objects.clear();
    storage.clear();

    int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();

    groupingThread gt(&storedMedoids, parser);

    gt.setAttributesData(&attributesData);

    qDebug() << "Attributes data set.";

    gt.initialize();

    double weightUpdateCoefficient =
        ui->lineEdit_weightModifier->text().toDouble();

    for(stepNumber = 0; stepNumber < stepsNumber; ++stepNumber)
    {
      updateWeights();
      //storage.updateWeights(weightUpdateCoefficient);

      qDebug() << "Performing a step";

      if(stepNumber == 290)
      {
        ((progressiveDistributionDataReader*)(reader.get()))->setNewSource(alternativeDistribution.get());
        //reader.reset(new progressiveDistributionDataReader(alternativeDistribution.get(), 0));
        //reader.reset(new progressiveDistributionDataReader(targetDistribution.get(), progressionSize));
        //algorithm = generateReservoirSamplingAlgorithm(reader.get(), parser.get());
      }

      qDebug() << "After if";

      algorithm->performSingleStep(&objects, stepNumber);

      // TR TODO: It's not working for biased algorithm
      std::shared_ptr<cluster> newCluster =
          std::shared_ptr<cluster>(new cluster(stepNumber, objects.back()));
      newCluster->setTimestamp(stepNumber);

      clusters.push_back(newCluster);
      storage.addCluster(newCluster, 0);

      qDebug() << "Reservoir size in step " << stepNumber
               << " is: " << clusters.size();
      qDebug() << "Storage size in step " << stepNumber
               << " is: " << storage.size();

      if(clusters.size() >= algorithm->getReservoidMaxSize())
      {
        std::vector<double> unsortedReducedEstimatorValuesOnClusters
            = countUnsortedReducedEstimatorValuesOnEstimatorClusters(estimator.get());

        positionalSecondGradeEstimator =
          countPositionalSecondGradeEstimator(&unsortedReducedEstimatorValuesOnClusters);

        findUncommonClusters(estimator.get());

        removeUnpromissingClusters();

        qDebug() << "Clusters size after reduction: " << clusters.size();

        if(clusters.size() < MEDOIDS_NUMBER) continue;

        std::vector<std::shared_ptr<cluster>> currentClusters
            = getClustersForEstimator();

        std::shared_ptr<weightedSilvermanSmoothingParameterCounter>
          smoothingParamCounter( new weightedSilvermanSmoothingParameterCounter(&currentClusters, 0));

        std::vector<double> smoothingParameters;

        smoothingParameters
          .push_back(smoothingParamCounter->countSmoothingParameterValue());

        estimator->setSmoothingParameters(smoothingParameters);

        std::shared_ptr<groupingThread> gThread(
              new groupingThread(&storedMedoids, parser)
        );

        //gThread->setAttributesData(&attributesData);
        //gThread->getClustersForGrouping(clusters);

        runningSubthreads.push_back(gThread);

        qDebug() << "Got objects for grouping.";

        //runningSubthreads.back()->start();

        gt.getClustersForGrouping(clusters);
        gt.run();

        objects.clear();
        clusters.clear();
        clustersForVDE.clear();

        qDebug() << "Objects cleared.";
      }

      updateClustersTemporalDerivativeTimesInARow();

      std::vector<std::shared_ptr<cluster>> currentClusters
          = getClustersForEstimator();

      estimator->setClusters(currentClusters);

      targetFunction.reset(generateTargetFunction(&means, &stDevs));

      drawPlots(estimator.get(), targetFunction.get());

      start = std::chrono::duration_cast< std::chrono::milliseconds >(
          std::chrono::system_clock::now().time_since_epoch()).count();

      delay(ui->lineEdit_milisecondsDelay->text().toInt());
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

  return 0;
}

void MainWindow::delay( int ms )
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

    // TODO TR: Ensure that this doesn't result in memory leaks
    ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, KERNEL_COLUMN_INDEX, new QComboBox());

    ((QComboBox*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, KERNEL_COLUMN_INDEX)))->insertItems(0, kernelTypes);

    // Add input box with validator for smoothing parameters
    ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX, new QLineEdit());

    // TODO TR: Ensure that this doesn't result in memory leaks
    ((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->setText("1.0");
    ((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->setValidator(smoothingParameterValidator);

    // Add input box for carrier restriction value
    ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX, new QLineEdit());

    // TODO TR: Ensure that this doesn't result in memory leaks
    ((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX)))->setText("None.");
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

    QTableWidget *targetFunctionTablePointer = (QTableWidget*)ui->tableWidget_targetFunctions,
                 *meansTablePointer, *stDevsTablePointer;

    for(int rowIndex = 0; rowIndex < numberOfRows; ++rowIndex)
    {
        // TODO TR: Ensure that this doesn't result in memory leaks
        targetFunctionTablePointer->setCellWidget(rowIndex, MEAN_COLUMN_INDEX, new QTableWidget());

        meansTablePointer = ((QTableWidget*)ui->tableWidget_targetFunctions->cellWidget(rowIndex, MEAN_COLUMN_INDEX));
        meansTablePointer->setRowCount(dimensionsNumber);
        meansTablePointer->setColumnCount(1);
        meansTablePointer->horizontalHeader()->hide();

        // TODO TR: Ensure that this doesn't result in memory leaks
        targetFunctionTablePointer->setCellWidget(rowIndex, STDEV_COLUMN_INDEX, new QTableWidget());

        stDevsTablePointer = (QTableWidget*)ui->tableWidget_targetFunctions->cellWidget(rowIndex, STDEV_COLUMN_INDEX);
        stDevsTablePointer->setRowCount(dimensionsNumber);
        stDevsTablePointer->setColumnCount(1);
        stDevsTablePointer->horizontalHeader()->hide();

        for(int dimensionNumber = 0; dimensionNumber < dimensionsNumber; ++dimensionNumber)
        {
           meansTablePointer->setCellWidget(dimensionNumber, 0, new QLineEdit());
           ((QLineEdit*)(meansTablePointer->cellWidget(dimensionNumber, 0)))->setText("0.0");
           ((QLineEdit*)(meansTablePointer->cellWidget(dimensionNumber, 0)))->setValidator(meanValidator);

           stDevsTablePointer->setCellWidget(dimensionNumber, 0, new QLineEdit());
           ((QLineEdit*)(stDevsTablePointer->cellWidget(dimensionNumber, 0)))->setText("1.0");
           ((QLineEdit*)(stDevsTablePointer->cellWidget(dimensionNumber, 0)))->setValidator(stDevValidator);
        }

        // TODO TR: Ensure that this doesn't result in memory leaks
        targetFunctionTablePointer->setCellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX, new QLineEdit());
        ((QLineEdit*)(targetFunctionTablePointer->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->setMaxLength(6);
        ((QLineEdit*)(targetFunctionTablePointer->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->setValidator(contributionValidator);
        QObject::connect(((QLineEdit*)(targetFunctionTablePointer->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX))), SIGNAL(textEdited(QString)), this, SLOT(updateLastContribution()));
    }

    // Disable last contribution cell, as it's filled automatically
    ((QLineEdit*)(targetFunctionTablePointer->cellWidget(numberOfRows -1, CONTRIBUTION_COLUMN_INDEX)))->setEnabled(false);

    uniformContributions();
}

void MainWindow::uniformContributions()
{
    int numberOfRows = ui->tableWidget_targetFunctions->rowCount(), lastRowIndex = numberOfRows - 1;

    for(int rowIndex = 0; rowIndex < lastRowIndex; ++rowIndex)
        ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->setText(QString::number(100.0/numberOfRows));

    ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(lastRowIndex, CONTRIBUTION_COLUMN_INDEX)))->setText(QString::number(countLastContribution()));
}

qreal MainWindow::countLastContribution()
{
    qreal result = 100.0;

    int lastRowIndex = ui->tableWidget_targetFunctions->rowCount()-1;

    for(int rowIndex = 0; rowIndex < lastRowIndex; ++rowIndex)
        result -= ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->text().toDouble();

    return result;
}

void MainWindow::updateLastContribution()
{
    int lastRowIndex = ui->tableWidget_targetFunctions->rowCount()-1;
    qreal lastContributionValue = countLastContribution();

    ((QLineEdit*)(ui
                  ->tableWidget_targetFunctions
                  ->cellWidget(lastRowIndex, CONTRIBUTION_COLUMN_INDEX)))
                  ->setText(QString::number(lastContributionValue));
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

    foreach(std::shared_ptr<QVector<qreal>> sample, samples)
        samplesColumn.append(sample.get()->at(rowNumber));

    //value = (counter.*methodFunctionPointer)();
    value = counter->countSmoothingParameterValue();

    ((QLineEdit*)(ui
                  ->tableWidget_dimensionKernels
                  ->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))
      ->setText(QString::number(value));
  }
}

smoothingParameterCounter *MainWindow::generateSmoothingParameterCounter(QVector<qreal> *samplesColumn)
{
  int smoothingParameterCounterID = ui
                                    ->comboBox_smoothingParameterCountingMethod
                                    ->currentIndex();

  switch (smoothingParameterCounterID)
  {
    case WEIGHTED_SILVERMAN:
      return new weightedSilvermanSmoothingParameterCounter(samplesColumn, NULL);
    break;
    case RANK_3_PLUG_IN:
      return new pluginSmoothingParameterCounter(samplesColumn, 3);
    break;
    case RANK_2_PLUG_IN:
    default:
      return new pluginSmoothingParameterCounter(samplesColumn, 2);
    break;
  }

  return NULL;
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
  double weightModifier = ui->lineEdit_weightModifier->text().toDouble();
  double weightDeletionThreshold = ui->lineEdit_deletionThreshold->text().toDouble();

  for(unsigned int level = 0; level < storedMedoids.size(); ++level)
  {
    for(unsigned int medoidNumber = 0; medoidNumber < storedMedoids[level].size(); ++medoidNumber)
    {
      storedMedoids[level][medoidNumber]->setWeight(
        weightModifier * storedMedoids[level][medoidNumber]->getWeight()
      );

      /*
      if(storedMedoids[level][medoidNumber]->getWeight() < weightDeletionThreshold)
        storedMedoids[level].erase(
          storedMedoids[level].begin() + medoidNumber
        );
      */
    }
  }

  for(unsigned int i = 0; i < clusters.size(); ++i)
  {
    clusters[i]->setWeight(weightModifier * clusters[i]->getWeight());
    if(clusters[i]->getWeight() < weightDeletionThreshold) clusters.erase(clusters.begin() + i);
  }

  qDebug() << "Weights updated.";
}

void MainWindow::moveTargetFunctionLeft()
{
  //QVector<std::shared_ptr<QVector<qreal>>> newMeans;
  //newMeans.push_back(std::shared_ptr<QVector<qreal>>(new QVector<qreal>({-1})));
  //means = newMeans;
  //means = QVector<std::shared_ptr<QVector<qreal>>>(std::shared_ptr<QVector<qreal>>());
}
