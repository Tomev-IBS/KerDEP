#include "DESDA.h"

#include <QTime>
#include <QCoreApplication>

DESDA::DESDA(std::shared_ptr<kernelDensityEstimator> estimator,
             std::shared_ptr<kernelDensityEstimator> estimatorDerivative,
             double weightModifier,
             weightedSilvermanSmoothingParameterCounter *smoothingParamCounter,
             reservoirSamplingAlgorithm *samplingAlgorithm,
             std::vector<std::shared_ptr<cluster> > *clusters,
             std::vector<std::vector<std::shared_ptr<cluster> > > *storedMedoids,
             double desiredRarity, groupingThread *gt, double v):
  _weightModifier(weightModifier), _samplingAlgorithm(samplingAlgorithm),
  _estimatorDerivative(estimatorDerivative), _estimator(estimator),
  _smoothingParamCounter(smoothingParamCounter), _clusters(clusters),
  _storedMedoids(storedMedoids), _desiredRarity(desiredRarity),
  _grpThread(gt), _v(v)
{
  _objects.clear();
}

void DESDA::performStep()
{
  updateWeights();
  _samplingAlgorithm->performSingleStep(&_objects, _stepNumber);

  std::shared_ptr<cluster> newCluster =
      std::shared_ptr<cluster>(new cluster(_stepNumber, _objects.back()));
  newCluster->setTimestamp(_stepNumber);

  _clusters->push_back(newCluster);

  _smoothingParamCounter->updateSmoothingParameterValue(
    _weightModifier,
    std::stod(_clusters->back()->getObject()->attributesValues["Val0"])
  );

  std::vector<double> smoothingParameters =
  {
    _smoothingParamCounter->getSmoothingParameterValue() * _smoothingParameterMultiplier
  };

  _estimator->setSmoothingParameters(smoothingParameters);
  _estimatorDerivative->setSmoothingParameters(smoothingParameters);

  auto currentClusters = getClustersForEstimator();

  qDebug() << "Reservoir size in step "
             << _stepNumber << " is: " << currentClusters.size();

  countKDEValuesOnClusters();

  if(currentClusters.size() >= _samplingAlgorithm->getReservoidMaxSize())
  {
    std::vector<std::shared_ptr<cluster>> clustersForGrouping;

    for(int cNum = 0; cNum < _numberOfClustersForGrouping; ++cNum)
    {
      clustersForGrouping.push_back((*_clusters)[0]);
      _clusters->erase(_clusters->begin(), _clusters->begin()+1);
    }

    _objects.erase(_objects.begin(),
                   _objects.begin() + (_numberOfClustersForGrouping - _medoidsNumber));
    _grpThread->getClustersForGrouping(clustersForGrouping);
    _grpThread->run();
    currentClusters = getClustersForEstimator();
    _clusters = &((*_storedMedoids)[0]); // Reassigment needed, as (probably) grouping algorithm resets it
  }

  //updateA();

  _estimator->setClusters(currentClusters);

  updatePrognosisParameters();
  countKDEDerivativeValuesOnClusters();

  ++_stepNumber;
}

void DESDA::updateWeights()
{
  for(unsigned int level = 0; level < _storedMedoids->size(); ++level)
  {
    for(unsigned int medoidNumber = 0; medoidNumber < (*_storedMedoids)[level].size(); ++medoidNumber)
    {
      (*_storedMedoids)[level][medoidNumber]->setWeight(
        _weightModifier * (*_storedMedoids)[level][medoidNumber]->getWeight()
      );
    }
  }
}

std::vector<std::shared_ptr<cluster> > DESDA::getClustersForEstimator()
{
  std::vector<std::shared_ptr<cluster>> consideredClusters;

  for(std::vector<std::shared_ptr<cluster>> level : (*_storedMedoids))
  {
    for(std::shared_ptr<cluster> c : level)
    {
      if(c->getWeight() >= _positionalSecondGradeEstimator)
        consideredClusters.push_back(c);
    }
  }

  return consideredClusters;
}

void DESDA::countKDEValuesOnClusters()
{
  QVector<qreal> x;

  auto consideredClusters = getClustersForEstimator();
  _estimator->setClusters(consideredClusters);

  for(std::shared_ptr<cluster> c : consideredClusters)
  {
    x.clear();
    x.push_back(std::stod(c->getRepresentative()->attributesValues["Val0"]));
    double estimatorValueOnCluster = _estimator->getValue(&x);
    c->_currentKDEValue = estimatorValueOnCluster;
  }
}

unsigned long long DESDA::findUncommonClusters()
{
  _uncommonClusters.clear();

  auto consideredClusters = getClustersForEstimator();

  for(auto c : consideredClusters)
  {
    if(c->_currentKDEValue < _maxEstimatorValueOnDomain * _a)
      _uncommonClusters.push_back(c);
  }

  return _uncommonClusters.size();
}

void DESDA::updateA()
{
  auto allConsideredClusters = getClustersForEstimator();
  double summaricClustersWeight = 0.0;

  for(std::shared_ptr<cluster> c : allConsideredClusters)
    summaricClustersWeight += c->getWeight();

  double currentUncommonClusterWeight = 0.0;

  for(std::shared_ptr<cluster> uc : _uncommonClusters)
    currentUncommonClusterWeight += uc->getWeight();

  currentUncommonClusterWeight /= summaricClustersWeight;

  _a += (_desiredRarity - currentUncommonClusterWeight)
      * _maxEstimatorValueOnDomain;

  if(_a > _MAX_A) _a = _MAX_A;
  if(_a < _MIN_A) _a = _MIN_A;

  _previousUncommonClustersWeight = currentUncommonClusterWeight;
}

void DESDA::updatePrognosisParameters()
{
  auto currentClusters = getClustersForEstimator();

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

void DESDA::countKDEDerivativeValuesOnClusters()
{
  std::vector<std::shared_ptr<cluster>> currentClusters
      = getClustersForEstimator();

  std::vector<double> prognosisCoefficients = {};

  for(auto c : currentClusters)
    prognosisCoefficients.push_back(c->predictionParameters[1]);

  if(prognosisCoefficients.size() == currentClusters.size())
  {
    _estimatorDerivative->setAdditionalMultipliers(prognosisCoefficients);
    _estimatorDerivative->setClusters(currentClusters);

    QVector<qreal> x;

    for(auto c : currentClusters)
    {
      QVector<qreal> pt;
      pt.push_back(std::stod(c->getRepresentative()->attributesValues["Val0"]));

      c->_KDEDerivativeValue = _estimatorDerivative->getValue(&pt);
    }
  }
}

QVector<double> DESDA::getKernelPrognosisDerivativeValues(const QVector<qreal> *X)
{
  std::vector<std::shared_ptr<cluster>> currentClusters
      = getClustersForEstimator();

  std::vector<double> prognosisCoefficients;

  QVector<double> kernelPrognosisDerivativeValues = {};

  prognosisCoefficients.clear();

  for(auto c : currentClusters)
    prognosisCoefficients.push_back(c->predictionParameters[1]);

  if(prognosisCoefficients.size() == currentClusters.size())
  {
    _estimatorDerivative->setAdditionalMultipliers(prognosisCoefficients);
    _estimatorDerivative->setClusters(currentClusters);

    for(qreal x: *X)
    {
      QVector<qreal> pt;
      pt.push_back(x);
      kernelPrognosisDerivativeValues.push_back(
        _estimatorDerivative->getValue(&pt) / pow(_v, 0.5)
      );
    }
  }

  return kernelPrognosisDerivativeValues;
}
