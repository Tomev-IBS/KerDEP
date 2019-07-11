#include "DESDA.h"

#include <QTime>
#include <QCoreApplication>
#include <algorithm>
#include <fstream>
#include <math.h>

#include "Reservoir_sampling/distributionDataSample.h"

DESDA::DESDA(std::shared_ptr<kernelDensityEstimator> estimator,
             std::shared_ptr<kernelDensityEstimator> estimatorDerivative,
             std::shared_ptr<kernelDensityEstimator> enchancedKDE,
             double weightModifier,
             weightedSilvermanSmoothingParameterCounter *smoothingParamCounter,
             reservoirSamplingAlgorithm *samplingAlgorithm,
             std::vector<std::shared_ptr<cluster> > *clusters,
             std::vector<std::vector<std::shared_ptr<cluster> > > *storedMedoids,
             double desiredRarity, groupingThread *gt, double v,
             double newWeightB, int mE, int kpssX, int lambda):
  _weightModifier(weightModifier), _samplingAlgorithm(samplingAlgorithm),
  _estimatorDerivative(estimatorDerivative), _estimator(estimator),
  _smoothingParamCounter(smoothingParamCounter), _clusters(clusters),
  _storedMedoids(storedMedoids), _desiredRarity(desiredRarity),
  _grpThread(gt), _v(v), _newWeightB(newWeightB),
  _enhancedKDE(enchancedKDE), _mE(mE), _lambda(lambda)
{
  _objects.clear();
  std::shared_ptr<sample> e1000Sample =
      std::make_shared<distributionDataSample>();
  emE = cluster(e1000Sample);
  emE._deactualizationParameter = w_E;

  _maxM = 2 * mE;

  _samplingAlgorithm->changeReservoirMaxSize(_maxM);

  _mE = 1000;
  _m = _maxM;

  _kpssM = 250;// 2;
  int l = kpssX * pow(_kpssM / 100, 0.25);

  _stepNumber = 1;

  stationarityTest.reset(new KPSSStationarityTest(_kpssM, avg, l));  
}

int round(int val){
  return ceil(val) - val >= 0.5 ? ceil(val) : floor(val);
}

int sgn(double val)
{
  return (val > 0) - (val < 0);
}

double sum(std::vector<double> vals)
{
  double valsSum = 0.0;

  for(auto val : vals) valsSum += val;

  return valsSum;
}

double absSum(std::vector<double> vals)
{
  double valsSum = 0.0;

  for(auto val : vals) valsSum += fabs(val);

  return valsSum;
}

double stDev(std::vector<double> vals)
{
  int m = vals.size();

  if(m < 2) return 1.0;

  double sigma = 0.0;

  sigma = sum(vals);
  sigma *= sigma / m;

  double sqrSum = 0;

  for(auto val : vals) sqrSum += val * val;

  sigma = sqrSum - sigma;

  sigma /= m - 1;

  return std::sqrt(sigma);
}

void DESDA::performStep()
{
  updateM();

  // If weights degrades geomatrically
  if(_shouldCluster) updateWeights();

  _samplingAlgorithm->performSingleStep(&_objects, _stepNumber);

  std::shared_ptr<cluster> newCluster =
      std::shared_ptr<cluster>(new cluster(_stepNumber, _objects.back()));
  newCluster->setTimestamp(_stepNumber);

  while(_clusters->size() >= _maxM && !_shouldCluster)
  {
    //_clusters->erase(_clusters->begin(), _clusters->begin()+1);
    _clusters->pop_back();
    _objects.erase(_objects.begin(), _objects.begin() + 1);
  }

  //_clusters->push_back(newCluster);
  _clusters->insert(_clusters->begin(), newCluster);

  // If weights degrades accodring to new formula
  if(!_shouldCluster) updateWeights();

  _smoothingParamCounter->updateSmoothingParameterValue(
    _weightModifier,
    std::stod(_clusters->front()->getObject()->attributesValues["Val0"])
  );

  _smoothingParamCounter->setClusters(_clusters, 0);

  std::vector<double> smoothingParameters =
  {
    _smoothingParamCounter->countSmoothingParameterValue() * _smoothingParameterMultiplier
  };

  _estimator->setSmoothingParameters(smoothingParameters);
  _estimatorDerivative->setSmoothingParameters(smoothingParameters);
  //smoothingParameters[0] = smoothingParameters[0];

  _enhancedKDE->setSmoothingParameters(smoothingParameters);

  auto currentClusters = getClustersForEstimator();

  qDebug() << "Reservoir size in step "
             << _stepNumber << " is: " << currentClusters.size();

  //qDebug () << "1";

  countKDEValuesOnClusters();

  //qDebug () << "2";

  //avg = getAverageOfFirstMSampleValues(_mE);
  avg = getNewEmEValue();

  /*
  if(_stepNumber % 10 == 0){
      saveWeightsToFile(std::to_string(_stepNumber) + ".txt");
      saveEmEWeightsToFile("EmE_" + std::to_string(_stepNumber) + ".txt");
  }
  */

  //qDebug () << "3";

  // Start at 0
  //if(_stepNumber >= _mE)
  stationarityTest->addNewSample(
      std::stod(_clusters->front()->getObject()->attributesValues["Val0"])
  );

  //qDebug () << "4";

  emE._currentKDEValue = avg;

  if(currentClusters.size() >= _samplingAlgorithm->getReservoidMaxSize()
     && _shouldCluster)
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
    // Reassigment needed, as (probably) grouping algorithm resets it
    _clusters = &((*_storedMedoids)[0]);
  }

  //updateA();

  //qDebug () << "Estimator";

  _estimator->setClusters(currentClusters);

  updatePrognosisParameters();

  if(emE.predictionParameters.size() == 0)
    emE.initializePredictionParameters(avg);
  else
    emE.updatePredictionParameters(avg);

  // Save to file
  /*
  std::string rowToSave =
    _clusters->front()->getObject()->attributesValues["Val0"] + ",";
  rowToSave += e1000.rowToSave;
  qDebug() << QString::fromStdString(rowToSave);

  std::ofstream experimentDataFile;

  experimentDataFile.open("d:\\Dysk Google\\Badania\\experimentData.csv", std::ios_base::app);
  experimentDataFile << rowToSave;

  experimentDataFile.close();
  */

  while(aemEVals.size() >= _mE)
  {
    emEVals.pop_back();
    aemEVals.pop_back();
  }

  emEVals.insert(emEVals.begin(), emE._currentKDEValue);
  aemEVals.insert(aemEVals.begin(), emE.predictionParameters[1]);

  countKDEDerivativeValuesOnClusters();

  ++_stepNumber;
}

void DESDA::updateWeights()
{
  if(_shouldCluster)
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
  else
  {
    for(int clusterNum = _clusters->size() - 1; clusterNum > -1; --clusterNum)
    {
      // In formula it's (i - 1), but indexes are from 1 not 0, thus no -1.
      double newWeight =
          1.0 - _newWeightB * (clusterNum) / _m;

      newWeight = std::max(0.0, newWeight);

      (*_clusters)[clusterNum]->setWeight(newWeight);
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

  if(!_shouldCluster)
  {
    while(consideredClusters.size() > _m){
      //qDebug() << "m = " << _m;
      consideredClusters.pop_back();
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
      //c->updateDeactualizationParameter(c->_currentKDEValue);
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

void DESDA::updateM()
{
  int m = 0;

  if(emE.predictionParameters.size() < 2) return;

  m = round(_maxM * ( 1 - _lambda * _u_i * fabs(emE.predictionParameters[1]) / getStdDevOfFirstMSampleValues(_mE)));

  m = std::max(m, _minM);
  _m = m;
}

double DESDA::getNewEmEValue()
{
    // Initialize EmEWeights vector
    _EmEWeights.clear();
    _EmEWeightsSum = 0;

    for(int i = 0; i < _mE; ++i){
        double val = 1.0 - _newWeightB * i / _mE;
        _EmEWeights.push_back(val);
        _EmEWeightsSum += val;
    }

    // Get weighted average
    double avg = 0;
    double mE = std::min((int)_clusters->size(), _mE);
    for(int i = 0; i < mE; ++i){
        avg += std::stod(_clusters->at(i)->getObject()->attributesValues["Val0"]) * _EmEWeights[i];
    }

    return avg / _EmEWeightsSum;
}

void DESDA::saveWeightsToFile(std::string fileName)
{


    std::string filePath = "D:\\Dysk Google\\TR Badania\\TEST reEksperyment 196 (w_0.98, m_E=m_Eta=500)\\" + fileName;
    std::string line = "";

    std::ofstream experimentDataFile;
    experimentDataFile.open(filePath, std::ios_base::app);

    experimentDataFile << "b = " << std::to_string(_newWeightB) << ", m = " << std::to_string(_m) << std::endl;

    for(int i = 0; i < _clusters->size(); ++i){
         line = std::to_string(i) + ". " + std::to_string((*_clusters)[i]->getWeight()) + "\n";
         experimentDataFile << line;
    }

    experimentDataFile.close();
}

void DESDA::saveEmEWeightsToFile(std::string fileName)
{
    std::string filePath = "D:\\Dysk Google\\TR Badania\\TEST reEksperyment 196 (w_0.98, m_E=m_Eta=500)\\" + fileName;
    std::string line = "";

    std::ofstream experimentDataFile;
    experimentDataFile.open(filePath, std::ios_base::app);

    experimentDataFile << "b = " << std::to_string(_newWeightB) << ", m = " << std::to_string(_m) << std::endl;

    for(int i = 0; i < _EmEWeights.size(); ++i){
         line = std::to_string(i) + ". " + std::to_string(_EmEWeights[i]) + "\n";
         experimentDataFile << line;
    }

    experimentDataFile.close();
}

QVector<double> DESDA::getKernelPrognosisDerivativeValues(const QVector<qreal> *X)
{
  std::vector<std::shared_ptr<cluster>> currentClusters
      = getClustersForEstimator();

  std::vector<double> prognosisCoefficients = {};
  QVector<double> kernelPrognosisDerivativeValues = {};

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
        _estimatorDerivative->getValue(&pt) * 1000 // For visibility
      );
    }
  }

  return kernelPrognosisDerivativeValues;
}

QVector<double> DESDA::getEnhancedKDEValues(const QVector<qreal> *X)
{
  //qDebug() << "EKDE values getting.";

  std::vector<std::shared_ptr<cluster>> currentClusters
      = getClustersForEstimator();

  std::vector<double> standardWeights = {};
  QVector<double> enhancedKDEValues = {};

  QVector<qreal> clustersXs = {};

  for(auto c : currentClusters){
    clustersXs.push_back(
      std::stod(c->getRepresentative()->attributesValues["Val0"])
    );
  }

  QVector<double> derivativeVal = getKernelPrognosisDerivativeValues(&clustersXs);

  _selectedVValues.clear();

  // Enhance weights of clusters
  double enhancedWeight = 0.0;
  double v_i = 0.0;

  _u_i = 0.0;

  // Count u_i
  if(_stepNumber >= 1000)
    _u_i = 1.0 / (1 + exp(- (_alpha * fabs(getStationarityTestValue()) - _beta)));

  _newWeightB = _u_i;

  double avgC2 = 0;
  double maxAParam = 0;
  //for(auto c : currentClusters)
  for(int i = 0; i < currentClusters.size(); ++i)
  {
    std::shared_ptr<cluster> c = currentClusters[i];
    enhancedWeight = c->getWeight();

    // Count v_i
    v_i = 2 * delta * ( 1.0 / (1.0 + exp(- gamma * derivativeVal[i])) - 0.5);

    maxAParam = std::max(derivativeVal[i] * 1000, maxAParam);

    enhancedWeight *= (1 + _u_i * v_i);

    standardWeights.push_back(c->getWeight());

    if(standardWeights.size() == 10  || standardWeights.size() == 50  ||
       standardWeights.size() == 200 /*|| standardWeights.size() == 300 ||
       standardWeights.size() == 500 || standardWeights.size() == 700*/)
      _selectedVValues.push_back(v_i);

    avgC2 += c->predictionParameters[1];
    c->setWeight(enhancedWeight);
  }

  avgC2 /= currentClusters.size();
  //qDebug() << "avgC2 = " << avgC2;

  _enhancedKDE->setClusters(currentClusters);

  for(qreal x: *X)
  {
    QVector<qreal> pt;
    pt.push_back(x);
    enhancedKDEValues.push_back(
      _enhancedKDE->getValue(&pt)
    );
  }

  // Restore weights
  for(unsigned int i = 0; i < currentClusters.size(); ++i)
    currentClusters[i]->setWeight(standardWeights[i]);

  //qDebug() << "Max a param: " << maxAParam;
  return enhancedKDEValues;
}

double DESDA::getAverageOfFirstMSampleValues(int M)
{
  double avg = 0;

  int m0 = std::min(M, (int)_clusters->size());

  for(int i = 0; i < m0; ++i)
    avg += std::stod(_clusters->at(i)->getObject()->attributesValues["Val0"]);

  return avg / m0;
}

double DESDA::getStdDevOfFirstMSampleValues(int M)
{
  int m0 = std::min(M, (int)_clusters->size());

  if(m0 == 1) return 1;

  double sumOfVals = 0, sumOfSquares = 0, val = 0;

  for(int i = 0; i < m0; ++i){
    val = std::stod(_clusters->at(i)->getObject()->attributesValues["Val0"]);
    sumOfVals += val;
    sumOfSquares += pow(val, 2);
  }

  sumOfVals *= sumOfVals;
  sumOfVals /= m0 * (m0 - 1);

  sumOfSquares /= m0 - 1;

  return pow(sumOfSquares - sumOfVals, 0.5);
}

cluster DESDA::getEmECluster()
{
  return emE;
}

double DESDA::getStationarityTestValue()
{
  return stationarityTest->getTestsValue();
}

double DESDA::emEStDev()
{
  return stDev(emEVals);
}

double DESDA::aemEAvg()
{
  return sum(aemEVals) / aemEVals.size();
}

double DESDA::aemEStDev()
{
  return stDev(aemEVals);
}

double DESDA::aemEVersor()
{
  double aesum = sum(aemEVals);
  double aeabssum = absSum(aemEVals);

  if(_stepNumber > 999 && _stepNumber < 1010){
    qDebug() << "Step  : " << _stepNumber;
    qDebug() << "Sum   = " << aesum;
    qDebug() << "abSum = " << aeabssum;
  }

  return aesum / aeabssum;
}


