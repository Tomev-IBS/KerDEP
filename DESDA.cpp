#include "DESDA.h"

#include <QTime>
#include <QCoreApplication>
#include <QRandomGenerator>
#include <algorithm>
#include <fstream>
#include <math.h>
#include <numeric>

#include "Reservoir_sampling/distributionDataSample.h"

DESDA::DESDA(std::shared_ptr<kernelDensityEstimator> estimator,
             std::shared_ptr<kernelDensityEstimator> estimatorDerivative,
             std::shared_ptr<kernelDensityEstimator> enchancedKDE,
             double weightModifier,
             weightedSilvermanSmoothingParameterCounter *smoothingParamCounter,
             reservoirSamplingAlgorithm *samplingAlgorithm,
             std::vector<std::shared_ptr<cluster> > *clusters,
             std::vector<std::shared_ptr<cluster> > *storedMedoids,
             double desiredRarity, groupingThread *gt,
             double newWeightB, int mE, int kpssX, int lambda):
  _weightModifier(weightModifier), _samplingAlgorithm(samplingAlgorithm),
  _estimatorDerivative(estimatorDerivative), _estimator(estimator),
  _smoothingParamCounter(smoothingParamCounter), _clusters(clusters),
  _storedMedoids(storedMedoids), _r(desiredRarity),
  _grpThread(gt), _newWeightB(newWeightB),
  _enhancedKDE(enchancedKDE), _mE(mE), _lambda(lambda)
{
  _objects.clear();
  std::shared_ptr<sample> e1000Sample =
      std::make_shared<distributionDataSample>();
  emE = cluster(e1000Sample);
  emE._deactualizationParameter = w_E;

  _maxM = 2 * mE;

  _samplingAlgorithm->changeReservoirMaxSize(_maxM);
  _clustersExaminedForAsIndices = {(unsigned int)(0.2 * _maxM),
                                   (unsigned int)(0.5 * _maxM),
                                   (unsigned int)(0.8 * _maxM)};

  _mE = 1000;
  _m = _maxM;
  _minM = _maxM / 5;
  _kpssM = _maxM;
  int l = round(kpssX * pow(_kpssM / 100, 0.25));

  _sgmKPSS = -1;

  _stepNumber = 1;

  stationarityTest.reset(new KPSSStationarityTest(_kpssM, avg, l));  
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

double sigmoid(double x){
    return 1.0 / (1.0 + exp(-x));
}

void DESDA::performStep()
{
  updateM();
  updateDelta();

  _samplingAlgorithm->performSingleStep(&_objects, _stepNumber);

  std::shared_ptr<cluster> newCluster =
      std::shared_ptr<cluster>(new cluster(_stepNumber, _objects.back()));
  newCluster->setTimestamp(_stepNumber);

  while(_clusters->size() >= _maxM)
  {
    _clusters->pop_back();
    _objects.erase(_objects.begin(), _objects.begin() + 1);
  }

  _clusters->insert(_clusters->begin(), newCluster);
  updateExaminedClustersIndices();

  _smoothingParamCounter->updateSmoothingParameterValue(
    _weightModifier,
    std::stod(_clusters->front()->getObject()->attributesValues["Val0"])
  );

  _smoothingParamCounter->setClusters(_clusters, 0);

  std::vector<double> smoothingParameters = { _smoothingParamCounter->countSmoothingParameterValue() * _smoothingParameterMultiplier };

  _estimator->setSmoothingParameters(smoothingParameters);
  _estimatorDerivative->setSmoothingParameters(smoothingParameters);
  _enhancedKDE->setSmoothingParameters(smoothingParameters);
  _h = smoothingParameters[0]; // For smaller domain counting

  qDebug() << "Reservoir size in step " << _stepNumber
           << " is: " << getClustersForEstimator().size() << ".";

  countKDEValuesOnClusters();
  updatePrognosisParameters();
  updateMaxAbsAVector();
  updateAverageMaxAbsAsInLastKPSSMSteps();
  updateAverageMaxAbsAsInLastMinMSteps();
  updateExaminedClustersAsVector();

  // Start at 0
  stationarityTest->addNewSample(
      std::stod(_clusters->front()->getObject()->attributesValues["Val0"])
  );

  _sgmKPSS = sigmoid(_psi * stationarityTest->getTestsValue() - 11.1); // sgmKPSS
  _r = 0.01 + 0.09 * _sgmKPSS;

  updateWeights();

  emE._currentKDEValue = getNewEmEValue();
  emE.updatePrediction();

  ++_stepNumber;
}

void DESDA::updateWeights()
{
  _examinedClustersWStar.clear();
  auto consideredClusters = getClustersForEstimator();
  auto m = consideredClusters.size();

  for(int i = 0; i < consideredClusters.size(); ++i){
    double newWeight = 2 * (1.0 - i * _sgmKPSS / m);
    consideredClusters[i]->setCWeight(newWeight);
    if(std::count(_examinedClustersIndices.begin(), _examinedClustersIndices.end(), i))
        _examinedClustersWStar.push_back(newWeight);
  }
}

/** DESDA::updateExaminedClustersIndices
 * @brief A function that updates indices of examined clusters. Note that it
 * should be called after m update.
 */
void DESDA::updateExaminedClustersIndices()
{
  auto desiredClustersLocations = {0.2, 0.5, 0.8};
  _examinedClustersIndices.clear();

  auto m = getClustersForEstimator().size();

  for(auto val : desiredClustersLocations)
      _examinedClustersIndices.push_back(int(val * m));
}

std::vector<std::shared_ptr<cluster> > DESDA::getClustersForEstimator()
{
  std::vector<std::shared_ptr<cluster> > consideredClusters = {};

  int i = 0;

  while (consideredClusters.size() < _m && i < _clusters->size()) {
    consideredClusters.push_back((*_clusters)[i]);
    ++i;
  }

  return consideredClusters;
}

std::vector<std::shared_ptr<cluster> > DESDA::getClustersForWindowedEstimator()
{
    std::vector<std::shared_ptr<cluster>> consideredClusters = {};

    for(auto c : *_clusters) consideredClusters.push_back(c);

    return consideredClusters;
}

/** DESDA::enhanceWeightsOfUncommonElements
 * @brief Enhance weights of atypical elements.
 *
 * The method is described in Kulczycki, Kus, Rybotycki 2020
 *
 */
void DESDA::enhanceWeightsOfUncommonElements()
{
  _examinedClustersIndicesInUncommonClustersVector.clear();
  qDebug() << "Enhancing weights of the uncommon.";

  auto uncommonElements = getAtypicalElements();

  std::vector<double> examinedClustersEnhancedWeights = {};

  for(int i = 0; i < uncommonElements.size(); ++i){
    auto ue = uncommonElements[i];
    double weightEnhancer = 2 * sigmoid(ue->predictionParameters[1] / _averageMaxPredictionAInLastKPSSMSteps) - 1;
    weightEnhancer *= _sgmKPSS;
    weightEnhancer += 1;
    if(std::count(_examinedClustersIndicesInUncommonClustersVector.begin(),
                  _examinedClustersIndicesInUncommonClustersVector.end(), i)){
      examinedClustersEnhancedWeights.push_back(weightEnhancer);
    }
    ue->setCWeight(ue->getCWeight() * weightEnhancer);
  }

  int i = 0;
  for(auto val : _examinedClustersIndicesInUncommonClustersVector){
    if(val == -1){
        _examinedClustersWStar3.push_back(1);
    }
    else {
        _examinedClustersWStar3.push_back(i);
        ++i;
    }
  }
}

void DESDA::countKDEValuesOnClusters()
{
  std::vector<double> x;

  auto consideredClusters = getClustersForEstimator();
  _estimator->setClusters(consideredClusters);

  for(std::shared_ptr<cluster> c : *_clusters)
  {
    x.clear();
    x.push_back(std::stod(c->getRepresentative()->attributesValues["Val0"]));
    double estimatorValueOnCluster = _estimator->getValue(&x);
    c->_currentKDEValue = estimatorValueOnCluster;
  }
}

void DESDA::updatePrognosisParameters()
{
  for(std::shared_ptr<cluster> c : *_clusters)
    c->updatePrediction();
}

void DESDA::updateM()
{
  if(emE.predictionParameters.size() < 2) return;
  if(_sgmKPSS /*sgmKPSS*/ < 0) return;

  _m = (_maxM - (_maxM - _minM) * _sgmKPSS);
  _m = _m < _clusters->size() + 1 ? _m : _clusters->size() + 1;
}

void DESDA::updateDelta()
{
    if(emE.predictionParameters.size() < 2) return;

    double sigmoidArg = _s;
    sigmoidArg += _mu * fabs(emE.predictionParameters[1]);
    delta = sigmoid(sigmoidArg);
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

/** DESDA::updateMaxAbsAVector
 * @brief Updates vector of absolute values of a. This vector should store _maxM
 * values at most and it's values should be counted on all stored clusters (as
 * all of them have their predictions updated.
 */
void DESDA::updateMaxAbsAVector()
{
  // Add new value
  _maxAbsAs.insert(_maxAbsAs.begin(), getCurrentMaxAbsA());

  // Ensure size is as expected
  while(_maxAbsAs.size() > _maxM)
      _maxAbsAs.pop_back();
}

/** DESDA::getCurrentMaxAbsA
* @brief Finds and returns current maximal value of abs(a) of all clusters.
* @return Current maximal values of abs(a) of all clusters.
*/
double DESDA::getCurrentMaxAbsA()
{
  if(_clusters->size() < 0) return -1;
  double maxA = fabs((*_clusters)[0]->predictionParameters[1]);
  for(auto c : *_clusters){
    double currentA = fabs(c->predictionParameters[1]);
    maxA = currentA > maxA ? currentA : maxA;
  }

  return maxA;
}

/** DESDA::updateAverageMaxAbsAsInLastKPSSMSteps
 * @brief Updates the value of average max abs(a) in last KPSS M number of steps.
 *
 * Thus should be called after max(abs(a)) vector update.
 */
void DESDA::updateAverageMaxAbsAsInLastKPSSMSteps()
{
    int consideredElementsNumber = _maxAbsAs.size() < _kpssM ? _maxAbsAs.size() : _kpssM;

    double sumOfConsideredAbsAs = 0;

    for(int i = 0; i < consideredElementsNumber; ++i){
      sumOfConsideredAbsAs += _maxAbsAs[i];
    }

    _averageMaxPredictionAInLastKPSSMSteps = sumOfConsideredAbsAs / consideredElementsNumber;
}

void DESDA::updateAverageMaxAbsAsInLastMinMSteps()
{
    double sumOfConsideredAbsAs = 0;
    int consideredElementsNumber = _maxAbsAs.size() < _minM ? _maxAbsAs.size() : _minM;

    for(int i = 0; i < consideredElementsNumber; ++i){
      sumOfConsideredAbsAs += _maxAbsAs[i];
    }

    _averageMaxPredictionAInLastMinMSteps = sumOfConsideredAbsAs / consideredElementsNumber;
}

void DESDA::updateExaminedClustersAsVector()
{
    _examinedClustersAs.clear();

    for(auto val : _clustersExaminedForAsIndices){
      if(_clusters->size() > val)
        _examinedClustersAs.push_back((*_clusters)[val]->predictionParameters[1]);
      else
        _examinedClustersAs.push_back(0);
    }
}

double DESDA::getDomainMinValue(const std::vector<clusterPtr> &clusters)
{
    if(clusters.size() == 0) return - 3 * _h;

    double domainMin = std::stod(clusters[0]->getRepresentative()->attributesValues["Val0"]);

    for(auto c : clusters){
        auto cVal = std::stod(c->getRepresentative()->attributesValues["Val0"]);
        domainMin = cVal < domainMin ? cVal : domainMin;
    }

    return domainMin - 3 * _h;
}

double DESDA::getDomainMaxValue(const std::vector<clusterPtr> &clusters)
{
    if(clusters.size() == 0) return 3 * _h;

    double domainMax = std::stod(clusters[0]->getRepresentative()->attributesValues["Val0"]);

    for(auto c : clusters){
        auto cVal = std::stod(c->getRepresentative()->attributesValues["Val0"]);
        domainMax = cVal > domainMax ? cVal : domainMax;
    }

    return domainMax + 3 * _h;
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

    double domainMinValue = getDomainMinValue(currentClusters);
    double domainMaxValue = getDomainMaxValue(currentClusters);

    for(qreal x: *X)
    {
      if(x > domainMinValue && x < domainMaxValue){
          std::vector<double> pt = {x};
          kernelPrognosisDerivativeValues.push_back(
            _estimatorDerivative->getValue(&pt) * 1000 // For visibility
          );
      } else {
          kernelPrognosisDerivativeValues.push_back(0);
      }
    }
  }

  return kernelPrognosisDerivativeValues;
}

QVector<double> DESDA::getEnhancedKDEValues(const QVector<qreal> *X)
{
  auto currentClusters = getClustersForEstimator();
  auto standardWeights = getClustersWeights(currentClusters);
  sigmoidallyEnhanceClustersWeights(&currentClusters);

  QVector<double> enhancedKDEValues = {};

  _enhancedKDE->setClusters(currentClusters);

  double domainMinValue = getDomainMinValue(currentClusters);
  double domainMaxValue = getDomainMaxValue(currentClusters);

  for(qreal x: *X)
  {
    if(x > domainMinValue && x < domainMaxValue){
        std::vector<double> pt = {x};
        enhancedKDEValues.push_back(_enhancedKDE->getValue(&pt));
    } else {
        enhancedKDEValues.push_back(0);
    }
  }

  // Restore weights
  for(unsigned int i = 0; i < currentClusters.size(); ++i)
    currentClusters[i]->setCWeight(standardWeights[i]);

  return enhancedKDEValues;
}

/** DESDA::getClustersWeights
 * @brief Gets CWeights from given clusters.
 * @param clusters -- clusters to get weights from
 * @return std::vector<double> of weights
 */
std::vector<double> DESDA::getClustersWeights(const std::vector<std::shared_ptr<cluster> > &clusters)
{
  std::vector<double> weights = {};

  for(auto c : clusters)
    weights.push_back(c->getCWeight());

  return weights;
}

/** DESDA::sigmoidallyEnhanceClustersWeights
 * @brief Enhancing weights of considered cluters based on prognosis.
 * @param clusters - clusters to enhance wieghts
 */
void DESDA::sigmoidallyEnhanceClustersWeights(std::vector<std::shared_ptr<cluster> > *clusters)
{
  _examinedClustersWStar2.clear();

  QVector<qreal> clustersXs = {};

  for(auto c : *clusters)
    clustersXs.push_back(std::stod(c->getRepresentative()->attributesValues["Val0"]));

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

  for(int i = 0; i < clusters->size(); ++i)
  {
    auto c = (*clusters)[i];
    enhancedWeight = c->getCWeight();

    // Count v_i
    v_i = 2 * delta * ( 1.0 / (1.0 + exp(- gamma * derivativeVal[i])) - 0.5);

    enhancedWeight *= (1 + v_i);

    if(std::count(_examinedClustersIndices.begin(), _examinedClustersIndices.end(), i))
        _examinedClustersWStar2.push_back(1 + v_i);

    c->setCWeight(enhancedWeight);
  }
}

QVector<double> DESDA::getWindowKDEValues(const QVector<qreal> *X)
{
    QVector<qreal> x;
    QVector<double> windowKDEValues = {};

    auto consideredClusters = getClustersForWindowedEstimator();
    _estimator->setClusters(consideredClusters);
    _estimator->_shouldConsiderWeights = false;

    double domainMinValue = getDomainMinValue(consideredClusters);
    double domainMaxValue = getDomainMaxValue(consideredClusters);

    for(qreal x: *X)
    {
      if(x > domainMinValue && x < domainMaxValue){
        std::vector<double> q = {x};
        windowKDEValues.push_back(_estimator->getValue(&q));
      } else {
        windowKDEValues.push_back(0);
      }
    }

    _estimator->_shouldConsiderWeights = true;    
    _estimator->setClusters(getClustersForEstimator());

    return windowKDEValues;
}

QVector<double> DESDA::getKDEValues(const QVector<qreal> *X)
{
    QVector<qreal> x;
    QVector<double> KDEValues = {};

    auto consideredClusters = getClustersForEstimator();
    _estimator->setClusters(consideredClusters);
    _estimator->_shouldConsiderWeights = false;

    double domainMinValue = getDomainMinValue(consideredClusters);
    double domainMaxValue = getDomainMaxValue(consideredClusters);

    for(qreal x: *X)
    {
      if(x > domainMinValue && x < domainMaxValue){
        std::vector<double> q = {x};
        KDEValues.push_back(_estimator->getValue(&q));
      } else {
        KDEValues.push_back(0);
      }
    }

    return KDEValues;
}

QVector<double> DESDA::getWeightedKDEValues(const QVector<qreal> *X)
{
    QVector<qreal> x;
    QVector<double> weightedKDEValues = {};

    auto consideredClusters = getClustersForEstimator();
    _estimator->setClusters(consideredClusters);
    _estimator->_shouldConsiderWeights = true;

    for(qreal x: *X)
    {
      std::vector<double> q = {x};
      weightedKDEValues.push_back(_estimator->getValue(&q));
    }

    return weightedKDEValues;
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

  double avgME = 0;

  // Counting average
  for(int i = 0; i < m0; ++i){
    avgME += std::stod(_clusters->at(i)->getObject()->attributesValues["Val0"]);
  }

  avgME /= m0;

  // Counting var
  double val, var = 0;

  for(int i = 0; i < m0; ++i){
    val = std::stod(_clusters->at(i)->getObject()->attributesValues["Val0"]);
    var += pow(val - avgME, 2);
  }

  var /= m0 - 1;

  _stDev = pow(var, 0.5);
  return _stDev;
}

cluster DESDA::getEmECluster()
{
  return emE;
}

double DESDA::getStationarityTestValue()
{
  return stationarityTest->getTestsValue();
}

/** DESDA::getAtypicalElements
 * @brief Finds and returns vector of atypical elements.
 *
 * Atypical elements are found using Kulczycki-Kruszewski method.
 *
 * @return Vector of atypical/rare/uncommon elements in _clusters.
 */
std::vector<clusterPtr> DESDA::getAtypicalElements()
{
  auto AKDEValues = getVectorOfAcceleratedKDEValuesOnClusters();
  auto sortedIndicesValues = getSortedAcceleratedKDEValues(AKDEValues);
  recountQuantileEstimatorValue(sortedIndicesValues);
  std::vector<clusterPtr> atypicalElements = {};

  for(int i = 0; i < sortedIndicesValues.size(); ++i){
    if(_quantileEstimator > sortedIndicesValues[i].second){
      atypicalElements.push_back((*_clusters)[sortedIndicesValues[i].first]);

      if(std::count(_examinedClustersIndices.begin(), _examinedClustersIndices.end(), i))
         _examinedClustersIndicesInUncommonClustersVector.push_back(atypicalElements.size() - 1);
    }
    else
      if(std::count(_examinedClustersIndices.begin(), _examinedClustersIndices.end(), i))
        _examinedClustersIndicesInUncommonClustersVector.push_back(-1);
  }

  return atypicalElements;
}

std::vector<double> DESDA::getVectorOfAcceleratedKDEValuesOnClusters()
{
  //qDebug() << "Accelerated KDE Values on clusters.";
  std::vector<double> x;
  auto consideredClusters = getClustersForEstimator();
  auto standardWeights = getClustersWeights(consideredClusters);

  //qDebug() << "Standard weights:" << standardWeights;

  sigmoidallyEnhanceClustersWeights(&consideredClusters);

  std::vector<double> AKDEValues = {};
  auto m = consideredClusters.size();

  for(auto i = 0; i < m; ++i){
    auto c = consideredClusters[0];
    consideredClusters.erase(consideredClusters.begin(), consideredClusters.begin() + 1);
    _enhancedKDE->setClusters(consideredClusters);

    x.push_back(std::stod(c->getRepresentative()->attributesValues["Val0"]));
    AKDEValues.push_back(_enhancedKDE->getValue(&x));

    x.clear();
    consideredClusters.push_back(c);
  }

  // Restore weights
  for(unsigned int i = 0; i < consideredClusters.size(); ++i)
    consideredClusters[i]->setCWeight(standardWeights[i]);

  return AKDEValues;
}

std::vector<std::pair<int, double> > DESDA::getSortedAcceleratedKDEValues(const std::vector<double> &AKDEValues)
{
  // Create pairs containing indexes and values of AKDE
  std::vector<std::pair<int, double>> indexesValues = {};

  for(int i = 0; i < AKDEValues.size(); ++i){
    indexesValues.push_back(std::pair<int, double>(i, AKDEValues[i]));
  }

  // Sort them
  std::sort(
    indexesValues.begin(), indexesValues.end(),
        [](std::pair<int, double> a, std::pair<int, double> b){
      return a.second > b.second;
    }
  );

  std::reverse(indexesValues.begin(), indexesValues.end());

  return indexesValues;
}

void DESDA::recountQuantileEstimatorValue(const std::vector<std::pair<int, double> > &sortedIndicesValues)
{
  int m = sortedIndicesValues.size();
  double mr = _r * m;

  if(mr < 0.5){
      _quantileEstimator = sortedIndicesValues[0].second;
      return;
  }

  int i = mr + 0.5;

  _quantileEstimator = (0.5 + i - mr) * sortedIndicesValues[i - 1].second; // Remembert that indices in the formulas start from 1.
  _quantileEstimator += (0.5 - i + mr) * sortedIndicesValues[i].second; // Remembert that indices in the formulas start from 1.
}

QVector<double> DESDA::getRareElementsEnhancedKDEValues(const QVector<qreal> *X)
{
  QVector<double> enhancedKDEValues = {};
  auto currentClusters = getClustersForEstimator();
  auto standardWeights = getClustersWeights(currentClusters);

  enhanceWeightsOfUncommonElements();
  sigmoidallyEnhanceClustersWeights(&currentClusters);

  _enhancedKDE->setClusters(currentClusters);

  _examinedClustersW.clear();
  for(auto val : _examinedClustersIndices){
      _examinedClustersW.push_back(currentClusters[val]->getCWeight());
  }

  double domainMinValue = getDomainMinValue(currentClusters);
  double domainMaxValue = getDomainMaxValue(currentClusters);

  for(qreal x: *X)
  {
    if(x > domainMinValue && x < domainMaxValue){
        std::vector<double> pt = {x};
        enhancedKDEValues.push_back(_enhancedKDE->getValue(&pt));
    } else {
        enhancedKDEValues.push_back(0);
    }
  }

  // Restore weights
  for(unsigned int i = 0; i < currentClusters.size(); ++i)
    currentClusters[i]->setCWeight(standardWeights[i]);

  //qDebug() << "Max a param: " << maxAParam;
  return enhancedKDEValues;
}

/** DESDA::getAtypicalElementsValuesAndDerivatives
 * @brief Returns vector of current atypical elements values and their derivatives.
 *
 * This method assumes 1 dimensional, specifically constructed objects (for now).
 *
 * @return Vector of pairs of atypical elements values and their derivatives.
 */
QVector<std::pair<double, double>> DESDA::getAtypicalElementsValuesAndDerivatives()
{
  QVector<std::pair<double, double>> atypicalElementsValuesAndDerivatives = {};
  auto atypicalElements = getAtypicalElements();

  for(auto a : atypicalElements){
    std::pair<double, double> valueDerivative = std::pair<double, double>(0, 0);
    valueDerivative.first = std::stod(a->getObject()->attributesValues["Val0"]);
    valueDerivative.second = a->predictionParameters[1];
    atypicalElementsValuesAndDerivatives.push_back(valueDerivative);
  }

  return atypicalElementsValuesAndDerivatives;
}

double DESDA::getMaxAbsAOnLastKPSSMSteps()
{
  if(_maxAbsAs.size() < _kpssM)
      return *(std::max_element(_maxAbsAs.begin(), _maxAbsAs.end()));

  return *(std::max_element(_maxAbsAs.begin() + _maxAbsAs.size() - _kpssM, _maxAbsAs.end()));
}
