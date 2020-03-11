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
  _clustersExaminedForAsIndices = {(unsigned int)(0.2 * _maxM),
                                   (unsigned int)(0.5 * _maxM),
                                   (unsigned int)(0.8 * _maxM)};

  _mE = 1000;
  _m = _maxM;

  _kpssM = _maxM * 0.2;
  int l = kpssX * pow(_kpssM / 100, 0.25);

  _sgmKPSS = -1;

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
  updateWeights();

  _smoothingParamCounter->updateSmoothingParameterValue(
    _weightModifier,
    std::stod(_clusters->front()->getObject()->attributesValues["Val0"])
  );

  _smoothingParamCounter->setClusters(_clusters, 0);

  std::vector<double> smoothingParameters = { _smoothingParamCounter->countSmoothingParameterValue() * _smoothingParameterMultiplier };

  _estimator->setSmoothingParameters(smoothingParameters);
  _estimatorDerivative->setSmoothingParameters(smoothingParameters);
  _enhancedKDE->setSmoothingParameters(smoothingParameters);

  countKDEValuesOnClusters();
  updatePrognosisParameters();
  updateMaxAbsAVector();
  updateAverageMaxAbsAsInLastKPSSMSteps();
  updateExaminedClustersAsVector();

  auto currentClusters = getClustersForEstimator();

  qDebug() << "Reservoir size in step "
             << _stepNumber << " is: " << currentClusters.size();

  //avg = getAverageOfFirstMSampleValues(_mE);
  avg = getNewEmEValue();

  // Start at 0
  stationarityTest->addNewSample(
      std::stod(_clusters->front()->getObject()->attributesValues["Val0"])
  );

  _sgmKPSS = sigmoid(_psi * stationarityTest->getTestsValue() - 11.1); // sgmKPSS

  emE._currentKDEValue = avg;

  _estimator->setClusters(currentClusters);

  emE.updatePrediction();

  ++_stepNumber;
}

void DESDA::updateWeights()
{
  _examinedClustersWStar.clear();

  // p weight update
  if(_stepNumber > _maxM){
      double weightsSum = 0;

      for(int clusterNum = _clusters->size() - 2; clusterNum > -1; --clusterNum){
        weightsSum += (*_clusters)[clusterNum]->getWeight();
      }

      for(int clusterNum = _clusters->size() - 2; clusterNum > -1; --clusterNum){
        (*_clusters)[clusterNum]->setWeight((_m - 1) * (*_clusters)[clusterNum]->getWeight() / weightsSum);
      }
  }

  double cWeightSum = 0;

  for(int clusterNum = _clusters->size() - 1; clusterNum > -1; --clusterNum)
  {
      // In formula it's (i - 1), but indexes are from 1 not 0, thus no -1.
      double newWeight =
          (*_clusters)[clusterNum]->getWeight() * (1.0 - _newWeightB * (clusterNum) / _m);

      newWeight = std::max(0.0, newWeight);

      (*_clusters)[clusterNum]->setCWeight(newWeight);

     cWeightSum += newWeight;
  }

  for(int clusterNum = _clusters->size() - 1; clusterNum > -1; --clusterNum){
      (*_clusters)[clusterNum]->setCWeight((*_clusters)[clusterNum]->getCWeight() * _m / cWeightSum);

      // Weights star
      if(std::count(_examinedClustersIndices.begin(), _examinedClustersIndices.end(), clusterNum)){
          this->_examinedClustersWStar.push_back((*_clusters)[clusterNum]->getCWeight());
      }
  }

  std::reverse(_examinedClustersWStar.begin(), _examinedClustersWStar.end());
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
  auto uncommonElements = getAtypicalElements();

  std::vector<double> _examinedClustersEnhancedWeights = {};

  for(int i = 0; i < uncommonElements.size(); ++i){
    auto ue = uncommonElements[i];
    double weightEnhancer = 2 * sigmoid(ue->predictionParameters[1] / _averageMaxPredictionAInLastKPSSMSteps) - 1;
    weightEnhancer *= _sgmKPSS;
    weightEnhancer += 1;
    if(std::count(_examinedClustersIndicesInUncommonClustersVector.begin(),
                  _examinedClustersIndicesInUncommonClustersVector.end(), i))
      _examinedClustersEnhancedWeights.push_back(weightEnhancer);
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
  int m = 0;

  if(emE.predictionParameters.size() < 2) return;
  if(_sgmKPSS /*sgmKPSS*/ < 0) return;

  // Old m
  //m = round(1.05 * _maxM * ( 1 - _lambda * _u_i * fabs(emE.predictionParameters[1]))); // getStdDevOfFirstMSampleValues(_mE)));
  // 8 X 2019 m
  // m = round(1.1 * _maxM * ( 1 - _lambda * fabs(emE.predictionParameters[1]))); // getStdDevOfFirstMSampleValues(_mE)));
  // 11 III 2020, article
  m = int(_maxM * (1 - 0.8 * _sgmKPSS));

  /* This is not needed in article's approach.
  m = std::max(m, _minM);
  _m = std::min(m, _maxM);
  */
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
  _maxAbsAs.push_back(getCurrentMaxAbsA());

  // Ensure size is as expected
  while(_maxAbsAs.size() > _maxM)
      _maxAbsAs.erase(_maxAbsAs.begin(), _maxAbsAs.begin() + 1);
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
    int offset = _maxAbsAs.size() - _kpssM;
    offset = offset < 0 ? 0 : offset;
    double sumOfConsideredAs =
            std::accumulate(_maxAbsAs.begin() + offset, _maxAbsAs.end(), 0);
    int consideredElementsNumber = _maxAbsAs.size() < _kpssM ? _maxAbsAs.size() : _kpssM;
    _averageMaxPredictionAInLastKPSSMSteps = sumOfConsideredAs / consideredElementsNumber;
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
      std::vector<double> pt;
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
  auto currentClusters = getClustersForEstimator();
  auto standardWeights = getClustersWeights(currentClusters);
  sigmoidallyEnhanceClustersWeights(&currentClusters);

  QVector<double> enhancedKDEValues = {};

  _enhancedKDE->setClusters(currentClusters);

  for(qreal x: *X) {
    std::vector<double> pt;
    pt.push_back(x);
    enhancedKDEValues.push_back(_enhancedKDE->getValue(&pt));
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

    for(qreal x: *X)
    {
      std::vector<double> q = {x};
      windowKDEValues.push_back(_estimator->getValue(&q));
    }

    _estimator->_shouldConsiderWeights = true;    
    _estimator->setClusters(getClustersForEstimator());

    return windowKDEValues;
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
  std::vector<double> x;
  auto consideredClusters = getClustersForEstimator();
  auto standardWeights = getClustersWeights(consideredClusters);
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

double DESDA::recountQuantileEstimatorValue(const std::vector<std::pair<int, double> > &sortedIndicesValues)
{
  int m = sortedIndicesValues.size();
  double mr = _r * m;

  if(mr <= 0.5) return sortedIndicesValues[0].second;
  if(mr >= m - 0.5) return sortedIndicesValues[m - 1].second;

  int i = mr + 0.5 - 1; // Indices of elements in Kruszewski starts from 1

  _quantileEstimator = (0.5 + i - mr) * sortedIndicesValues[i].second;
  _quantileEstimator += (0.5 - i + mr) * sortedIndicesValues[i + 1].second;
}

QVector<double> DESDA::getRareElementsEnhancedKDEValues(const QVector<qreal> *X)
{
  //qDebug() << "EKDE values getting.";

  QVector<double> enhancedKDEValues = {};
  auto currentClusters = getClustersForEstimator();
  auto standardWeights = getClustersWeights(currentClusters);

  sigmoidallyEnhanceClustersWeights(&currentClusters);
  enhanceWeightsOfUncommonElements();
  _enhancedKDE->setClusters(currentClusters);

  _examinedClustersW.clear();
  for(auto val : _examinedClustersIndices){
      _examinedClustersW.push_back(currentClusters[val]->getCWeight());
  }

  for(qreal x: *X)
  {
    std::vector<double> pt;
    pt.push_back(x);
    enhancedKDEValues.push_back(
      _enhancedKDE->getValue(&pt)
    );
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
