#include "DESDA.h"
#include "KDE/pluginsmoothingparametercounter.h"

#include <QTime>
#include <QCoreApplication>
#include <QRandomGenerator>
#include <algorithm>
#include <fstream>
#include <math.h>

# define M_PI           3.14159265358979323846  /* pi */

DESDA::DESDA(std::shared_ptr<kernelDensityEstimator> estimator,
             std::shared_ptr<kernelDensityEstimator> estimatorDerivative,
             std::shared_ptr<kernelDensityEstimator> enchancedKDE,
             reservoirSamplingAlgorithm *samplingAlgorithm,
             std::vector<std::shared_ptr<cluster> > *clusters,
             double desiredRarity, double pluginRank) :
    _samplingAlgorithm(samplingAlgorithm),
    _estimatorDerivative(estimatorDerivative), _estimator(estimator),
    _clusters(clusters), _r(desiredRarity),
    _enhancedKDE(enchancedKDE), _pluginRank(pluginRank) {
  _objects.clear();

  _maxM = _samplingAlgorithm->getReservoidMaxSize(); // _maxM should be like 1000 + 100 for every mode

  _m = _maxM;

  _minM = _maxM / 10; // This works for both 2D and 1D experiments with default settings.
  max_prognosis_error_clusters_ = 600;
  _kpssM = 600; // This is independent of maxM. Normally 500.

  _sgmKPSS = -1;
  _sgmKPSSPercent = 30;
  _stepNumber = 1;
  _smoothingParameterEnhancer = 1;

  prognosis_cluster_ = cluster(-1);
  error_domain_points_number_ = 500;

  for(int i = 0; i < estimator->getDimension(); ++i) {
    stationarityTests.push_back(std::make_shared<KPSSStationarityTest>(_kpssM));
    prognosis_clusters_.push_back(cluster(-1));
    prognosis_errors_.push_back(std::vector<double>());
  }

}

int sgn(double val) {
  return (val > 0) - (val < 0);
}

double sum(std::vector<double> vals) {
  double valsSum = 0.0;

  for(auto val : vals) valsSum += val;

  return valsSum;
}

double sigmoid(double x) {
  return 1.0 / (1.0 + exp(-x));
}

double average(const std::vector<double> &values) {
  double average = 0;

  for(auto val : values)
    average += val;

  average /= values.size();

  return average;
}

double var(const std::vector<double> &v){
  auto n = v.size();

  if(n < 2) { return 0; }

  double var = 0;

  for(auto val : v){ var += val; }

  var *= -var / n;

  for(auto val : v){ var += val * val; }

  return var / (n - 1);
}

double stdev(const std::vector<double> &v){
  return std::pow(var(v), 0.5);
}

void DESDA::performStep() {

  // Making place for new cluster
  while(_clustersForWindowed.size() >= _maxM) {
    _clustersForWindowed.pop_back();
    _objects.erase(_objects.begin(), _objects.begin() + 1);
  }

  while(_clusters->size() > _m) {
    _clusters->pop_back();
  }

  // Reservoir movement
  _samplingAlgorithm->performSingleStep(&_objects, _stepNumber);

  std::shared_ptr<cluster> newCluster =
      std::shared_ptr<cluster>(new cluster(_stepNumber, _objects.back()));
  newCluster->setTimestamp(_stepNumber);

  for(int i = 0; i < stationarityTests.size(); ++i) {
    std::string attribute = (*newCluster->getObject()->attirbutesOrder)[i];
    stationarityTests[i]->addNewSample(std::stod(newCluster->getObject()->attributesValues[attribute]));
  }

  _sgmKPSS = sigmoid(_sgmKPSSParameters[_sgmKPSSPercent][0] * getStationarityTestValue()
                     - _sgmKPSSParameters[_sgmKPSSPercent][1]);
  // Beta0 update
  _beta0 = 2.0 / 3 * _sgmKPSS; // According to formula from 13 IV 2020

  _clusters->insert(_clusters->begin(), newCluster);
  _clustersForWindowed.insert(_clustersForWindowed.begin(), newCluster);

  // M update
  updateM();
  updateExaminedClustersIndices(); // For labels update

  // DEBUG
  //ms.push_back(_m);
  //qDebug() << "ms: " << ms;
  // DEBUG

  _v =  1.0 - 1.0 / _minM;
  cluster::_deactualizationParameter = _v;

  // Calculate smoothing parameters
  compute_weighted_plugin = false;
  _windowedSmoothingParametersVector = calculateH(*_clusters);
  auto currentClusters = getClustersForEstimator();

  // Update weights
  updateWeights();
  compute_weighted_plugin = true;
  _smoothingParametersVector = calculateH(currentClusters);

  // DEBUG
  /*
  for(auto c : currentClusters){
    qDebug() << "Cluster: " << std::stod(c->getRepresentative()->attributesValues["Val0"]) << ", w = " << c->getCWeight();
  }
  //*/
  // DEBUG

  qDebug() << "Reservoir size in step " << _stepNumber
           << " is: " << currentClusters.size() << ".";

  // Update clusters prognosis
  countKDEValuesOnClusters();

  for(size_t i = 0; i < stationarityTests.size(); ++i){
    std::string attr_key = "Val" + std::to_string(i);
    prognosis_clusters_[i]._currentKDEValue = std::stod(_objects.back()->attributesValues[attr_key]);
    // DEBUG //
      //vals.push_back(std::stod(_objects.back()->attributesValues[attr_key]));
    // DEBUG //
  }

  updatePrognosisParameters();
  countDerivativeValuesOnClusters();
  updateMaxAbsDerivativeInCurrentStep();

  statistics_.clear();

  for(size_t i = 0; i < prognosis_errors_.size(); ++i){

    while(prognosis_errors_[i].size() >= max_prognosis_error_clusters_){
      prognosis_errors_[i].pop_back();
    }

    if(_stepNumber != 1) {
      prognosis_errors_[i].insert(prognosis_errors_[i].begin(),
                                  prognosis_clusters_[i]._currentKDEValue - prognosis_clusters_[i].getLastPrediction());
    }

    // DEBUG //
      // qDebug() << prognosis_clusters_[i]._currentKDEValue << " - " << prognosis_clusters_[i].getLastPrediction() << " = " << prognosis_clusters_[i]._currentKDEValue - prognosis_clusters_[i].getLastPrediction();
      // errors.push_back(prognosis_clusters_[i]._currentKDEValue - prognosis_clusters_[i].getLastPrediction());
      // progs.push_back(prognosis_clusters_[i].getLastPrediction());
    // DEBUG //

    e_ = prognosis_errors_[i].empty() ? 0 : prognosis_errors_[i][0];
    statistics_.push_back(ComputeStatistics(prognosis_errors_[i]));
    // DEBUG //
      // stats.push_back(statistics_[0]);
    // DEBUG //

    avg = average(prognosis_errors_[i]);
    std = stdev(prognosis_errors_[i]);

    //prognosis_cluster_.updatePrediction();
    prognosis_clusters_[i].updatePrediction();
  }

  _examinedClustersDerivatives.clear();
  for(auto index : _examinedClustersIndices) {
    if(index < 0) {
      _examinedClustersDerivatives.push_back(0);
    } else {
      _examinedClustersDerivatives.push_back(currentClusters[index]->_currentDerivativeValue);
    }
  }

  // Update uncommon elements
  _r = 0.01 + 0.09 * _sgmKPSS;

  // DEBUG //
    //qDebug() << "Progs: " << progs;
    //qDebug() << "Vals: " << vals;
    //qDebug() << "Stats: " << stats;
    //qDebug() << "\n\nErr: " << prognosis_errors_[0];
  // DEBUG //

  ++_stepNumber;
}

void DESDA::updateWeights() {
  _examinedClustersWStar.clear();
  auto consideredClusters = getClustersForEstimator();
  auto m = consideredClusters.size();

  for(int j = 0; j < std::count(_examinedClustersIndices.begin(), _examinedClustersIndices.end(), -1); ++j) {
    _examinedClustersWStar.push_back(0);
  }

  for(int i = 0; i < consideredClusters.size(); ++i) {
    double newWeight = 2 * (1.0 - i * _sgmKPSS / m);
    consideredClusters[i]->setCWeight(newWeight);
    for(int j = 0; j < std::count(_examinedClustersIndices.begin(), _examinedClustersIndices.end(), i); ++j)
      _examinedClustersWStar.push_back(newWeight);
  }
}

/** DESDA::updateExaminedClustersIndices
 * @brief A function that updates indices of examined clusters. Note that it
 * should be called after m update.
 */
void DESDA::updateExaminedClustersIndices() {
  auto desiredClustersLocations = {0.2, 0.5, 0.8};
  _examinedClustersIndices.clear();

  auto m = getClustersForEstimator().size();

  for(auto val : desiredClustersLocations)
    _examinedClustersIndices.push_back(round(val * m) - 1);
}

std::vector<std::shared_ptr<cluster> > DESDA::getClustersForEstimator() {
  std::vector<std::shared_ptr<cluster> > consideredClusters = {};
  int i = 0;

  while(consideredClusters.size() < _m && i < _clusters->size()) {
    consideredClusters.push_back((*_clusters)[i]);
    ++i;
  }

  return consideredClusters;
}

std::vector<std::shared_ptr<cluster> > DESDA::getClustersForWindowedEstimator() {
  return _clustersForWindowed;
}

/** DESDA::enhanceWeightsOfUncommonElements
 * @brief Enhance weights of atypical elements.
 *
 * The method is described in Kulczycki, Kus, Rybotycki 2020
 *
 */
void DESDA::enhanceWeightsOfUncommonElements() {
  auto uncommonElements = getAtypicalElements();
  std::vector<double> examinedClustersEnhancedWeights = {};
  std::vector<double> weightsEnhancers = {};

  for(int i = 0; i < uncommonElements.size(); ++i) {
    auto ue = uncommonElements[i];
    double weightEnhancer = ue->_currentDerivativeValue /
                            _maxAbsDerivativeValueInCurrentStep;
    weightEnhancer *= _sgmKPSS;
    weightEnhancer += 1;
    weightsEnhancers.push_back(weightEnhancer);
    ue->setCWeight(ue->getCWeight() * weightEnhancer);
  }

  // For labels
  _examinedClustersWStar3.clear();
  auto clusters = getClustersForEstimator();
  std::vector<std::string> examinedClustersIds = {};

  for(auto idx : _examinedClustersIndices) {
    if(idx < 0) examinedClustersIds.push_back("");
    else examinedClustersIds.push_back(clusters[idx]->getClustersId());
  }

  bool wasExaminedClusterUncommon = false;

  for(int i = 0; i < _examinedClustersIndices.size(); ++i) {

    wasExaminedClusterUncommon = false;

    for(int j = 0; j < uncommonElements.size(); ++j) {
      if(examinedClustersIds[i].compare(uncommonElements[j]->getClustersId()))
        continue;
      wasExaminedClusterUncommon = true;
      _examinedClustersWStar3.push_back(weightsEnhancers[j]);
      break;
    }
    if(!wasExaminedClusterUncommon)
      _examinedClustersWStar3.push_back(1);
  }
}

void DESDA::countKDEValuesOnClusters() {
  std::vector<double> x;
  auto consideredClusters = getClustersForEstimator();
  _estimator->setClusters(consideredClusters);

  _estimator->setSmoothingParameters(_smoothingParametersVector);

  for(std::shared_ptr<cluster> c : *_clusters) {
    x.clear();
    for(auto attribute: *c->getRepresentative()->attirbutesOrder) {
      x.push_back(std::stod(c->getRepresentative()->attributesValues[attribute]));
    }
    double estimatorValueOnCluster = _estimator->getValue(&x);
    c->_currentKDEValue = estimatorValueOnCluster;
  }
}

void DESDA::updatePrognosisParameters() {
  for(std::shared_ptr<cluster> c : *_clusters)
    c->updatePrediction();
}

/** @brief DESDA::countDerivativeValuesOnClusters
 *  Calculates derivative values on points represented by current clusters and
 *  assings them to clusters.
 */
void DESDA::countDerivativeValuesOnClusters() {
  // Get the domain. Formally only m would be needed, but it will not hurt
  // to count on whole domain.
  QVector<std::vector<double>> domain = {};

  for(auto c : *_clusters){
    domain.push_back(std::vector<double>());

    domain.back().push_back(std::stod(c->getObject()->attributesValues["Val0"]));

    if(stationarityTests.size() > 1){
      domain.back().push_back(std::stod(c->getObject()->attributesValues["Val1"]));
    }
  }

  auto derivativeValues = getKernelPrognosisDerivativeValues(&domain);

  for(int i = 0; i < _clusters->size(); ++i)
    (*_clusters)[i]->_currentDerivativeValue = derivativeValues[i];
}

void DESDA::updateM() {
  if(_sgmKPSS < 0) return;

  //_m = round(1.1 * _maxM - (1.1 * _maxM - 0.9 * _minM) * _sgmKPSS);
  _m = round(1.1 * _maxM * (1.0 - _sgmKPSS));
  _m = _m < _minM ? _minM : _m;
  _m = _clusters->size() < _m ? _clusters->size() : _m;
  _m = _m > _maxM ? _maxM : _m;
}

void DESDA::updateMaxAbsDerivativeInCurrentStep() {
  auto consideredClusters = getClustersForEstimator();
  _maxAbsDerivativeValueInCurrentStep = 0;
  for(auto c : consideredClusters) {
    auto cAbsDerivativeValue = fabs(c->_currentDerivativeValue);
    _maxAbsDerivativeValueInCurrentStep =
        _maxAbsDerivativeValueInCurrentStep > cAbsDerivativeValue ?
        _maxAbsDerivativeValueInCurrentStep
                                                                  : cAbsDerivativeValue;
  }
}

std::vector<double> DESDA::getAttributesValuesFromClusters(std::vector<std::shared_ptr<cluster> > clusters,
                                                           int dimension) {
  std::string attribute = (*_samplingAlgorithm->getAttributesList())[dimension];
  std::vector<double> values = {};

  for(auto c: clusters) {
    values.push_back(
        std::stod(c->getRepresentative()->attributesValues[attribute])
                    );
  }

  return values;
}

double DESDA::getDomainMinValue(const std::vector<double> &values, double h) {
  if(values.size() == 0) return -5 * h;

  double domainMin = values[0];

  for(auto val : values) {
    domainMin = val < domainMin ? val : domainMin;
  }

  return domainMin - 5 * h;
}

double DESDA::getDomainMaxValue(const std::vector<double> &values, double h) {
  if(values.size() == 0) return 5 * h;

  double domainMax = values[0];

  for(auto value : values) {
    domainMax = value > domainMax ? value : domainMax;
  }

  return domainMax + 5 * h;
}

QVector<double> DESDA::getErrorDomain(int dimension) {
  std::vector<std::shared_ptr<cluster>> currentClusters =
      getClustersForEstimator();
  std::vector<double> attributesValues =
      getAttributesValuesFromClusters(currentClusters, dimension);
  double domainMinValue = getDomainMinValue(attributesValues, _smoothingParametersVector[dimension]);
  double domainMaxValue = getDomainMaxValue(attributesValues, _smoothingParametersVector[dimension]);
  QVector<double> domain = {};
  double stepSize = (domainMaxValue - domainMinValue) / (error_domain_points_number_);

  for(auto val = domainMinValue; val < domainMaxValue; val += stepSize) {
    domain.push_back(val);
  }

  return domain;
}

QVector<double> DESDA::getWindowedErrorDomain(int dimension) {
  std::vector<std::shared_ptr<cluster>> currentClusters =
      getClustersForWindowedEstimator();
  std::vector<double> attributesValues =
      getAttributesValuesFromClusters(currentClusters, dimension);
  double domainMinValue =
      getDomainMinValue(attributesValues, _windowedSmoothingParametersVector[dimension]);
  double domainMaxValue =
      getDomainMaxValue(attributesValues, _windowedSmoothingParametersVector[dimension]);
  QVector<double> domain = {};
  double stepSize = (domainMaxValue - domainMinValue) / (error_domain_points_number_);

  for(auto val = domainMinValue; val <= domainMaxValue; val += stepSize) {
    domain.push_back(val);
  }

  return domain;
}

std::vector<double> DESDA::calculateH(const std::vector<clusterPtr> &clusters) {
  int dimensionsNumber = _estimator->getDimension();
  std::vector<double> smoothingParameters = {};

  if(clusters.size() < 2) {
    for(int i = 0; i < dimensionsNumber; ++i) {
      smoothingParameters.push_back(1);
    }
    return smoothingParameters;
  }

  // Plugin method.
  QVector<qreal> samples = {};
  QVector<qreal> weights = {};

  for(auto attribute: *_samplingAlgorithm->getAttributesList()) {
    samples.clear();
    for(auto c: clusters) {

      samples.append(std::stod(c->getRepresentative()->attributesValues[attribute]));

      if(compute_weighted_plugin){
        weights.append(c->getCWeight());
      } else {
        weights.append(1);
      }
    }

    pluginSmoothingParameterCounter counter(&samples, _pluginRank, &weights);
    smoothingParameters.push_back(counter.countSmoothingParameterValue()
                                  * _smoothingParameterEnhancer);
  }

  return smoothingParameters;
}

QVector<double> DESDA::getKernelPrognosisDerivativeValues(const QVector<std::vector<double>> *X, int dimension) {
  std::vector<std::shared_ptr<cluster>> currentClusters
      = getClustersForEstimator();
  std::vector<double> prognosisCoefficients = {};
  QVector<double> kernelPrognosisDerivativeValues = {};

  for(auto c : currentClusters)
    prognosisCoefficients.push_back(c->predictionParameters[1]);

  if(prognosisCoefficients.size() == currentClusters.size()) {
    _estimatorDerivative->setAdditionalMultipliers(prognosisCoefficients);
    _estimatorDerivative->setSmoothingParameters({_smoothingParametersVector});
    _estimatorDerivative->setClusters(currentClusters);

    if(stationarityTests.size() == 1) {

      std::string attributeKey =
          (*_samplingAlgorithm->getAttributesList())[dimension];
      std::vector<double> attributesValues =
          getAttributesValuesFromClusters(currentClusters, dimension);
      double domainMinValue =
          getDomainMinValue(attributesValues, _windowedSmoothingParametersVector[0]);
      double domainMaxValue =
          getDomainMaxValue(attributesValues, _windowedSmoothingParametersVector[0]);

      for(auto x: *X) {
        if(x[0] > domainMinValue && x[0] < domainMaxValue) {
          kernelPrognosisDerivativeValues.push_back(
              _estimatorDerivative->getValue(&x) * 1000);// 1000 for visibility
        }
        else {
          kernelPrognosisDerivativeValues.push_back(0);
        }
      }
    }
    else {
      for(auto x: *X) {
        kernelPrognosisDerivativeValues.push_back(_estimatorDerivative->getValue(&x) * 1000); // 1000 for visibility
      }
    }
  }

  return kernelPrognosisDerivativeValues;
}

std::vector<double> DESDA::getEnhancedKDEValues(const std::vector<std::vector<double>> *X, int dimension) {
  auto currentClusters = getClustersForEstimator();
  auto standardWeights = getClustersWeights(currentClusters);
  sigmoidallyEnhanceClustersWeights(&currentClusters);

  std::vector<double> enhancedKDEValues = {};

  _enhancedKDE->setClusters(currentClusters);
  _enhancedKDE->setSmoothingParameters({_smoothingParametersVector});

  std::vector<double> attributesValues =
      getAttributesValuesFromClusters(currentClusters, dimension);
  double domainMinValue =
      getDomainMinValue(attributesValues, _windowedSmoothingParametersVector[0]);
  double domainMaxValue =
      getDomainMaxValue(attributesValues, _windowedSmoothingParametersVector[0]);

  for(auto x: *X) {
    if(x[dimension] > domainMinValue && x[dimension] < domainMaxValue) {
      enhancedKDEValues.push_back(_enhancedKDE->getValue(&x));
    }
    else {
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
std::vector<double> DESDA::getClustersWeights(const std::vector<std::shared_ptr<cluster> > &clusters) {
  std::vector<double> weights = {};

  for(auto c : clusters)
    weights.push_back(c->getCWeight());

  return weights;
}

/** DESDA::sigmoidallyEnhanceClustersWeights
 * @brief Enhancing weights of considered cluters based on prognosis.
 * @param clusters - clusters to enhance wieghts
 */
void DESDA::sigmoidallyEnhanceClustersWeights(std::vector<std::shared_ptr<cluster> > *clusters) {
  _examinedClustersWStar2.clear();

  for(auto index : _examinedClustersIndices)
    if(index < 0) _examinedClustersWStar2.push_back(0);

  for(int i = 0; i < clusters->size(); ++i) {
    auto c = (*clusters)[i];
    double beta = _beta0 * c->_currentDerivativeValue /
                  _maxAbsDerivativeValueInCurrentStep;
    if(_maxAbsDerivativeValueInCurrentStep == 0)
      beta = 0;

    double weightEnhancement = 1 + _sgmKPSS * beta;
    c->setCWeight(c->getCWeight() * weightEnhancement);

    for(int j = 0; j < std::count(_examinedClustersIndices.begin(), _examinedClustersIndices.end(), i); ++j)
      _examinedClustersWStar2.push_back(weightEnhancement);
  }
}

vector<double> DESDA::getWindowKDEValues(const vector<vector<qreal>> *X, int dimension) {
  vector<double> windowKDEValues = {};
  auto currentClusters = getClustersForWindowedEstimator();

  _estimator->setClusters(currentClusters);
  _estimator->setSmoothingParameters({_windowedSmoothingParametersVector});
  _estimator->_shouldConsiderWeights = false;

  std::vector<double> attributesValues =
      getAttributesValuesFromClusters(currentClusters, dimension);
  double domainMinValue =
      getDomainMinValue(attributesValues, _windowedSmoothingParametersVector[0]);
  double domainMaxValue =
      getDomainMaxValue(attributesValues, _windowedSmoothingParametersVector[0]);

  for(auto x: *X) {
    if(x[dimension] > domainMinValue && x[dimension] < domainMaxValue) {
      windowKDEValues.push_back(_estimator->getValue(&x));
    }
    else {
      windowKDEValues.push_back(0);
    }
  }

  _estimator->_shouldConsiderWeights = true;
  _estimator->setClusters(getClustersForEstimator());

  return windowKDEValues;
}

std::vector<double> DESDA::getKDEValues(const vector<vector<double>> *X, int dimension) {
  std::vector<double> KDEValues = {};
  auto currentClusters = getClustersForEstimator();
  _estimator->setClusters(currentClusters);
  _estimator->setSmoothingParameters({_smoothingParametersVector});
  _estimator->_shouldConsiderWeights = false;

  std::vector<double> attributesValues =
      getAttributesValuesFromClusters(currentClusters, dimension);
  double domainMinValue =
      getDomainMinValue(attributesValues, _windowedSmoothingParametersVector[0]);
  double domainMaxValue =
      getDomainMaxValue(attributesValues, _windowedSmoothingParametersVector[0]);

  for(auto x: *X) {
    if(x[dimension] > domainMinValue && x[dimension] < domainMaxValue) {
      KDEValues.push_back(_estimator->getValue(&x));
    }
    else {
      KDEValues.push_back(0);
    }
  }

  return KDEValues;
}

std::vector<double> DESDA::getWeightedKDEValues(const vector<vector<double>> *X, int dimension) {
  std::vector<double> weightedKDEValues = {};
  auto currentClusters = getClustersForEstimator();
  _estimator->setClusters(currentClusters);
  _estimator->setSmoothingParameters({_smoothingParametersVector});
  _estimator->_shouldConsiderWeights = true;

  std::vector<double> attributesValues =
      getAttributesValuesFromClusters(currentClusters, dimension);
  double domainMinValue =
      getDomainMinValue(attributesValues, _windowedSmoothingParametersVector[0]);
  double domainMaxValue =
      getDomainMaxValue(attributesValues, _windowedSmoothingParametersVector[0]);

  for(auto x: *X) {
    if(x[dimension] > domainMinValue && x[dimension] < domainMaxValue) {
      weightedKDEValues.push_back(_estimator->getValue(&x));
    }
    else {
      weightedKDEValues.push_back(0);
    }
  }

  return weightedKDEValues;
}

double DESDA::getStationarityTestValue() {
  std::vector<double> stationarityTestsValues = {};
  for(auto test: stationarityTests) {
    stationarityTestsValues.push_back(test->getTestsValue());
  }
  // Max
  return *std::max_element(stationarityTestsValues.begin(), stationarityTestsValues.end());
}

void DESDA::prepareEstimatorForContourPlotDrawing() {
  auto currentClusters = getClustersForEstimator();
  _unmodifiedCWeightsOfClusters = getClustersWeights(*_clusters);

  sigmoidallyEnhanceClustersWeights(&currentClusters);
  enhanceWeightsOfUncommonElements();

  _estimator->setClusters(currentClusters);
  _estimator->setSmoothingParameters({_smoothingParametersVector});
}

void DESDA::restoreClustersCWeights() {
  for(int i = 0; i < _clusters->size(); ++i) {
    (*_clusters)[i]->setCWeight(_unmodifiedCWeightsOfClusters[i]);
  }
}

/** DESDA::getAtypicalElements
 * @brief Finds and returns vector of atypical elements.
 *
 * Atypical elements are found using Kulczycki-Kruszewski method.
 *
 * @return Vector of atypical/rare/uncommon elements in _clusters.
 */
std::vector<clusterPtr> DESDA::getAtypicalElements() {
  auto AKDEValues = getVectorOfAcceleratedKDEValuesOnClusters();
  auto sortedIndicesValues = getSortedAcceleratedKDEValues(AKDEValues);
  recountQuantileEstimatorValue(sortedIndicesValues);
  auto considered_clusters = getClustersForEstimator();
  std::vector<clusterPtr> atypicalElements = {};

  _trendsNumber = 0;

  // qDebug() << "Quantile estimator:" << _quantileEstimator;

  for(int i = 0; i < sortedIndicesValues.size(); ++i) {

    //qDebug() << sortedIndicesValues[i].second;

    if(_quantileEstimator > sortedIndicesValues[i].second) {
      atypicalElements.push_back((considered_clusters)[sortedIndicesValues[i].first]);
      if((considered_clusters)[sortedIndicesValues[i].first]->_currentDerivativeValue > 0){
        ++_trendsNumber;
      }
    }
    else{
      break;
    }
  }

  _rareElementsNumber = atypicalElements.size();

  return atypicalElements;
}

std::vector<double> DESDA::getVectorOfAcceleratedKDEValuesOnClusters() {
  std::vector<double> x;
  auto consideredClusters = getClustersForEstimator();
  auto standardWeights = getClustersWeights(consideredClusters);

  if(consideredClusters.size() == 1)
    return {consideredClusters[0]->_currentKDEValue};

  sigmoidallyEnhanceClustersWeights(&consideredClusters);

  std::vector<double> AKDEValues = {};
  auto m = consideredClusters.size();

  _enhancedKDE->setSmoothingParameters({_smoothingParametersVector});
  _enhancedKDE->_shouldConsiderWeights = true;

  for(auto i = 0; i < m; ++i) {
    auto c = consideredClusters[0];
    consideredClusters.erase(consideredClusters.begin(), consideredClusters.begin() + 1);
    _enhancedKDE->setClusters(consideredClusters);

    for(auto attribute : *c->getRepresentative()->attirbutesOrder) {
      x.push_back(std::stod(c->getRepresentative()->attributesValues[attribute]));
    }
    AKDEValues.push_back(_enhancedKDE->getValue(&x));

    x.clear();
    consideredClusters.push_back(c);
  }

  // Restore weights
  for(unsigned int i = 0; i < consideredClusters.size(); ++i)
    consideredClusters[i]->setCWeight(standardWeights[i]);

  return AKDEValues;
}

std::vector<std::pair<int, double> > DESDA::getSortedAcceleratedKDEValues(const std::vector<double> &AKDEValues) {
  // Create pairs containing indexes and values of AKDE
  std::vector<std::pair<int, double>> indexesValues = {};

  for(int i = 0; i < AKDEValues.size(); ++i) {
    indexesValues.push_back(std::pair<int, double>(i, AKDEValues[i]));
  }

  // Sort them
  std::sort(
      indexesValues.begin(), indexesValues.end(),
      [](std::pair<int, double> a, std::pair<int, double> b) {
        return a.second > b.second;
      }
           );

  std::reverse(indexesValues.begin(), indexesValues.end());

  //qDebug() << "VALS";
  //for(auto val : indexesValues){
  //  qDebug() << val;
  //}

  return indexesValues;
}

void DESDA::recountQuantileEstimatorValue(const std::vector<std::pair<int, double> > &sortedIndicesValues) {
  int m = sortedIndicesValues.size();
  double mr = _r * m;
  bool originalShouldConsiderWeights = _estimator->_shouldConsiderWeights;
  _estimator->_shouldConsiderWeights = false;

  if(mr < 0.5) {
    _quantileEstimator = sortedIndicesValues[0].second;

    if(_quantileEstimator < 1e-15) {
      //qDebug() << "mr = " << mr;
      //qDebug() << "Sorted indices values (using 0):";
      for(auto pair: sortedIndicesValues) {
        auto attrVals = _clusters->at(pair.first)->getRepresentative()->attributesValues;
        //std::vector<double> pt = {std::stod(attrVals["Val0"]), std::stod(attrVals["Val1"])};
        std::vector<double> pt = {std::stod(attrVals["Val0"])};
        //qDebug() << "\ti: " << pair.first << ", x: " << pt[0] << ", y: " << pt[1] << ", remembered value: " << pair.second << ", estimator value: " << _estimator->getValue(&pt);
      }
    }

    _estimator->_shouldConsiderWeights = originalShouldConsiderWeights;
    return;
  }

  int i = mr + 0.5;

  _quantileEstimator =
      (0.5 + i - mr) * sortedIndicesValues[i - 1].second; // Remember that indices in the formulas start from 1.
  _quantileEstimator +=
      (0.5 - i + mr) * sortedIndicesValues[i].second; // Remember that indices in the formulas start from 1.

  if(_quantileEstimator < 1e-15) {
    qDebug() << "mr = " << mr;
    qDebug() << "Sorted indices values (using " << i - 1 << "and" << i << "):";
    for(auto pair: sortedIndicesValues) {
      auto attrVals = _clusters->at(pair.first)->getRepresentative()->attributesValues;
      //std::vector<double> pt = {std::stod(attrVals["Val0"]), std::stod(attrVals["Val1"])};
      std::vector<double> pt = {std::stod(attrVals["Val0"])};
      qDebug() << "\ti: " << pair.first << ", x: " << pt[0] << ", y: " << pt[1]
               << ", remembered value: " << pair.second
               << ", estimator value: " << _estimator->getValue(&pt);
    }
  }
  _estimator->_shouldConsiderWeights = originalShouldConsiderWeights;
}

std::vector<double> DESDA::getRareElementsEnhancedKDEValues(const std::vector<std::vector<double>> *X, int dimension) {
  std::vector<double> enhancedKDEValues = {};
  auto currentClusters = getClustersForEstimator();
  auto standardWeights = getClustersWeights(currentClusters);

  sigmoidallyEnhanceClustersWeights(&currentClusters);
  enhanceWeightsOfUncommonElements();

  _enhancedKDE->setClusters(currentClusters);
  _enhancedKDE->setSmoothingParameters({_smoothingParametersVector});

  _examinedClustersW.clear();

  for(auto index : _examinedClustersIndices) {
    if(index < 0)
      _examinedClustersW.push_back(0);
    else
      _examinedClustersW.push_back(currentClusters[index]->getCWeight());
  }

  std::vector<double> attributesValues =
      getAttributesValuesFromClusters(currentClusters, dimension);
  double domainMinValue =
      getDomainMinValue(attributesValues, _windowedSmoothingParametersVector[0]);
  double domainMaxValue =
      getDomainMaxValue(attributesValues, _windowedSmoothingParametersVector[0]);

  for(auto x: *X) {
    if(x[dimension] > domainMinValue && x[dimension] < domainMaxValue) {
      enhancedKDEValues.push_back(_enhancedKDE->getValue(&x));
    }
    else {
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
QVector<std::pair<std::vector<double>, double>> DESDA::getAtypicalElementsValuesAndDerivatives() {
  QVector<std::pair<std::vector<double>, double>> atypicalElementsValuesAndDerivatives = {};
  auto atypicalElements = getAtypicalElements();


  /*
  qDebug() << "ATYPICAL ELEMENTS KDE VALUES";
  for(auto v : atypicalElements){
    qDebug() << "\t" << v->_currentKDEValue;
  }
  //*/

  int negative_derivatives = 0;
  int positive_derivatives = 0;

  for(auto a : atypicalElements) {
    std::pair<std::vector<double>, double> point_derivative = std::pair<std::vector<double>, double>({}, 0);

    for(size_t i = 0; i < this->stationarityTests.size(); ++i) {
      std::string attribute_name = "Val" + std::to_string(i);
      point_derivative.first.push_back(std::stod(a->getObject()->attributesValues[attribute_name]));
    }

    point_derivative.second = a->_currentDerivativeValue;
    atypicalElementsValuesAndDerivatives.push_back(point_derivative);

    // DEBUG
    /*
    if(a->_currentDerivativeValue > 0){
      ++positive_derivatives;
    } else {
      ++negative_derivatives;
    }
    // DEBUG
    //*/
  }

  return atypicalElementsValuesAndDerivatives;
}

vector<double> DESDA::GetPrognosisErrors() {

  int last_examined_cluster_number = std::min(int(_clustersForWindowed.size()), max_prognosis_error_clusters_);

  vector<double> errors;

  for(int i = 1; i < last_examined_cluster_number; ++i){
    errors.push_back((_clustersForWindowed)[i]->getLastPrediction() - (_clustersForWindowed)[i]->_currentKDEValue);
  }

  return errors;
}

double DESDA::ComputePrognosisError(const vector<double> &errors) const {
  return average(errors);
}

double DESDA::ComputeStatistics(const std::vector<double> &errors) const {
  if(errors.size() > 1){
    // DEBUG //
      // qDebug() << "Err:" << errors;
      // qDebug() << "Avg: " << average(errors);
      // qDebug() << "Std: " << stdev(errors);
    // DEBUG //
    return average(errors) / stdev(errors);
  }
  return 0;
}




