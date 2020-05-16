#ifndef DESDA_H
#define DESDA_H

#include <QDebug>
#include <memory>
#include <map>

#include "KDE/kerneldensityestimator.h"
#include "KDE/weightedSilvermanSmoothingParameterCounter.h"
#include "Distributions/distributions.h"
#include "Reservoir_sampling/reservoirSamplingAlgorithm.h"
#include "groupingThread/groupingThread.h"
#include "StationarityTests/kpssstationaritytest.h"

class DESDA
{
  public:
    DESDA(std::shared_ptr<kernelDensityEstimator> estimator,
          std::shared_ptr<kernelDensityEstimator> estimatorDerivative,
          std::shared_ptr<kernelDensityEstimator> enchancedKDE,
          double weightModifier,
          reservoirSamplingAlgorithm *samplingAlgorithm,
          std::vector<std::shared_ptr<cluster>> *clusters,
          std::vector<std::shared_ptr<cluster>> *storedMedoids,
          double desiredRarity, groupingThread *gt, double newWeightB);

    void performStep();
    QVector<double> getKernelPrognosisDerivativeValues(const QVector<qreal> *X);
    QVector<double> getEnhancedKDEValues(const QVector<qreal> *X);
    std::vector<double> getClustersWeights(const std::vector<std::shared_ptr<cluster>> &clusters);
    void sigmoidallyEnhanceClustersWeights(std::vector<std::shared_ptr<cluster>> *clusters);
    QVector<double> getWindowKDEValues(const QVector<qreal> *X);
    QVector<double> getKDEValues(const QVector<qreal> *X);
    QVector<double> getWeightedKDEValues(const QVector<qreal> *X);
    double getAverageOfFirstMSampleValues(int M);
    double getStdDevOfFirstMSampleValues(int M);
    cluster getEmECluster();
    double getStationarityTestValue();

    stationarityTestPtr stationarityTest;

    double _u_i = 0.0;

    // Prediction
    double _avg = 0;
    double _stDev = 1;

    // New prediction
    double _beta0 = 0.7;
    double _v = 0.999; // deactualization w of clusters

    // Sizes
    int _maxM = 1000;
    int _minM = 200;
    int _m = 0; // Cardinality of objects to build KDE

    double _newWeightB = 0;
    int _kpssM = 0; // m_Eta

    // From 31 XI 2019 mail
    // Delta update parameters
    double _s = -4.0;
    double _mu = 1000;

    // Rare elements
    double _r = 0.05;
    double _quantileEstimator = 0;
    int _rareElementsNumber = 0;
    QVector<double> getRareElementsEnhancedKDEValues(const QVector<qreal> *X);
    std::vector<std::shared_ptr<cluster>> getAtypicalElements();
    QVector<std::pair<double, double> > getAtypicalElementsValuesAndDerivatives();

    // General purpose
    int _sgmKPSSPercent = 25;
    double _sgmKPSS = 0;

    // Analysis
    double getMaxAbsAOnLastKPSSMSteps();

    /** According to 2020 first article **/
    std::vector<int> _examinedClustersIndices = {};
    std::vector<double> _maxAbsAs = {};
    std::vector<double> _maxAbsDerivatives = {};
    std::vector<double> _examinedClustersAs = {};
    std::vector<double> _examinedClustersDerivatives = {};
    std::vector<double> _examinedClustersWStar = {};
    std::vector<double> _examinedClustersWStar2 = {};
    std::vector<double> _examinedClustersWStar3 = {};
    std::vector<double> _examinedClustersW = {};

    QVector<double> getErrorDomain();
    QVector<double> getWindowedErrorDomain();

    double _averageMaxDerivativeValueInLastMASteps = 0;
    double _maxAbsDerivativeValueInCurrentStep = 0;
    int _mA = 100;

    double _h;
    double _hWindowed;

  protected:

    int _stepNumber = 0;
    int _numberOfClustersForGrouping = 100;
    int _medoidsNumber = 50;

    double _weightModifier = 0.0;
    double _smoothingParameterMultiplier = 1.0;
    double _positionalSecondGradeEstimator = 0.0;
    double _maxEstimatorValueOnDomain = 0.0;
    double _a = 0.0;
    double _previousUncommonClustersWeight = 0.0;

    std::shared_ptr<kernelDensityEstimator> _estimator;
    std::shared_ptr<kernelDensityEstimator> _estimatorDerivative;
    std::shared_ptr<kernelDensityEstimator> _enhancedKDE;
    std::shared_ptr<distribution> _targetDistribution;
    weightedSilvermanSmoothingParameterCounter *_smoothingParamCounter;
    reservoirSamplingAlgorithm *_samplingAlgorithm;
    std::vector<std::shared_ptr<cluster>> *_clusters;
    std::vector<std::shared_ptr<sample>> _objects;

    std::vector<std::shared_ptr<cluster>> *_storedMedoids;
    std::vector<std::shared_ptr<cluster>> _uncommonClusters;

    groupingThread *_grpThread;

    // Random deletion
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dist;
    double _d = 0;

    int randomizeIndexToDelete();
    void updateWeights();
    void updateExaminedClustersIndices();
    std::vector<std::shared_ptr<cluster>> getClustersForEstimator();
    std::vector<std::shared_ptr<cluster>> getClustersForWindowedEstimator();
    void countKDEValuesOnClusters();
    void updatePrognosisParameters();
    void countDerivativeValuesOnClusters();
    void updateM();
    void updateDelta();
    void updateMaxAbsAVector();
    double getCurrentMaxAbsA();
    void updateMaxAbsDerivativeVector();
    double getCurrentMaxAbsDerivativeValue();
    void updateAverageMaxAbsDerivativeInLastMASteps();
    void updateMaxAbsDerivativeInCurrentStep();
    void updateExaminedClustersAsVector();

    // Domain reduction
    double getDomainMinValue(const std::vector<clusterPtr> &clusters, double h);
    double getDomainMaxValue(const std::vector<clusterPtr> &clusters, double h);

    double calculateH(const std::vector<clusterPtr> &clusters);

    void enhanceWeightsOfUncommonElements();
    std::vector<double> getVectorOfAcceleratedKDEValuesOnClusters();
    std::vector<std::pair<int, double> > getSortedAcceleratedKDEValues(const std::vector<double> &AKDEValues);
    void recountQuantileEstimatorValue(const std::vector<std::pair<int, double> > &sortedIndicesValues);

    // Only one of these parameters should be left, but I need this for
    // more accurate experiments preparation.
    std::map<int, std::vector<double>> _sgmKPSSParameters = {
      {10, {0.298, 2.417}}, {11, {0.329, 2.435}}, {12, {0.360, 2.466}},
      {13, {0.393, 2.487}}, {14, {0.423, 2.513}}, {15, {0.457, 2.535}},
      {16, {0.490, 2.559}}, {17, {0.523, 2.583}}, {18, {0.556, 2.608}},
      {19, {0.592, 2.634}}, {20, {0.625, 2.659}}, {21, {0.661, 2.685}},
      {22, {0.696, 2.711}}, {23, {0.732, 2.732}}, {24, {0.773, 2.768}},
      {25, {0.804, 2.781}}, {26, {0.842, 2.819}}, {27, {0.878, 2.846}},
      {28, {0.916, 2.874}}, {29, {0.954, 2.902}}, {30, {0.995, 2.932}}
    };

    //_sgmKPSS = sigmoid(0.995 * stationarityTest->getTestsValue() - 2.932); // 30 percent
    //_sgmKPSS = sigmoid(0.954 * stationarityTest->getTestsValue() - 2.902); // 29 percent
    //_sgmKPSS = sigmoid(0.916 * stationarityTest->getTestsValue() - 2.874); // 28 percent
};

#endif // DESDA_H
