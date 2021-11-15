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
          reservoirSamplingAlgorithm *samplingAlgorithm,
          std::vector<std::shared_ptr<cluster>> *clusters,
          double desiredRarity, double pluginRank=2);

    void performStep();
    QVector<double> getKernelPrognosisDerivativeValues(const QVector<qreal> *X, int dimension=0);
    std::vector<double> getEnhancedKDEValues(const std::vector<std::vector<double>> *X, int dimension=0);
    std::vector<double> getWeightedKDEValues(const vector<vector<double>> *X, int dimension= 0);
    std::vector<double> getClustersWeights(const std::vector<std::shared_ptr<cluster>> &clusters);
    void sigmoidallyEnhanceClustersWeights(std::vector<std::shared_ptr<cluster>> *clusters);
    vector<double> getWindowKDEValues(const vector<vector<double>> *X, int dimension=0);
    std::vector<double> getKDEValues(const std::vector<std::vector<double>> *X, int dimension=0);

    double getStationarityTestValue();

    // 2D Plot changes
    void prepareEstimatorForContourPlotDrawing();
    std::vector<double> _unmodifiedCWeightsOfClusters = {};
    void restoreClustersCWeights();

    std::vector<stationarityTestPtr> stationarityTests = {};

    // New prediction
    double _beta0 = 0.7;
    double _v = 0.999; // deactualization w of clusters

    // Sizes
    int _maxM = 1000;
    int _minM = 100;
    int _m = 0; // Cardinality of objects to build KDE

    int _kpssM = 0; // m_Eta

    // Rare elements
    double _r = 0.05;
    double _quantileEstimator = 0;
    int _rareElementsNumber = 0;
    int _trendsNumber = 0;
    std::vector<double> getRareElementsEnhancedKDEValues(const std::vector<vector<double>> *X, int dimension=0);
    std::vector<std::shared_ptr<cluster>> getAtypicalElements();
    QVector<std::pair<std::vector<double>, double> > getAtypicalElementsValuesAndDerivatives();

    // General purpose
    int _sgmKPSSPercent = 25;
    double _sgmKPSS = 0;

    /** According to 2020 first article **/
    std::vector<int> _examinedClustersIndices = {};
    std::vector<double> _examinedClustersDerivatives = {};
    std::vector<double> _examinedClustersWStar = {};
    std::vector<double> _examinedClustersWStar2 = {};
    std::vector<double> _examinedClustersWStar3 = {};
    std::vector<double> _examinedClustersW = {};

    QVector<double> getErrorDomain(int dimension=0);
    QVector<double> getWindowedErrorDomain(int dimension=0);

    double _maxAbsDerivativeValueInCurrentStep = 0;

    std::vector<double> _smoothingParametersVector;
    std::vector<double> _windowedSmoothingParametersVector;
    double _smoothingParameterEnhancer;

    // Prognosis error
    double e_ = 0;
    std::vector<double> statistics_ = {};
    double avg = 0;
    double std = 0;
    // DEBUG
      //QVector<double> vals;
      //QVector<double> progs;
      //QVector<double> stats;
      //QVector<double> errors;
    // DEBUG

    const double r010_ = 1.64;
    const double r005_ = 1.96;
    const double r001_ = 2.58;


    // Weighted C parameter
    bool compute_weighted_plugin = true;

    std::vector<std::shared_ptr<cluster>> getClustersForEstimator();

  protected:

    int _stepNumber = 0;
    int _pluginRank = 2;
    int error_domain_points_number_ = 500;

    std::shared_ptr<kernelDensityEstimator> _estimator;
    std::shared_ptr<kernelDensityEstimator> _estimatorDerivative;
    std::shared_ptr<kernelDensityEstimator> _enhancedKDE;
    reservoirSamplingAlgorithm *_samplingAlgorithm;
    std::vector<std::shared_ptr<cluster>> *_clusters;
    std::vector<std::shared_ptr<cluster>> _clustersForWindowed;
    std::vector<std::shared_ptr<sample>> _objects;

    void updateWeights();
    void updateExaminedClustersIndices();
    std::vector<std::shared_ptr<cluster>> getClustersForWindowedEstimator();
    void countKDEValuesOnClusters();
    void updatePrognosisParameters();
    void countDerivativeValuesOnClusters();
    void updateM();
    void updateMaxAbsDerivativeInCurrentStep();

    // Prognosis error
    std::vector<cluster> prognosis_clusters_ = {};
    cluster prognosis_cluster_;
    std::vector<std::vector<double>> prognosis_errors_ = {};
    int max_prognosis_error_clusters_ = 1;
    vector<double> GetPrognosisErrors();
    double ComputePrognosisError(const std::vector<double> &errors) const;
    double ComputeStatistics(const std::vector<double> &errors) const;

    // Domain reduction
    std::vector<double> getAttributesValuesFromClusters(std::vector<std::shared_ptr<cluster>> clusters, int dimension=0);
    double getDomainMinValue(const std::vector<double> &values, double h);
    double getDomainMaxValue(const std::vector<double> &values, double h);

    std::vector<double> calculateH(const std::vector<clusterPtr> &clusters);

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
      {28, {0.916, 2.874}}, {29, {0.954, 2.902}}, {30, {0.995, 2.932}},
      {31, {1.033, 2.960}}, {32, {1.074, 2.991}}, {33, {1.115, 3.021}},
      {34, {1.156, 3.051}}, {35, {1.196, 3.081}}, {36, {1.240, 3.113}},
      {37, {1.283, 3.145}}, {38, {1.327, 3.178}}, {39, {1.370, 3.209}},
      {40, {1.416, 3.243}}
    };
};

#endif // DESDA_H
