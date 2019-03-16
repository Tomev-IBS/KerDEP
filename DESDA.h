#ifndef DESDA_H
#define DESDA_H

#include <QDebug>
#include <memory>

#include "KDE/kerneldensityestimator.h"
#include "KDE/weightedSilvermanSmoothingParameterCounter.h"
#include "Distributions/distributions.h"
#include "Reservoir_sampling/reservoirSamplingAlgorithm.h"
#include "groupingThread/groupingThread.h"

class DESDA
{
  public:
    DESDA(std::shared_ptr<kernelDensityEstimator> estimator,
          std::shared_ptr<kernelDensityEstimator> estimatorDerivative,
          std::shared_ptr<kernelDensityEstimator> enchancedKDE,
          double weightModifier,
          weightedSilvermanSmoothingParameterCounter *_smoothingParamCounter,
          reservoirSamplingAlgorithm *samplingAlgorithm,
          std::vector<std::shared_ptr<cluster>> *clusters,
          std::vector<std::vector<std::shared_ptr<cluster>>> *storedMedoids,
          double desiredRarity, groupingThread *gt,
          double v, double newWeightB, int mE);

    void performStep();
    QVector<double> getKernelPrognosisDerivativeValues(const QVector<qreal> *X);
    QVector<double> getEnhancedKDEValues(const QVector<qreal> *X);
    double getAverageOfFirstMSampleValues(int M);
    double getStdDevOfFirstMSampleValues(int M);
    cluster getEmECluster();

    double emEStDev();
    double aemEAvg();
    double aemEStDev();
    double aemEVersor();

    double _u_i = 0.0;
    std::vector<double> _selectedVValues = {};

    double beta = 20000.0, alpha = 0.0005, delta = 0.5, gamma = 50000.0;
    double w_E = 0.98;

    std::vector<double> aemEVals = {};
    std::vector<double> emEVals = {};

  protected:

    const double _MAX_A = 1.5;
    const double _MIN_A = 0.01;

    int _stepNumber = 0;
    int _numberOfClustersForGrouping = 100;
    int _medoidsNumber = 50;
    int _mE = 0; // Cardinality of observed new objects

    double _weightModifier = 0.0;
    double _smoothingParameterMultiplier = 1.0;
    double _positionalSecondGradeEstimator = 0.0;
    double _maxEstimatorValueOnDomain = 0.0;
    double _a = 0.0;
    double _desiredRarity = 0.01;
    double _previousUncommonClustersWeight = 0.0;
    double _v = 1.0;

    double _newWeightB = 0;

    cluster emE;
    bool _shouldCluster = false;

    std::shared_ptr<kernelDensityEstimator> _estimator;
    std::shared_ptr<kernelDensityEstimator> _estimatorDerivative;
    std::shared_ptr<kernelDensityEstimator> _enhancedKDE;
    std::shared_ptr<distribution> _targetDistribution;
    weightedSilvermanSmoothingParameterCounter *_smoothingParamCounter;
    reservoirSamplingAlgorithm *_samplingAlgorithm;
    std::vector<std::shared_ptr<cluster>> *_clusters;
    std::vector<std::shared_ptr<sample>> _objects;

    std::vector<std::vector<std::shared_ptr<cluster>>> *_storedMedoids;
    std::vector<std::shared_ptr<cluster>> _uncommonClusters;

    groupingThread *_grpThread;

    void updateWeights();
    std::vector<std::shared_ptr<cluster>> getClustersForEstimator();
    void countKDEValuesOnClusters();
    unsigned long long findUncommonClusters();
    void updateA();
    void updatePrognosisParameters();
    void countKDEDerivativeValuesOnClusters();

};

#endif // DESDA_H
