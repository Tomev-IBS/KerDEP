#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleValidator>
#include <QStringList>
#include <vector>
#include <unordered_map>
#include <chrono>

#include "QCustomPlot/qcustomplot.h"

#include "Reservoir_sampling/reservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/sample.h"
#include "Functions/Kernels/kernels.h"
#include "KDE/kerneldensityestimator.h"
#include "KDE/smoothingParameterCounter.h"
#include "Functions/function.h"
#include "groupingThread/kMedoidsAlgorithm/attributeData.h"

enum positionalSecondGradeEstimatorCountingMethods
{
  STANDARD = 0,
  WEIGHTED = 1
};

namespace Ui
{
    class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

  public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

  protected:
    void keyPressEvent(QKeyEvent* event);

    int stepNumber = 0;

    std::vector<std::shared_ptr<cluster>> uncommonClusters;

  private:
    void setupValidators();
    void setupPlot();
    void setupKernelsTable();

    long long start;

    double _longestStepExecutionInSecs = 0;

    // Uncommon clusters dynamic parameter
    const double MAX_A                = 1.5;
    const double MIN_A                = 0.01;
    double _a                         = 1;
    double _maxEstimatorValueOnDomain = 0;
    double _previousUncommonClustersWeight = 0.0;
    void updateA();

    // Errors
    double _errorEJ = 0.0, _errorEJP = 0.0;
    std::vector<qreal> KDEValues, lastKDEValues, KDEPValues, ModelValues;
    QVector<qreal> WKDEValues, lastWKDEValues;

    // Default settings
    const qreal MAX_X                 = 999.0;
    const qreal MIN_X                 = -999.0;
    const qreal MAX_Y                 = 99.0;
    const qreal MIN_Y                 = -99.0;
    const qreal MIN_SMOOTHING_P       = 0.0;
    const qreal MAX_SMOOTHING_P       = 2.0;
    const int   DECIMAL_NUMBERS       = 3;
    const qreal DEFAULT_MIN_X         = -5;
    const qreal DEFAULT_MAX_X         = 15;
    const qreal DEFAULT_MIN_Y         = -0.3;
    const qreal DEFAULT_MAX_Y         = 0.5;

    const unsigned int MEDOIDS_NUMBER = 10;

    std::vector<std::vector<std::shared_ptr<cluster>>> storedMedoids;

    Ui::MainWindow *ui;

    QVector<std::shared_ptr<QVector<qreal>>> samples;
    std::vector<std::shared_ptr<sample>> objects;
    std::vector<std::shared_ptr<cluster>> *clusters;

    QVector<qreal> KDETemporalDerivativeY;

    std::unordered_map<std::string, attributeData*> attributesData;

    // Tests
    void testNewFunctionalities();

    QVector<qreal> KDEEstimationY;

    // Prediction
    QVector<double> _kernelPrognosisDerivativeValues;

    // Enhancement
    QVector<double> _sigmoidallyEnhancedPlotY;
    QVector<double> _lastSigmoidallyEnhancedPlotY;
    double  _summaricKDEError1 = 0.0, _summaricKDEPError1 = 0.0, _summaricKDESError1 = 0.0,
            _summaricKDEError2 = 0.0, _summaricKDEPError2 = 0.0, _summaricKDESError2 = 0.0,
            _summaricKDEErrorSup = 0.0, _summaricKDEPErrorSup = 0.0, _summaricKDESErrorSup = 0.0,
            _summaricKDEErrorMod = 0.0, _summaricKDEPErrorMod = 0.0, _summaricKDESErrorMod = 0.0,
            _summaricKDEsError1 = 0.0, _summaricKDEPsError1 = 0.0, _summaricKDESsError1 = 0.0,
            _summaricKDEsError2 = 0.0, _summaricKDEPsError2 = 0.0,_summaricKDESsError2 = 0.0,
            _summaricKDEsErrorSup = 0.0, _summaricKDEPsErrorSup = 0.0,
            _summaricKDEsErrorMod = 0.0, _summaricKDEPsErrorMod = 0.0,
            _summaricKDESsErrorMod = 0.0, _summaricKDESsErrorSup = 0.0;

    double  modelExtrema = 0.0, KDEExtrema = 0.0, KDEPExtrema = 0.0, WKDEExtrema = 0.0,
            lastModelExtrema = 0.0, lastKDEExtrema = 0.0, lastKDEPExtrema = 0.0, lastWKDEExtrema = 0.0;

    QVector<double> maxAs = {};

    double positionalSecondGradeEstimator = 0.0;

    std::shared_ptr<dataParser> parser;
    std::shared_ptr<dataReader> reader;

    QVector<qreal> oldKerernelY;
    QVector<qreal> newKernelY;

    std::shared_ptr<kernelDensityEstimator> kernelPrognoser;
    std::shared_ptr<kernelDensityEstimator> _enchancedKDE;

    QVector<std::shared_ptr<QVector<qreal>>> means, stDevs;

    int positionalSecondGradeEstimatorCountingMethod = WEIGHTED;

    QStringList kernelTypes;

    unsigned long long insertObjectsBetweenIntervals(unsigned int objectsNumber);
    unsigned long long generateInterIntervalObjects(std::vector<std::shared_ptr<sample>> *interIntervalObjects,
                                                    unsigned int objectsNumber);
    unsigned long long selectDesiredNumberOfInterIntervalObjects(std::vector<std::shared_ptr<sample>> *interIntervalObjects);
    unsigned long long insertClustersFromInterIntervalObjects(std::vector<std::shared_ptr<sample>> *interIntervalObjects);
    double setInterIntervalClustersWeights(std::vector<std::shared_ptr<cluster>> *newClusters);
    double countInterIntervalClustersWeight();

    unsigned long long insertMassiveData();
    unsigned long long generateMassiveData(std::vector<std::shared_ptr<sample> > *dataContainer);

    void drawPlots(kernelDensityEstimator* estimator, function* targetFunction);
    void clearPlot();
    void resizePlot();
    void addModelPlot(const QVector<qreal> *X, const QVector<qreal> *Y);
    void addEstimatedPlot(const QVector<qreal> *X, const QVector<qreal> *Y);
    void addWeightedEstimatorPlot(const QVector<qreal> *X, const QVector<qreal> *Y);
    double countNewtonianDerivative(int i, const QVector<qreal> *Y);
    void addKernelPrognosisDerivativePlot(const QVector<qreal> *X);
    void addSigmoidallyEnhancedEstimationPlot(const QVector<qreal> *X, kernelDensityEstimator *estimator);
    int updateClusterPredictionParameter(std::shared_ptr<cluster> c, double KDEValue);
    int initializeClusterPredictionParameter(std::shared_ptr<cluster> c, double KDEValue);
    unsigned long long markUncommonClusters();
    void markNewTrends();
    void markClustersWithNegativeDerivative();

    void addTemporalDerivativePlot(const QVector<qreal> *X, const QVector<qreal> *Y);

    void fillStandardDeviations(QVector<std::shared_ptr<QVector<qreal> > > *stDevs);
    void fillMeans(QVector<std::shared_ptr<QVector<qreal> > > *means);

  private slots:

    void refreshKernelsTable();
    void addKernelToTable(int rowNumber, QDoubleValidator *smoothingParameterValidator);
    void refreshTargetFunctionTable();
    void uniformContributions();
    qreal countLastContribution();
    void updateLastContribution();

    void fillDomain(QVector<std::shared_ptr<point> > *domain, std::shared_ptr<point> *prototypePoint);
    distribution *generateTargetDistribution(QVector<std::shared_ptr<QVector<qreal> > > *means,
                                             QVector<std::shared_ptr<QVector<qreal> > > *stDevs);
    reservoirSamplingAlgorithm *generateReservoirSamplingAlgorithm(dataReader *reader,
                                                                   dataParser *parser);
    kernelDensityEstimator *generateKernelDensityEstimator(int dimensionsNumber);
    function *generateTargetFunction(QVector<std::shared_ptr<QVector<qreal> > > *means, QVector<std::shared_ptr<QVector<qreal> > > *stDevs);

    void on_pushButton_start_clicked();
    int canAnimationBePerformed(int dimensionsNumber);
    void clusterMassiveData(std::vector<std::shared_ptr<sample>> *objects,
                            std::vector<std::vector<std::shared_ptr<cluster>>> *storage);
    std::vector<std::shared_ptr<cluster>> getClustersForEstimator();
    void countKDEValuesOnClusters(std::shared_ptr<kernelDensityEstimator> estimator);
    unsigned long long findUncommonClusters();

    QString formatNumberForDisplay(double number);


    void delay(int ms);

    void on_spinBox_dimensionsNumber_editingFinished();

    void on_pushButton_addTargetFunction_clicked();

    void on_pushButton_removeTargetFunction_clicked();

    double numericIntegral(const QVector<qreal> *Y);
};

enum kernelSettingsColumns
{
  KERNEL_COLUMN_INDEX                 = 0,
  SMOOTHING_PARAMETER_COLUMN_INDEX    = 1,
  CARRIER_RESTRICTION_COLUMN_INDEX    = 2
};

enum targetFunctionSettingsColumns
{
  MEAN_COLUMN_INDEX           = 0,
  STDEV_COLUMN_INDEX          = 1,
  CONTRIBUTION_COLUMN_INDEX   = 2
};

enum reservoirSamplingAlgorithms
{
  BASIC_RESERVOIR_SAMPLING_ALGORITHM = 0,
  BIASED_RESERVOIR_SAMPLING_ALGORITHM = 1
};

#endif // MAINWINDOW_H
