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
    int stepNumber = 0;

    QVector<std::shared_ptr<point>> _domain;
    QVector<double> _drawableDomain;

  private:
    const QPen _MODEL_PLOT_PEN       = QPen(Qt::red);
    const QPen _WINDOWED_PLOT_PEN    = QPen(Qt::black);
    const QPen _KDE_PLOT_PEN         = QPen(Qt::blue);
    const QPen _WEIGHTED_PLOT_PEN    = QPen(Qt::cyan);
    const QPen _DERIVATIVE_PLOT_PEN  = QPen(Qt::green);
    const QPen _DESDA_KDE_PLOT_PEN   = QPen(Qt::darkRed);
    const QPen _DESDA_RARE_ELEMENTS_KDE_PLOT_PEN = QPen(Qt::magenta);

    std::vector<QCPAbstractItem *> _linesOnPlot;

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
    std::vector<qreal> ModelValues;

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

    std::vector<std::shared_ptr<cluster>> storedMedoids;

    Ui::MainWindow *ui;

    vector<std::shared_ptr<vector<double>>>  samples;
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
    QVector<double> _rareElementsEnhancedPlotY;
    QVector<double> _lessElementsEstimatorY;
    QVector<double> _weightedEstimatorY;
    QVector<double> _windowedEstimatorY;
    QVector<double> _lastSigmoidallyEnhancedPlotY;
    double _summaricWindowKDEError1 = 0, _summaricKDEError1 = 0, _summaricKDEPError1 = 0, _summaricKDESError1 = 0, _summaricKDENError1 = 0,
    _summaricWindowKDEError2 = 0, _summaricKDEError2 = 0, _summaricKDEPError2 = 0, _summaricKDESError2 = 0, _summaricKDENError2 = 0,
    _summaricWindowKDEErrorSup = 0, _summaricKDEErrorSup = 0, _summaricKDEPErrorSup = 0, _summaricKDESErrorSup = 0, _summaricKDENErrorSup = 0,
    _summaricWindowKDEErrorMod = 0, _summaricKDEErrorMod = 0, _summaricKDEPErrorMod = 0, _summaricKDESErrorMod = 0, _summaricKDENErrorMod = 0;

    double  modelExtrema = 0.0, windowKDEExtrema = 0.0, KDEExtrema = 0.0, KDEPExtrema = 0.0, WKDEExtrema = 0.0, REESEExtrema = 0,
            lastModelExtrema = 0.0, lastKDEExtrema = 0.0, lastKDEPExtrema = 0.0, lastWKDEExtrema = 0.0;

    QVector<double> maxAs = {};
    QVector<std::pair<double, double>> _atypicalElementsValuesAndDerivatives = {};
    double _quantileEstimatorValue = 0;

    double positionalSecondGradeEstimator = 0.0;

    std::shared_ptr<dataParser> parser;
    std::shared_ptr<dataReader> reader;

    QVector<qreal> oldKerernelY;
    QVector<qreal> newKernelY;

    std::shared_ptr<kernelDensityEstimator> kernelPrognoser;
    std::shared_ptr<kernelDensityEstimator> _enchancedKDE;

    vector<std::shared_ptr<vector<double>>> means, stDevs;

    int positionalSecondGradeEstimatorCountingMethod = WEIGHTED;

    QStringList kernelTypes;

    void drawPlots(kernelDensityEstimator* estimator, function* targetFunction);
    void clearPlot();
    void addPlot(const QVector<qreal> *Y, const QPen &pen);
    void resizePlot();
    void addModelPlot(const QVector<qreal> *X, const QVector<qreal> *Y);
    void addWindowedEstimatorPlot(const QVector<qreal> *X);
    void addEstimatedPlot(const QVector<qreal> *X, const QVector<qreal> *Y);
    void addWeightedEstimatorPlot(const QVector<qreal> *X, const QVector<qreal> *Y);
    void addKernelPrognosisDerivativePlot(const QVector<qreal> *X);
    void addSigmoidallyEnhancedEstimationPlot(const QVector<qreal> *X, kernelDensityEstimator *estimator);
    int updateClusterPredictionParameter(std::shared_ptr<cluster> c, double KDEValue);
    int initializeClusterPredictionParameter(std::shared_ptr<cluster> c, double KDEValue);
    unsigned long long markUncommonClusters();

    void addTemporalDerivativePlot(const QVector<qreal> *X, const QVector<qreal> *Y);

    void fillStandardDeviations(vector<std::shared_ptr<vector<double> > > *stDevs);
    void fillMeans(vector<std::shared_ptr<vector<double> > > *means);

  private slots:

    void refreshKernelsTable();
    void addKernelToTable(int rowNumber, QDoubleValidator *smoothingParameterValidator);
    void refreshTargetFunctionTable();
    void uniformContributions();
    qreal countLastContribution();
    void updateLastContribution();

    void fillDomain(QVector<std::shared_ptr<point> > *domain, std::shared_ptr<point> *prototypePoint);
    distribution* generateTargetDistribution(
      vector<std::shared_ptr<vector<double>>> *means,
      vector<std::shared_ptr<vector<double>>> *stDevs);
    reservoirSamplingAlgorithm *generateReservoirSamplingAlgorithm(dataReader *reader,
                                                                   dataParser *parser);
    kernelDensityEstimator *generateKernelDensityEstimator(int dimensionsNumber);
    function* generateTargetFunction(
        vector<std::shared_ptr<vector<double>>>* means,
        vector<std::shared_ptr<vector<double>>>* stDevs);

    void on_pushButton_start_clicked();
    int canAnimationBePerformed(int dimensionsNumber);

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
