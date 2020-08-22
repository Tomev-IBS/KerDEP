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

#include "DESDA.h"

#include <QtDataVisualization>
#include "UI/plot.h"

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
    QVector<double> _windowedErrorDomain;
    QVector<double> _errorDomain;

  private:
    const QPen _MODEL_PLOT_PEN                     = QPen(Qt::red);
    const QPen _WINDOWED_PLOT_PEN                  = QPen(Qt::black);
    const QPen _KDE_PLOT_PEN                       = QPen(Qt::blue);
    const QPen _WEIGHTED_PLOT_PEN                  = QPen(Qt::cyan);
    const QPen _DERIVATIVE_PLOT_PEN                = QPen(QColor(255, 165, 0)); // Orange
    const QPen _STANDARIZED_DERIVATIVE_PLOT_PEN    = QPen(QColor(255, 220, 0)); // Yellow
    const QPen _DESDA_KDE_PLOT_PEN                 = QPen(Qt::magenta);
    const QPen _DESDA_RARE_ELEMENTS_KDE_PLOT_PEN   = QPen(Qt::green);

    // 3D plots
    //QtDataVisualization::Q3DSurface surface;

    // Contour plots
    Plot *contourPlot = nullptr;

    std::vector<QCPAbstractItem *> _linesOnPlot;

    void setupValidators();
    void setupPlot();
    void setupKernelsTable();

    long long start;

    double _longestStepExecutionInSecs = 0;

    int screenGenerationFrequency = 1;

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
    QVector<double> _modelPlotY;
    QVector<double> _sigmoidallyEnhancedPlotY;
    QVector<double> _rareElementsEnhancedPlotY;
    QVector<double> _lessElementsEstimatorY;
    QVector<double> _weightedEstimatorY;
    QVector<double> _windowedEstimatorY;

    QVector<double> _windowedModelPlotY;
    QVector<double> _windowedEstimatorErrorY;
    QVector<double> _modelPlotErrorY;
    QVector<double> _lessElementsEstimatorErrorY;
    QVector<double> _weightedEstimatorErrorY;
    QVector<double> _sigmoidallyEnhancedErrorPlotY;
    QVector<double> _rareElementsEnhancedErrorPlotY;

    QVector<double> getTargetFunctionValuesOnDomain(QVector<double> *domain);


    // Errors
    double _L1_w = 0, _L1_m = 0, _L1_d = 0, _L1_p = 0, _L1_n = 0,
    _L2_w = 0, _L2_m = 0, _L2_d = 0, _L2_p = 0, _L2_n = 0,
    _sup_w = 0, _sup_m = 0, _sup_d = 0, _sup_p = 0, _sup_n = 0,
    _mod_w = 0, _mod_m = 0, _mod_d = 0, _mod_p = 0, _mod_n = 0;

    double  modelExtrema = 0.0, windowKDEExtrema = 0.0, KDEExtrema = 0.0, KDEPExtrema = 0.0, WKDEExtrema = 0.0, REESEExtrema = 0,
            lastModelExtrema = 0.0, lastKDEExtrema = 0.0, lastKDEPExtrema = 0.0, lastWKDEExtrema = 0.0;

    double calculateL1Error(const QVector<double> &model, const QVector<double> estimated, const double &domainLength);
    double calculateL2Error(const QVector<double> &model, const QVector<double> estimated, const double &domainLength);
    double calculateSupError(const QVector<double> &model, const QVector<double> estimated);
    double findExtrema(const QVector<double> &values, const QVector<double> domain);
    point find2DExtrema(const QVector<double> &values, const QVector<point> &points);

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
    std::shared_ptr<function> targetFunction;

    void drawPlots(DESDA* DESDAAlgorithm);
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
    void on_pushButton_clicked();

    void resizeEvent(QResizeEvent* event);

    // 2D Plot
    QVector<std::vector<double>> generate2DPlotErrorDomain(DESDA *DESDAAlgorithm);
    double calculate2DDomainArea(const QVector<std::vector<double>>& domain);
    QVector<double> getFunctionsValueOnDomain(function *func, QVector<point> domain);

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
