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

#include "UI/plot.h"

enum class PositionalSecondGradeEstimatorCountingMethods {
  kStandard = 0,
  kWeighted = 1
};

enum class KernelSettingsColumns {
  kKernelColumnIndex = 0,
  kSmoothingParameterColumnIndex = 1,
  kCarrierRestrictionColumnIndex = 2
};

enum class TargetFunctionSettingsColumns {
  kMeanColumnIndex = 0,
  kStandardDeviationColumnIndex = 1,
  kContributionColumnIndex = 2
};

enum class ReservoirSamplingAlgorithms {
  kBasicReservoirSamplingAlgorithm = 0,
  kBiasedReservoirSamplingAlgorithm = 1
};

namespace Ui {
  class MainWindow;
}

class MainWindow : public QMainWindow {
  Q_OBJECT

  public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;

  protected:
    int step_number_ = 0;
    QVector<std::shared_ptr<point>> domain_;
    QVector<double> drawable_domain_;
    QVector<double> windowed_error_domain_;
    QVector<double> error_domain_;

  private:
    // Pens for 1d plot
    const static QPen model_plot_pen_;
    const static QPen windowed_plot_pen_;
    const static QPen kde_plot_pen_;
    const static QPen weighted_plot_pen_;
    const static QPen derivative_plot_pen_;
    const static QPen standardized_derivative_plot_pen_;
    const static QPen desda_kde_plot_pen_;
    const static QPen desda_rare_elements_kde_plot_pen;
    // Contour plots
    Plot *contourPlot = nullptr;
    std::vector<QCPAbstractItem *> _linesOnPlot;
    void setupValidators();
    void setupPlot();
    void setupKernelsTable();
    long long start;
    int screenGenerationFrequency = 1;
    // Default settings
    const qreal MAX_X = 999.0;
    const qreal MIN_X = -999.0;
    const qreal MAX_Y = 99.0;
    const qreal MIN_Y = -99.0;
    const qreal MIN_SMOOTHING_P = 0.0;
    const qreal MAX_SMOOTHING_P = 2.0;
    const int DECIMAL_NUMBERS = 3;
    const qreal DEFAULT_MIN_X = -5;
    const qreal DEFAULT_MAX_X = 15;
    const qreal DEFAULT_MIN_Y = -0.3;
    const qreal DEFAULT_MAX_Y = 0.5;
    std::vector<std::shared_ptr<cluster>> storedMedoids;
    Ui::MainWindow *ui;
    vector<std::shared_ptr<vector<double>>> samples;
    std::vector<std::shared_ptr<sample>> objects;
    std::vector<std::shared_ptr<cluster>> *clusters;
    std::unordered_map<std::string, attributeData *> attributesData;
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
    double calculateL1Error(const QVector<double> &model,
                            const QVector<double> estimated,
                            const double &domainLength);
    double calculateL2Error(const QVector<double> &model,
                            const QVector<double> estimated,
                            const double &domainLength);
    double calculateSupError(const QVector<double> &model,
                             const QVector<double> estimated);
    double findExtrema(const QVector<double> &values,
                       const QVector<double> domain);
    point find2DExtrema(const QVector<double> &values,
                        const QVector<point> &points);
    QVector<double> maxAs = {};
    QVector<std::pair<double, double>>
        _atypicalElementsValuesAndDerivatives = {};
    double _quantileEstimatorValue = 0;
    double positionalSecondGradeEstimator = 0.0;
    std::shared_ptr<dataParser> parser;
    std::shared_ptr<dataReader> reader;
    QVector<qreal> oldKernelY;
    std::shared_ptr<kernelDensityEstimator> kernelPrognoser;
    std::shared_ptr<kernelDensityEstimator> enhanced_KDE_;
    vector<std::shared_ptr<vector<double>>> means, stDevs;
    QStringList kernelTypes;
    std::shared_ptr<function> targetFunction;
    void drawPlots(DESDA *DESDAAlgorithm);
    void clearPlot();
    void addPlot(const QVector<qreal> *Y, const QPen &pen);
    void resizePlot();
    unsigned long long markUncommonClusters();
    void fillStandardDeviations(
        vector<std::shared_ptr<vector<double> > > *stDevs);
    void fillMeans(vector<std::shared_ptr<vector<double> > > *means);

  private slots:

    void refreshKernelsTable();
    void addKernelToTable(int rowNumber,
                          QDoubleValidator *smoothingParameterValidator);
    void refreshTargetFunctionTable();
    void uniformContributions();
    qreal countLastContribution();
    void updateLastContribution();
    void fillDomain(QVector<std::shared_ptr<point> > *domain,
                    std::shared_ptr<point> *prototypePoint);
    distribution *generateTargetDistribution(
        vector<std::shared_ptr<vector<double>>> *means,
        vector<std::shared_ptr<vector<double>>> *stDevs);
    reservoirSamplingAlgorithm *generateReservoirSamplingAlgorithm(
        dataReader *reader,
        dataParser *parser);
    kernelDensityEstimator *generateKernelDensityEstimator(
        int dimensionsNumber);
    function *generateTargetFunction(
        vector<std::shared_ptr<vector<double>>> *means,
        vector<std::shared_ptr<vector<double>>> *stDevs);
    void on_pushButton_start_clicked();
    int canAnimationBePerformed(int dimensionsNumber);
    QString formatNumberForDisplay(double number);
    void delay(int ms);
    void on_spinBox_dimensionsNumber_editingFinished();
    void on_pushButton_addTargetFunction_clicked();
    void on_pushButton_removeTargetFunction_clicked();
    void on_pushButton_clicked();
    void resizeEvent(QResizeEvent *event);
    // 2D Plot
    QVector<std::vector<double>> generate2DPlotErrorDomain(
        DESDA *DESDAAlgorithm);
    double calculate2DDomainArea(const QVector<std::vector<double>> &domain);
    QVector<double> getFunctionsValueOnDomain(function *func,
                                              QVector<point> domain);

};



#endif // MAINWINDOW_H
