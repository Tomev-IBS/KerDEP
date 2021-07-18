#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QMainWindow>
#include <QDoubleValidator>
#include <QStringList>
#include <vector>
#include <unordered_map>
#include <chrono>

#include "UI/plotLabel.h"

#include "QCustomPlot/qcustomplot.h"

#include "Reservoir_sampling/reservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/sample.h"
#include "Functions/Kernels/kernels.h"
#include "KDE/kerneldensityestimator.h"
#include "KDE/smoothingParameterCounter.h"
#include "Functions/function.h"
#include "groupingThread/kMedoidsAlgorithm/attributeData.h"
#include "ClusterKernelWrappers/enhancedClusterKernelAlgorithm.h"

#include "Benchmarking/errorsCalculator.h"

#include "DESDA.h"
#include "kerDepCcWde.h"
#include "SOMKEAlgorithm.h"

#include "UI/plot.h"

enum class KernelSettingsColumns : int {
  kKernelColumnIndex = 0,
  kSmoothingParameterColumnIndex = 1,
  kCarrierRestrictionColumnIndex = 2
};

enum class TargetFunctionSettingsColumns : int {
  kMeanColumnIndex = 0,
  kStandardDeviationColumnIndex = 1,
  kContributionColumnIndex = 2
};

enum class ReservoirSamplingAlgorithms : int {
  kBasicReservoirSamplingAlgorithm = 0,
  kBiasedReservoirSamplingAlgorithm = 1
};

namespace Ui {
  class MainWindow;
}

typedef std::shared_ptr<double> double_ptr;
typedef std::shared_ptr<int> int_ptr;

class MainWindow : public QMainWindow {
  Q_OBJECT

  public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;

  protected:
    int step_number_ = 0;
    QVector<std::shared_ptr<point>> domain_;
    QVector<double> drawable_domain_;

    QVector<std::shared_ptr<plotLabel>> plot_labels_;
    double label_vertical_offset_ = 0;
    double label_horizontal_offset_ = 0;
    const double label_vertical_offset_step_ = 0.03;

    void AddErrorLabelsToPlot(const QVector<QString> &labels, const QVector<double_ptr> &values);
    void AddErrorLabelToPlot(const QString &label, double *value);
    void AddDoubleLabelToPlot(const QString &label, double *value);
    void AddIntLabelToPlot(const QString &label, int *value);
    void AddConstantLabelToPlot(const QString &label);
    void AddL1ErrorsToSum(QVector<ErrorsCalculator*> &errors_calculators, QVector<double> &errors_sums);
    void AddL2ErrorsToSum(QVector<ErrorsCalculator*> &errors_calculators, QVector<double> &errors_sums);
    void AddSupErrorsToSum(QVector<ErrorsCalculator*> &errors_calculators, QVector<double> &errors_sums);
    void AddModErrorsToSum(QVector<ErrorsCalculator*> &errors_calculators, QVector<double> &errors_sums);
    void AddColorsLegendToPlot();


  private:
    // Pens for 1d plot
    const QPen model_plot_pen_ = QPen(Qt::red);
    //const QPen windowed_plot_pen_ = QPen(QColor(255, 220, 0));
    const QPen windowed_plot_pen_ = QPen(QColor(155, 135, 11));
    const QPen kde_plot_pen_ = QPen(QColor(0, 255, 0));
    const QPen weighted_plot_pen_ = QPen(QColor(0, 255, 255));
    const QPen desda_kde_plot_pen_ = QPen(QColor(0, 0, 255));
    const QPen desda_rare_elements_kde_plot_pen_ = QPen(Qt::black, 2);

    const QPen derivative_plot_pen_ = QPen(QColor(255, 165, 0)); // Orange
    const QPen standardized_derivative_plot_pen_ = QPen(QColor(115, 65, 45)); // Yellow

    void FillErrorIndicesColors();

    QVector<QColor> error_indices_colors = {};

    // Contour plots
    Plot *contour_plot_ = nullptr;
    std::vector<QCPAbstractItem *> lines_on_plot_;

    long long start = 0;
    int screen_generation_frequency_ = 1;
    // Default settings
    const qreal kMaxX = 999.0;
    const qreal kMinX = -999.0;
    const qreal kMaxY = 99.0;
    const qreal kMinY = -99.0;
    const qreal kMinSmoothingParameter = 0.0;
    const qreal kMaxSmoothingParameter = 2.0;
    const int kDecimalNumbers = 3;
    const qreal kDefaultMinX = -5;
    const qreal kDefaultMaxX = 15;
    const qreal kDefaultMinY = -0.3;
    const qreal kDefaultMaxY = 0.5;
    std::vector<std::shared_ptr<cluster>> stored_medoids_;
    Ui::MainWindow *ui;
    vector<std::shared_ptr<vector<double>>> samples_;
    std::vector<std::shared_ptr<sample>> objects_;
    std::vector<std::shared_ptr<cluster>> *clusters_ = nullptr;
    std::unordered_map<std::string, attributeData *> attributes_data_;
    // Tests
    /*static*/ void testNewFunctionalities();
    void SetupValidators();
    void SetupPlot();
    void SetupKernelsTable();
    // Prediction
    QVector<double> kernel_prognosis_derivative_values_;
    // Errors
    double l1_w_ = 0, l1_n_ = 0,
           l2_w_ = 0, l2_n_ = 0,
           sup_w_ = 0, sup_n_ = 0,
           mod_w_ = 0, mod_n_ = 0;
    QVector<std::pair<double, double>>
        atypical_elements_values_and_derivatives_ = {};
    double quantile_estimator_value_ = 0;
    std::shared_ptr<dataParser> parser_;
    std::shared_ptr<dataReader> reader_;
    std::shared_ptr<kernelDensityEstimator> derivative_estimator_;
    std::shared_ptr<kernelDensityEstimator> enhanced_kde_;
    vector<std::shared_ptr<vector<double>>> means_, standard_deviations_;
    QStringList kernel_types_;
    std::shared_ptr<function> target_function_;
    void DrawPlots(DESDA *DESDAAlgorithm);
    void DrawPlots(EnhancedClusterKernelAlgorithm *CKAlgorithm);
    void DrawPlots(KerDEP_CC_WDE *WDEAlgorithm);
    void DrawPlots(SOMKEAlgorithm *somke_algorithm);
    void ClearPlot();
    void AddPlot(const QVector<qreal> *Y, const QPen &pen);
    void ResizePlot();
    unsigned long long MarkUncommonClusters();
    void FillStandardDeviations(
        vector<std::shared_ptr<vector<double>>> *stDevs);
    void FillMeans(vector<std::shared_ptr<vector<double>>> *means);

  private slots:

    void RefreshKernelsTable();
    void AddKernelToTable(int rowNumber,
                          QDoubleValidator *smoothingParameterValidator);
    void RefreshTargetFunctionTable();
    void UniformContributions();
    qreal CountLastContribution();
    void UpdateLastContribution();
    void FillDomain(QVector<std::shared_ptr<point> > *domain,
                    std::shared_ptr<point> *prototypePoint);
    distribution *GenerateTargetDistribution(
        vector<std::shared_ptr<vector<double>>> *means,
        vector<std::shared_ptr<vector<double>>> *stDevs);
    reservoirSamplingAlgorithm *GenerateReservoirSamplingAlgorithm(
        dataReader *reader,
        dataParser *parser);
    kernelDensityEstimator *GenerateKernelDensityEstimator(
        int dimensionsNumber);
    function *GenerateTargetFunction(
        vector<std::shared_ptr<vector<double>>> *means,
        vector<std::shared_ptr<vector<double>>> *stDevs);
    static int CanAnimationBePerformed(int dimensionsNumber);
    static QString FormatNumberForDisplay(double number);
    void on_pushButton_start_clicked();
    void Run1DExperimentWithDESDA();
    void Run1DExperimentWithClusterKernels();
    void Run1DExperimentWithWDE();
    void Run1DExperimentWithSOMKE();
    void on_spinBox_dimensionsNumber_editingFinished();
    void on_pushButton_addTargetFunction_clicked();
    void on_pushButton_removeTargetFunction_clicked();
    void on_pushButton_clicked();
    void resizeEvent(QResizeEvent *event) override;
    // 2D Plot
    static std::vector<std::vector<double>> Generate2DPlotErrorDomain(DESDA *DESDAAlgorithm);
    static std::vector<std::vector<double>> Generate1DPlotErrorDomain(DESDA *DESDAAlgorithm);
    static std::vector<std::vector<double>> Generate1DWindowedPlotErrorDomain(DESDA *DESDAAlgorithm);
    static double Calculate2DDomainArea(const std::vector<std::vector<double>> &domain);
    static std::vector<double> GetFunctionsValueOnDomain(function *func, const std::vector<std::vector<double>> &domain);

};



#endif // MAIN_WINDOW_H
