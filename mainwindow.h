#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleValidator>
#include <QStringList>
#include <vector>
#include <unordered_map>

#include "Reservoir_sampling/reservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/sample.h"
#include "Functions/Kernels/kernels.h"
#include "KDE/kerneldensityestimator.h"
#include "KDE/smoothingparametercounter.h"
#include "Functions/function.h"
#include "groupingThread/kMedoidsAlgorithm/attributeData.h"

namespace Ui
{
    class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

  public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

  protected:
    void keyPressEvent(QKeyEvent* event);

  private:
    void setupValidators();
    void setupPlot();
    void setupKernelsTable();

    const qreal MAX_X           = 999.0;
    const qreal MIN_X           = -999.0;
    const qreal MAX_Y           = 99.0;
    const qreal MIN_Y           = -99.0;
    const qreal MIN_SMOOTHING_P = 0.0;
    const qreal MAX_SMOOTHING_P = 2.0;
    const int   DECIMAL_NUMBERS = 3;
    const qreal DEFAULT_MIN_X   = -5;
    const qreal DEFAULT_MAX_X   = 5;
    const qreal DEFAULT_MIN_Y   = -0.05;
    const qreal DEFAULT_MAX_Y   = 1;

    std::vector<std::vector<std::shared_ptr<cluster>>> storedMedoids;

    Ui::MainWindow *ui;

    QVector<std::shared_ptr<QVector<qreal>>> samples;
    std::vector<std::shared_ptr<sample>> objects;
    std::vector<std::shared_ptr<cluster>> clusters;

    std::unordered_map<std::string, attributeData*> attributesData;

    std::shared_ptr<dataParser> parser;
    std::shared_ptr<dataReader> reader;

    QStringList kernelTypes;

    int insertObjectsBetweenIntervals(unsigned int objectsNumber);
      int generateInterIntervalObjects(std::vector<std::shared_ptr<sample>> *interIntervalObjects,
                                       unsigned int objectsNumber);
      int selectDesiredNumberOfInterIntervalObjects(std::vector<std::shared_ptr<sample>> *interIntervalObjects);

    void drawPlots(kernelDensityEstimator* estimator, function* targetFunction);
      void clearPlot();
      void resizePlot();
      void addPlot(const QVector<qreal> *X, const QVector<qreal> *Y);
      void addModelPlot(const QVector<qreal> *X, const QVector<qreal> *Y);
      void addEstimatedPlot(const QVector<qreal> *X, const QVector<qreal> *Y);

    void fillStandardDeviations(QVector<std::shared_ptr<QVector<qreal> > > *stDevs);
    void fillMeans(QVector<std::shared_ptr<QVector<qreal> > > *means);

  private slots:

    void refreshKernelsTable();
      void addKernelToTable(int rowNumber, QDoubleValidator *smoothingParameterValidator);
    void refreshTargetFunctionTable();
      void uniformContributions();
        qreal countLastContribution();
    void updateLastContribution();

    void on_pushButton_generate_clicked();
      void fillDomain(QVector<std::shared_ptr<point> > *domain, std::shared_ptr<point> *prototypePoint);
      void generateSamples(QVector<std::shared_ptr<QVector<qreal> > > *means,
                           QVector<std::shared_ptr<QVector<qreal> > > *stDevs);
        distribution *generateTargetDistribution(QVector<std::shared_ptr<QVector<qreal> > > *means,
                                                 QVector<std::shared_ptr<QVector<qreal> > > *stDevs);
        reservoirSamplingAlgorithm *generateReservoirSamplingAlgorithm(dataReader *reader,
                                                                       dataParser *parser);
      kernelDensityEstimator *generateKernelDensityEstimator(int dimensionsNumber);
      function *generateTargetFunction(QVector<std::shared_ptr<QVector<qreal> > > *means, QVector<std::shared_ptr<QVector<qreal> > > *stDevs);
      QColor getRandomColor();
      void testKDE(kernelDensityEstimator* KDE, function* targetFunction);
        int testKDEError(kernelDensityEstimator *KDE, function *targetFunction);
        int testRareElementsDetector(kernelDensityEstimator *KDE);
        void fillTestDomain(QVector<point *> *domain, point* prototypePoint);

    void on_pushButton_animate_clicked();
      int canAnimationBePerformed(int dimensionsNumber);
      void clusterMassiveData(std::vector<std::shared_ptr<sample>> *objects,
                              std::vector<std::vector<std::shared_ptr<cluster>>> *storage);
      void delay(int ms);

    void on_pushButton_clear_clicked();

    void on_spinBox_dimensionsNumber_editingFinished();

    void on_pushButton_countSmoothingParameters_clicked();
      smoothingParameterCounter *generateSmoothingParameterCounter(QVector<qreal> *samplesColumn);

    void on_pushButton_addTargetFunction_clicked();

    void on_pushButton_removeTargetFunction_clicked();

    void updateWeights();

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

enum smoothingParameterCountingMethods
{
  RANK_2_PLUG_IN      = 0,
  RANK_3_PLUG_IN      = 1,
  WEIGHTED_SILVERMAN  = 2
};

enum reservoirSamplingAlgorithms
{
    BASIC_RESERVOIR_SAMPLING_ALGORITHM = 0,
    BIASED_RESERVOIR_SAMPLING_ALGORITHM = 1
};

#endif // MAINWINDOW_H
