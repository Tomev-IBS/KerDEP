#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleValidator>
#include <QStringList>
#include <vector>

#include "Reservoir_sampling/reservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/sample.h"
#include "Functions/Kernels/kernels.h"
#include "KDE/kerneldensityestimator.h"

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
        const int DECIMAL_NUMBERS   = 3;
        const qreal DEFAULT_MIN_X   =-5;
        const qreal DEFAULT_MAX_X   =5;
        const qreal DEFAULT_MIN_Y   =-0.05;
        const qreal DEFAULT_MAX_Y   =1;

        Ui::MainWindow *ui;

        QVector<QVector<qreal>*> samples;
        std::vector<sample*> objects;

        QStringList kernelTypes;

        void drawPlots(kernelDensityEstimator* estimator, function* targetFunction);
            void clearPlot();
            void resizePlot();
            void addPlot(const QVector<qreal> *X, const QVector<qreal> *Y);

        void fillStandardDeviations(QVector<QVector<qreal> *> *stDevs);
        void fillMeans(QVector<QVector<qreal> *> *means);


private slots:

        void refreshKernelsTable();
            void addKernelToTable(int rowNumber, QDoubleValidator *smoothingParameterValidator);
        void refreshTargetFunctionTable();
            void uniformContributions();
                qreal countLastContribution();
        void updateLastContribution();

        void on_pushButton_generate_clicked();
            void fillDomain(QVector<point *> *domain, point* prototypePoint);
            void generateSamples(QVector<QVector<qreal> *> *means, QVector<QVector<qreal> *> *stDevs);
                distribution *generateTargetDistribution(QVector<QVector<qreal> *> *means, QVector<QVector<qreal> *> *stDevs);
                reservoirSamplingAlgorithm *generateReservoirSamplingAlgorithm(dataReader *reader, dataParser *parser);
            kernelDensityEstimator *generateKernelDensityEstimator(int dimensionsNumber);
            function *generateTargetFunction(QVector<QVector<qreal> *> *means, QVector<QVector<qreal> *> *stDevs);
            QColor getRandomColor();
            void testKDE(kernelDensityEstimator* KDE, function* targetFunction);
                void fillTestDomain(QVector<point *> *domain, point* prototypePoint);

        void on_pushButton_animate_clicked();

        void on_pushButton_clear_clicked();

        void on_spinBox_dimensionsNumber_editingFinished();

        void on_pushButton_countSmoothingParameters_clicked();

        void on_pushButton_addTargetFunction_clicked();

        void on_pushButton_removeTargetFunction_clicked();

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
    RANK_2_PLUG_IN    = 0,
    RANK_3_PLUG_IN    = 1
};

enum reservoirSamplingAlgorithms
{
    BASIC_RESERVOIR_SAMPLING_ALGORITHM = 0,
    BIASED_RESERVOIR_SAMPLING_ALGORITHM = 1
};

#endif // MAINWINDOW_H
