#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "Functions/Kernels/kernels.h"
#include "KDE/kerneldensityestimator.h"

namespace Ui {
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
    const int DECIMAL_NUMBERS   = 3;
    const qreal DEFAULT_MIN_X   =-5;
    const qreal DEFAULT_MAX_X   =5;
    const qreal DEFAULT_MIN_Y   =-0.05;
    const qreal DEFAULT_MAX_Y   =1;

private slots:

    void refreshKernelsTable();
    void refreshTargetFunctionTable();
        void uniformContributions();
            qreal countLastContribution();
    void updateLastContribution();

    void on_pushButton_generate_clicked();
        void clearPlot();
        void addPlot(const QVector<qreal> *X, const QVector<qreal> *Y);
        void fillDomain(QVector<point *> *domain, point* prototypePoint);
        void generateSamples(QVector<QVector<qreal> *> *means, QVector<QVector<qreal> *> *stDevs);
        QColor getRandomColor();
        void testKDE(kernelDensityEstimator* KDE, function* targetFunction);
            void fillTestDomain(QVector<point *> *domain, point* prototypePoint);

    void on_pushButton_clear_clicked();

    void on_spinBox_dimensionsNumber_editingFinished();


        void on_pushButton_countSmoothingParameters_clicked();

        void on_pushButton_addTargetFunction_clicked();

        void on_pushButton_removeTargetFunction_clicked();

private:
    Ui::MainWindow *ui;

    QVector<QVector<qreal>*> samples;
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
