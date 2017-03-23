#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_pushButton_generate_clicked();
        void clearPlot();
        void generateSamples();
        QColor getRandomColor();

    void on_pushButton_clear_clicked();

    void on_spinBox_dimensionsNumber_editingFinished();
        void refreshKernelsTable();

        void on_pushButton_countSmoothingParameters_clicked();

private:
    Ui::MainWindow *ui;

    QVector<qreal> samples;
};

enum kernelSettingsColumns
{
    KERNEL_COLUMN_INDEX                 = 0,
    SMOOTHING_PARAMETER_COLUMN_INDEX    = 1,
    CARRIER_RESTRICTION_COLUMN_INDEX    = 2
};

enum smoothingParameterCountingMethods
{
    RANK_2_PLUG_IN    = 0,
    RANK_3_PLUG_IN    = 1
};

#endif // MAINWINDOW_H
