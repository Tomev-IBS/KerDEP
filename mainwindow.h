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
        qreal countNormalDistributionDensityValue(qreal x);
        qreal countKDEEstimationValue(qreal x);

private:
    Ui::MainWindow *ui;

    qreal   mean, standardDeviation;
    int     sampleSize, seed;
};

#endif // MAINWINDOW_H
