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
        QColor getRandomColor();
        qreal countKDEEstimationValue(qreal x);

    void on_pushButton_clear_clicked();

private:
    Ui::MainWindow *ui;

    qreal   mean, standardDeviation;
    int     sampleSize, seed;
};

#endif // MAINWINDOW_H
