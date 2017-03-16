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

    void on_pushButton_clear_clicked();

    void on_spinBox_dimensionsNumber_editingFinished();
        void refreshKernelsTable();

private:
    Ui::MainWindow *ui;
};

enum kernelSettingsColumns
{
    KERNEL_COLUMN_INDEX                = 0,
    SMOOTHING_PARAMETER_COLUMN_INDEX   = 1
};

#endif // MAINWINDOW_H
