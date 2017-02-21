#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "random"
#include "QDebug"

#include "climits"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Set validators for lineEdits to max int
    QIntValidator* lineEditsValidator = new QIntValidator(1, std::numeric_limits<int>::max(), this);

    ui->lineEdit_sampleSize->setValidator(lineEditsValidator);
    ui->lineEdit_seed->setValidator(lineEditsValidator);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_generate_clicked()
{
    // Log that application started generating plot
    qDebug() << "Plot generation started.";
    qDebug() << "Seed: " + ui->lineEdit_seed->text() +
                ", Sample size: " + ui->lineEdit_sampleSize->text();


    // Create random engine generator
    int seed = ui->lineEdit_seed->text().toInt();
    std::default_random_engine generator(seed);

    std::normal_distribution<qreal> distribution(0.5, 0.17);

    // Generate a vector of values from normal distribution
    int sampleSize = ui->lineEdit_sampleSize->text().toInt();

    QList<qreal> samples;

    for(int sampleNumber = 0; sampleNumber < sampleSize; ++sampleNumber)
    {
        qreal number = distribution(generator);
        samples.append(number);
    }

    // Generate plot using QCustomPlot



    qDebug() << samples;
}
