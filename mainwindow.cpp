#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "math.h"

#include "random"
#include "QDebug"

#include "climits"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QLocale locale = QLocale::English;
    locale.setNumberOptions(QLocale::c().numberOptions());

    // Set validators for lineEdits
    const QIntValidator* seedAndSizeValidator = new QIntValidator(1, std::numeric_limits<int>::max(), this);
    QDoubleValidator* meanAndDeviationValidator = new QDoubleValidator(-5.0, 5.0, 3, this);
    meanAndDeviationValidator->setLocale(locale);
    meanAndDeviationValidator->setNotation(QDoubleValidator::StandardNotation);
    QDoubleValidator* xAxisValidator = new QDoubleValidator(-10.0, 10.0, 3, this);
    xAxisValidator->setLocale(locale);
    xAxisValidator->setNotation(QDoubleValidator::StandardNotation);
    QDoubleValidator* yAxisValidator = new QDoubleValidator(-0.05, 1.0, 3, this);
    yAxisValidator->setLocale(locale);
    yAxisValidator->setNotation(QDoubleValidator::StandardNotation);

    ui->lineEdit_sampleSize->setValidator(seedAndSizeValidator);
    ui->lineEdit_seed->setValidator(seedAndSizeValidator);

    ui->lineEdit_mean->setValidator(meanAndDeviationValidator);
    ui->lineEdit_stdDeviation->setValidator(meanAndDeviationValidator);

    ui->lineEdit_minX->setValidator(xAxisValidator);
    ui->lineEdit_maxX->setValidator(xAxisValidator);

    ui->lineEdit_minY->setValidator(yAxisValidator);
    ui->lineEdit_maxY->setValidator(yAxisValidator);

    // Setup plot

    ui->widget_plot->xAxis->setLabel("x");
    ui->widget_plot->yAxis->setLabel("f(x)");

    ui->widget_plot->xAxis->setRange(-5, 5);
    ui->widget_plot->yAxis->setRange(-0.05, 0.8);

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

    // Resize plot
    qreal   minX = ui->lineEdit_minX->text().toDouble(),
            maxX = ui->lineEdit_maxX->text().toDouble(),
            minY = ui->lineEdit_minY->text().toDouble(),
            maxY = ui->lineEdit_maxY->text().toDouble();

    // Check if sizes are entered correctly
    if(minX < maxX)
    {
        // If so change
        ui->widget_plot->xAxis->setRange(minX, maxX);
    }
    else
    {
        // If not log it and do nothing
        qDebug() << "Minimal x value cannot be lower than it's maximal value.";
    }

    if(minY < maxY)
    {
        // If so change
        ui->widget_plot->yAxis->setRange(minY, maxY);
    }
    else
    {
        // If not log it and do nothing
        qDebug() << "Minimal y value cannot be lower than it's maximal value.";
    }


    // Create random engine generator
    seed                = ui->lineEdit_seed->text().toInt();
    mean                = ui->lineEdit_mean->text().toDouble();
    standardDeviation   = ui->lineEdit_stdDeviation->text().toDouble();
    std::default_random_engine generator(seed);

    std::normal_distribution<qreal> distribution(mean, standardDeviation);

    // Generate a vector of values from normal distribution
    int sampleSize = ui->lineEdit_sampleSize->text().toInt();

    QList<qreal> samples;

    QVector<double> X;
    QVector<double> normalDistributionY;

    for(int x = -500; x < 500; ++x)
    {
        X.append(x/100.0);
        normalDistributionY.append(countNormalDistributionDensityValue(x/100.0));
    }

    // Generate plot of normal distribution using QCustomPlot
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(0)->setData(X, normalDistributionY);
    ui->widget_plot->replot();

    for(int sampleNumber = 0; sampleNumber < sampleSize; ++sampleNumber)
    {
        qreal number = distribution(generator);
        samples.append(number);
    }
}

qreal MainWindow::countNormalDistributionDensityValue(qreal x)
{
    qreal result = exp(- pow((x - mean), 2) / (2 * pow(standardDeviation, 2)));
    result /= (standardDeviation * sqrt(2 * M_PI));

    return result;
}
