#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "math.h"
#include "Functions/gaussianprobabilitydensityfunction.h"

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

    // Add

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

    // Check if prior plots should be saved
    if(!ui->checkBox_keepPriorPlots->isChecked())
    {
        // If not clear plot
        clearPlot();
    }

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
        // If not log it and correct
        qDebug() << "Minimal x value cannot be lower than it's maximal value.";
        minX = ui->widget_plot->xAxis->range().minRange;
        maxX = ui->widget_plot->xAxis->range().maxRange;
    }

    if(minY < maxY)
    {
        // If so change
        ui->widget_plot->yAxis->setRange(minY, maxY);
    }
    else
    {
        // If not log it and correct
        qDebug() << "Minimal y value cannot be lower than it's maximal value.";
        minY = ui->widget_plot->yAxis->range().minRange;
        maxY = ui->widget_plot->yAxis->range().maxRange;
    }


    // Create random engine generator
    seed                = ui->lineEdit_seed->text().toInt();
    mean                = ui->lineEdit_mean->text().toDouble();
    standardDeviation   = ui->lineEdit_stdDeviation->text().toDouble();
    std::default_random_engine generator(seed);

    std::normal_distribution<qreal> distribution(mean, standardDeviation);

    // Generate a vector of values from normal distribution
    function* gaussianProbabilityDensityFunc = new gaussianProbabilityDensityFunction(mean, standardDeviation);

    QVector<qreal> X;
    QVector<qreal> normalDistributionY;

    for(int x = minX*100; x < maxX*100; ++x)
    {
        X.append(x/100.0);
        normalDistributionY.append(gaussianProbabilityDensityFunc->getValue(new QVector<qreal>{x/100.0}));
    }

    // Generate plot of normal distribution using QCustomPlot
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(X, normalDistributionY);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(getRandomColor()));

    // Generate a vector of values from selected KDE
    int sampleSize = ui->lineEdit_sampleSize->text().toInt();

    QVector<qreal> KDEEstimationY;

    foreach(qreal x, X) KDEEstimationY.append(countKDEEstimationValue(x));

    // Generate a plot of KDE
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(X, KDEEstimationY);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(getRandomColor()));

    // Draw plots
    ui->widget_plot->replot();

    QList<qreal> samples;

    for(int sampleNumber = 0; sampleNumber < sampleSize; ++sampleNumber)
    {
        qreal number = distribution(generator);
        samples.append(number);
    }
}

void MainWindow::clearPlot()
{
    while(ui->widget_plot->graphCount() != 0)
        ui->widget_plot->removeGraph(0);

    ui->widget_plot->replot();
}

QColor MainWindow::getRandomColor()
{
    return QColor(rand()%110 + 50, rand()%110 + 50, rand()%110 + 50);
}

qreal MainWindow::countKDEEstimationValue(qreal x)
{
    qreal result = 0;

    return result;
}

void MainWindow::on_pushButton_clear_clicked()
{
    clearPlot();
}
