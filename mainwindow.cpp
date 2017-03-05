#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "math.h"

#include "Functions/Kernels/kernels.h"
#include "Distributions/distributions.h"
#include "KDE/kerneldensityestimator.h"

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

    // Setup kernels table

    ui->tableWidget_dimensionKernels->setRowCount(1);
    ui->tableWidget_dimensionKernels->setColumnCount(2);

    QStringList dimensionKernelsTableHeaders;
    dimensionKernelsTableHeaders << "#" << "Kernel type";

    ui->tableWidget_dimensionKernels->setHorizontalHeaderLabels(dimensionKernelsTableHeaders);
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
    mean                = ui->lineEdit_mean->text().toDouble();
    standardDeviation   = ui->lineEdit_stdDeviation->text().toDouble();
    seed                = ui->lineEdit_seed->text().toDouble();

    // Generate a vector of values from normal distribution
    kernel* gaussianProbabilityDensityFunc = new normalKernel(mean, standardDeviation);

    QVector<qreal> X;
    QVector<qreal> normalDistributionY;

    QVector<qreal>* tempValueHolder = new QVector<qreal>();

    for(int x = minX*100; x < maxX*100; ++x)
    {
        X.append(x/100.0);
        tempValueHolder->clear();
        tempValueHolder->append(x/100.0);
        normalDistributionY.append(gaussianProbabilityDensityFunc->getValue(tempValueHolder));
    }

    // Generate plot of normal distribution using QCustomPlot
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(X, normalDistributionY);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(getRandomColor()));

    // Generate a vector of values from selected KDE
    int kernelFunctionID = ui->comboBox_kernel->currentIndex();
    kernelDensityEstimator* estimator;
    function* kernel;
    distribution* targetDistribution = new normalDistribution(seed, mean, standardDeviation);

    switch(kernelFunctionID)
    {
        case NORMAL:
            kernel = new normalKernel(mean, standardDeviation);
        break;

        case TRIANGLE:
            kernel = new triangleKernel();
        break;

        case EPANECZNIKOW:
            kernel = new epanecznikowKernel();
        break;
        case DULL:
        default:
            kernel = new dullKernel();
        break;
    }

    estimator = new kernelDensityEstimator(
                ui->lineEdit_sampleSize->text().toInt(),
                ui->lineEdit_smoothingParam->text().toDouble(),
                kernel,
                targetDistribution
    );

    QVector<qreal> KDEEstimationY;

    // TODO: Place counting in another thread

    foreach(qreal x, X) KDEEstimationY.append(estimator->getValue(x));

    // Generate a plot of KDE
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(X, KDEEstimationY);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(getRandomColor()));

    // Draw plots
    ui->widget_plot->replot();
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

void MainWindow::on_pushButton_clear_clicked()
{
    clearPlot();
}
