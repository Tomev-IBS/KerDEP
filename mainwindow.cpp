#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "math.h"

#include "Distributions/distributions.h"
#include "KDE/kerneldensityestimator.h"
#include "KDE/pluginsmoothingparametercounter.h"

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
    QDoubleValidator* xAxisValidator = new QDoubleValidator(-10.0, 10.0, 3, this);
    xAxisValidator->setLocale(locale);
    xAxisValidator->setNotation(QDoubleValidator::StandardNotation);
    QDoubleValidator* yAxisValidator = new QDoubleValidator(-0.05, 1.0, 3, this);
    yAxisValidator->setLocale(locale);
    yAxisValidator->setNotation(QDoubleValidator::StandardNotation);

    ui->lineEdit_sampleSize->setValidator(seedAndSizeValidator);
    ui->lineEdit_seed->setValidator(seedAndSizeValidator);

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
    ui->tableWidget_dimensionKernels->horizontalHeader()->setStretchLastSection(true);

    refreshKernelsTable();
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

    // Generate a vector of values from normal distribution
    kernel* gaussianProbabilityDensityFunc = new normalKernel();

    // TODO TR: Need a multidimensional vector of attributes
    // To keep things simple let's consider only these domains wherein each dimension has equal size
    QVector<point*> domain;
    int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();
    QVector<qreal> X;
    QVector<qreal> normalDistributionY;

    // Fill domain with points
    fillDomain(&domain, NULL);

    foreach(auto x, domain)
    {
        normalDistributionY.append(gaussianProbabilityDensityFunc->getValue(x));
        X.append(x->at(0));
    }

    // Generate plot of normal distribution using QCustomPlot
    if(dimensionsNumber == 1)
    {
        addPlot(&X, &normalDistributionY);
    }

    // Generate samples
    generateSamples();

    // Generate KDE
    QVector<int> kernelsIDs;
    QVector<qreal> smoothingParameters;
    QVector<QString> carriersRestrictions;

    for(int rowNumber = 0; rowNumber < dimensionsNumber; ++rowNumber)
    {
        kernelsIDs.append(((QComboBox*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, KERNEL_COLUMN_INDEX)))->currentIndex());
        smoothingParameters.append(((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->text().toDouble());
        carriersRestrictions.append(((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX)))->text());
    }

    kernelDensityEstimator* estimator = new kernelDensityEstimator(
                                            &samples,
                                            &smoothingParameters,
                                            &carriersRestrictions,
                                            PRODUCT,
                                            &kernelsIDs
    );

    // Generate a vector of values from selected KDE

    QVector<qreal> KDEEstimationY;

    // TODO: Place counting in another thread

    foreach(QVector<qreal>* x, domain)
    {
        KDEEstimationY.append(estimator->getValue(x));
    }

    // Generate a plot of KDE
    if(dimensionsNumber == 1)
    {
        addPlot(&X, &KDEEstimationY);

        // Draw plots
        ui->widget_plot->replot();
    }
}

void MainWindow::clearPlot()
{
    while(ui->widget_plot->graphCount() != 0)
        ui->widget_plot->removeGraph(0);

    ui->widget_plot->replot();
}

void MainWindow::addPlot(const QVector<qreal> *X, const QVector<qreal> *Y)
{
    ui->widget_plot->addGraph();
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setData(*X, *Y);
    ui->widget_plot->graph(ui->widget_plot->graphCount()-1)->setPen(QPen(getRandomColor()));
}

void MainWindow::fillDomain(QVector<point*>* domain, point *prototypePoint)
{
    // Check if domain is nullpointer
    if(domain == NULL) return;

   // Check if prototype is a null pointer
    if(prototypePoint == NULL)
    {
        // If so make it a point pointer
        prototypePoint = new point();
    }

    qreal val = ui->lineEdit_minX->text().toDouble();

    while(val <= ui->lineEdit_maxX->text().toDouble())
    {
        prototypePoint->append(val);

        if(prototypePoint->size() == ui->spinBox_dimensionsNumber->value())
        {
            domain->append(new point());

            foreach(qreal dimensionVal, *prototypePoint)
            {
                domain->last()->append(dimensionVal);
            }
        }
        else
        {
            fillDomain(domain, prototypePoint);
        }

        prototypePoint->removeLast();

        val += ui->lineEdit_domainDensity->text().toDouble();
    }
}

void MainWindow::generateSamples()
{
    samples.clear();

    qreal seed = ui->lineEdit_seed->text().toDouble();

    // TODO TR: May be selectable in the future.
    distribution* targetDistribution = new normalDistribution(seed);
    int sampleSize = ui->lineEdit_sampleSize->text().toInt();

    if(sampleSize < 1)
    {
        qDebug() << "Sample size < 1.";
        return;
    }

    for(int sampleNumber = 0; sampleNumber < sampleSize; ++sampleNumber)
        samples.append(targetDistribution->getValue());
}

QColor MainWindow::getRandomColor()
{
    return QColor(rand()%110 + 50, rand()%110 + 50, rand()%110 + 50);
}

void MainWindow::on_pushButton_clear_clicked()
{
    clearPlot();
}

void MainWindow::on_spinBox_dimensionsNumber_editingFinished()
{
    refreshKernelsTable();
}

void MainWindow::refreshKernelsTable()
{
    // Get new number of rows
    int newNumberOfRows = ui->spinBox_dimensionsNumber->value();

    // If new number of rows is equal to current number of rows do nothing
    if(newNumberOfRows == ui->tableWidget_dimensionKernels->rowCount())
    {
        return;
    }

    // Set combo box options
    QStringList comboBoxOptions;
    comboBoxOptions << "Normal" << "Triangle" << "Epanecznikow" << "Dull";

    // Set new row count
    ui->tableWidget_dimensionKernels->setRowCount(newNumberOfRows);

    QLocale locale = QLocale::English;
    locale.setNumberOptions(QLocale::c().numberOptions());

    QDoubleValidator* smoothingParameterValidator = new QDoubleValidator(-2.0, 2.0, 3, this);
    smoothingParameterValidator->setLocale(locale);
    smoothingParameterValidator->setNotation(QDoubleValidator::StandardNotation);

    for(int rowNumber = 0; rowNumber < newNumberOfRows; ++rowNumber)
    {
        // Add combobox with kernels

        // TODO TR: Ensure that this doesn't result in memory leaks
        ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, KERNEL_COLUMN_INDEX, new QComboBox());

        ((QComboBox*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, KERNEL_COLUMN_INDEX)))->insertItems(0, comboBoxOptions);

        // Add input box with validator for smoothing parameters
        ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX, new QLineEdit());

        // TODO TR: Ensure that this doesn't result in memory leaks
        ((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->setText("1.0");
        ((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))->setValidator(smoothingParameterValidator);

        // Add input box for carrier restriction value
        ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX, new QLineEdit());

        // TODO TR: Ensure that this doesn't result in memory leaks
        ((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX)))->setText("0.0");
    }
}

void MainWindow::on_pushButton_countSmoothingParameters_clicked()
{
    generateSamples();

    pluginSmoothingParameterCounter counter(&samples);

    qreal value;

    switch(ui->comboBox_smoothingParameterCountingMethod->currentIndex())
    {
        case RANK_3_PLUG_IN:
            value = counter.count3rdRankPluginSmoothingParameter();
        break;
        case RANK_2_PLUG_IN:
        default:
            value = counter.count2ndRankPluginSmoothingParameter();
        break;
    }

    int numberOfRows = ui->tableWidget_dimensionKernels->rowCount();

    for(int rowNumber = 0; rowNumber < numberOfRows; ++rowNumber)
    {
        ((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))
                ->setText(QString::number(value));
    }
}
