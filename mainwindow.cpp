#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "math.h"

#include "Functions/multivariatenormalprobabilitydensityfunction.h"
#include "Distributions/distributions.h"
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
    refreshTargetFunctionTable();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_generate_clicked()
{
    // Log that application started generating KDE
    qDebug() << "KDE generation started.";
    qDebug() << "Seed: " + ui->lineEdit_seed->text() +
                ", Sample size: " + ui->lineEdit_sampleSize->text();

    int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

    QVector<qreal> means, stDevs;

    for(int i = 0; i < dimensionsNumber; ++i)
    {
        means.append(0);
        stDevs.append(1);
    }

    // Generate a vector of values from normal distribution
    function* normalDistributionProbabilityDensityFunction = new multivariateNormalProbabilityDensityFunction(&means, &stDevs);

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

    // Test estimator
    testKDE(estimator, normalDistributionProbabilityDensityFunction);

    // Run plot related tasks if dimension number is equal to 1
    if(dimensionsNumber == 1)
    {
        // Check if prior plots should be saved
        if(!ui->checkBox_keepPriorPlots->isChecked())
        {
            // If not clear plot
            clearPlot();
        }

        // Resize plot
        // TR TODO: Place it in another method
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

        QVector<point*> domain;

        QVector<qreal> X;
        QVector<qreal> normalDistributionY;

        // Fill domain with points
        // To keep things simple let's consider only these domains wherein each dimension has equal size
        fillDomain(&domain, NULL);

        foreach(auto x, domain)
        {
            normalDistributionY.append(normalDistributionProbabilityDensityFunction->getValue(x));
            X.append(x->at(0));
        }

        // Generate plot of normal distribution using QCustomPlot
        addPlot(&X, &normalDistributionY);

        // Generate a vector of values from selected KDE
        QVector<qreal> KDEEstimationY;

        // TODO: Place counting in another thread
        foreach(QVector<qreal>* x, domain)
        {
            KDEEstimationY.append(estimator->getValue(x));
        }

        // Generate a plot of KDE
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
    int dimensionsNum = ui->spinBox_dimensionsNumber->value();
    QVector<qreal> means;
    QVector<qreal> stDevs;

    while(means.size() != dimensionsNum)
    {
        means.append(0);
        stDevs.append(1);
    }

    distribution* targetDistribution = new normalDistribution(seed, &means, &stDevs);

    int sampleSize = ui->lineEdit_sampleSize->text().toInt();

    if(sampleSize < 1)
    {
        qDebug() << "Sample size < 1.";
        return;
    }

    for(int sampleNumber = 0; sampleNumber < sampleSize; ++sampleNumber)
    {
        samples.append(new QVector<qreal>());
        targetDistribution->getValue(samples.last());
    }
}

QColor MainWindow::getRandomColor()
{
    return QColor(rand()%110 + 50, rand()%110 + 50, rand()%110 + 50);
}

void MainWindow::testKDE(kernelDensityEstimator *KDE, function *targetFunction)
{
    QVector<point *> testDomain;

    fillTestDomain(&testDomain, NULL);

    qreal error = 0;

    foreach(auto arg, testDomain)
    {
        qDebug()    << "Point: " << *arg
                    << "Target: " << targetFunction->getValue(arg)
                    << "Estimated: " << KDE->getValue(arg)
                    << "Difference: " << qAbs(targetFunction->getValue(arg) - KDE->getValue(arg));
        error += qAbs(targetFunction->getValue(arg) - KDE->getValue(arg));
    }

    qDebug() << "Error: " << error;
    qDebug() << "Average error: " << error/testDomain.size();
}

void MainWindow::fillTestDomain(QVector<point *> *domain, point *prototypePoint)
{
    // Check if domain is nullpointer
    if(domain == NULL) return;

   // Check if prototype is a null pointer
    if(prototypePoint == NULL)
    {
        // If so make it a point pointer
        prototypePoint = new point();
    }

    for(int i = -1; i <= 1; ++i)
    {
        prototypePoint->append(i);

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
            fillTestDomain(domain, prototypePoint);
        }

        prototypePoint->removeLast();
    }
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
        ((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, CARRIER_RESTRICTION_COLUMN_INDEX)))->setText("None.");
    }
}

void MainWindow::refreshTargetFunctionTable()
{
    int numberOfRows = ui->tableWidget_targetFunctions->rowCount();

    // Ensure that rows number is at least 1
    if(numberOfRows == 0)
        numberOfRows = 1;

    // Set row count
    ui->tableWidget_targetFunctions->setRowCount(numberOfRows);

    QLocale locale = QLocale::English;
    locale.setNumberOptions(QLocale::c().numberOptions());

    QIntValidator* meanValidator = new QIntValidator(-10, 10, this);

    QDoubleValidator* stDevValidator = new QDoubleValidator(-5.0, 5.0, 3, this);
    stDevValidator->setLocale(locale);
    stDevValidator->setNotation(QDoubleValidator::StandardNotation);

    QDoubleValidator* contributionValidator = new QDoubleValidator(0.0, 1.0, 3, this);
    contributionValidator->setLocale(locale);
    contributionValidator->setNotation(QDoubleValidator::StandardNotation);

    for(int rowIndex = 0; rowIndex < numberOfRows; ++rowIndex)
    {
        // TODO TR: Ensure that this doesn't result in memory leaks
        ui->tableWidget_targetFunctions->setCellWidget(rowIndex, MEAN_COLUMN_INDEX, new QLineEdit());
        ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(rowIndex, MEAN_COLUMN_INDEX)))->setText("0");
        ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(rowIndex, MEAN_COLUMN_INDEX)))->setValidator(meanValidator);

        // TODO TR: Ensure that this doesn't result in memory leaks
        ui->tableWidget_targetFunctions->setCellWidget(rowIndex, STDEV_COLUMN_INDEX, new QLineEdit());
        ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(rowIndex, STDEV_COLUMN_INDEX)))->setText("1.0");
        ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(rowIndex, STDEV_COLUMN_INDEX)))->setValidator(stDevValidator);

        // TODO TR: Ensure that this doesn't result in memory leaks
        ui->tableWidget_targetFunctions->setCellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX, new QLineEdit());
        ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->setValidator(contributionValidator);
    }

    uniformContributions();
}

void MainWindow::uniformContributions()
{

}

void MainWindow::on_pushButton_countSmoothingParameters_clicked()
{
    generateSamples();

    // Check which method was selected

    qreal(pluginSmoothingParameterCounter::*(methodFunctionPointer))(void);

    switch(ui->comboBox_smoothingParameterCountingMethod->currentIndex())
    {
        case RANK_3_PLUG_IN:
            methodFunctionPointer = &pluginSmoothingParameterCounter::count3rdRankPluginSmoothingParameter;
        break;
        case RANK_2_PLUG_IN:
        default:
            methodFunctionPointer = &pluginSmoothingParameterCounter::count2ndRankPluginSmoothingParameter;
        break;
    }

    // Count smoothing parameter for each dimension

    int numberOfRows = ui->tableWidget_dimensionKernels->rowCount();
    QVector<qreal> samplesColumn;
    pluginSmoothingParameterCounter counter(&samplesColumn);
    qreal value;

    for(int rowNumber = 0; rowNumber < numberOfRows; ++rowNumber)
    {
        // Create vector that consists of variables inside this dimension

        samplesColumn.clear();

        foreach(QVector<qreal>* sample, samples)
            samplesColumn.append(sample->at(rowNumber));

        value = (counter.*methodFunctionPointer)();

        ((QLineEdit*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, SMOOTHING_PARAMETER_COLUMN_INDEX)))
                ->setText(QString::number(value));
    }
}

