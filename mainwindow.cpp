#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "math.h"

#include "Functions/multivariatenormalprobabilitydensityfunction.h"
#include "Functions/complexfunction.h"

#include "Distributions/distributions.h"
#include "KDE/pluginsmoothingparametercounter.h"

#include "Reservoir_sampling/biasedReservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/basicReservoirSamplingAlgorithm.h"

#include "Reservoir_sampling/distributiondataparser.h"
#include "Reservoir_sampling/distributiondatareader.h"
#include "Reservoir_sampling/progressivedistributiondatareader.h"

#include "QDebug"

#include "climits"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    kernelTypes << "Normal" << "Triangle" << "Epanecznikow" << "Dull";

    ui->setupUi(this);

    setupValidators();
    setupPlot();
    setupKernelsTable();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::setupValidators()
{
    QLocale locale = QLocale::English;
    locale.setNumberOptions(QLocale::c().numberOptions());

    const QIntValidator* seedAndSizeValidator = new QIntValidator(1, std::numeric_limits<int>::max(), this);
    QDoubleValidator* xAxisValidator = new QDoubleValidator(MIN_X, MAX_X, DECIMAL_NUMBERS, this);
    xAxisValidator->setLocale(locale);
    xAxisValidator->setNotation(QDoubleValidator::StandardNotation);
    QDoubleValidator* yAxisValidator = new QDoubleValidator(MIN_Y, MAX_Y, DECIMAL_NUMBERS, this);
    yAxisValidator->setLocale(locale);
    yAxisValidator->setNotation(QDoubleValidator::StandardNotation);

    ui->lineEdit_sampleSize->setValidator(seedAndSizeValidator);
    ui->lineEdit_seed->setValidator(seedAndSizeValidator);

    ui->lineEdit_minX->setValidator(xAxisValidator);
    ui->lineEdit_maxX->setValidator(xAxisValidator);

    ui->lineEdit_minY->setValidator(yAxisValidator);
    ui->lineEdit_maxY->setValidator(yAxisValidator);
}

void MainWindow::setupPlot()
{
    ui->widget_plot->xAxis->setLabel("x");
    ui->widget_plot->yAxis->setLabel("f(x)");

    ui->widget_plot->xAxis->setRange(DEFAULT_MIN_X, DEFAULT_MAX_X);
    ui->widget_plot->yAxis->setRange(DEFAULT_MIN_Y, DEFAULT_MAX_Y);
}

void MainWindow::setupKernelsTable()
{
    ui->tableWidget_dimensionKernels->horizontalHeader()->setStretchLastSection(true);

    refreshKernelsTable();
    refreshTargetFunctionTable();
}

void MainWindow::on_pushButton_generate_clicked()
{
    // Log that application started generating KDE
    qDebug() << "KDE generation started.";
    qDebug() << "Seed: " + ui->lineEdit_seed->text() +
                ", Sample size: " + ui->lineEdit_sampleSize->text();

    int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount(),
        targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

    QVector<QVector<qreal>*> means, stDevs;
    QVector<qreal> contributions;
    QVector<function*> elementalFunctions;

    // Generate samples
    generateSamples(&means, &stDevs);

    // Check if contributions are set correctly. If they are, then last contribution is >= 0;
    if(((QLineEdit*)(ui
                     ->tableWidget_targetFunctions
                     ->cellWidget(targetFunctionElementsNumber -1, CONTRIBUTION_COLUMN_INDEX))
                    )
            ->text().toDouble() <= 0)
    {
        // If not then uniform distributions and log error
        qDebug() << "Contributions aren't set correctly. Uniforming contributions.";
        uniformContributions();
    }

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        contributions.append
        (
            ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(functionIndex, CONTRIBUTION_COLUMN_INDEX)))
            ->text().toDouble()
        );

        elementalFunctions.append(new multivariateNormalProbabilityDensityFunction(means.last(), stDevs.last()));
    }

    // Generate a vector of values from normal distribution
    function* targetFunction = new complexFunction(&contributions, &elementalFunctions);

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
    testKDE(estimator, targetFunction);

    // Run plot related tasks if dimension number is equal to 1
    if(dimensionsNumber == 1)
    {
        // Check if prior plots should be saved
        if(!ui->checkBox_keepPriorPlots->isChecked())
        {
            // If not clear plot
            clearPlot();
        }

        resizePlot();

        QVector<point*> domain;

        QVector<qreal> X;
        QVector<qreal> normalDistributionY;

        // Fill domain with points
        // To keep things simple let's consider only these domains wherein each dimension has equal size
        fillDomain(&domain, NULL);

        foreach(auto x, domain)
        {
            normalDistributionY.append(targetFunction->getValue(x));
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

void MainWindow::resizePlot()
{
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

void MainWindow::generateSamples(QVector<QVector<qreal> *> *means, QVector<QVector<qreal> *> *stDevs)
{
    samples.clear();
    objects.clear();

    int sampleSize = ui->lineEdit_sampleSize->text().toInt();

    if(sampleSize < 1)
    {
        qDebug() << "Sample size < 1.";
        return;
    }

    qreal seed = ui->lineEdit_seed->text().toDouble();
    srand(seed);

    // TODO TR: May be selectable in the future.
    int dimensionsNumber = ui->spinBox_dimensionsNumber->value(),
        targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

    QVector<qreal> contributions;
    QVector<distribution*> elementalDistributions;

    for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex)
    {
        means->append(new QVector<qreal>());
        stDevs->append(new QVector<qreal>());

        contributions.append
        (
            ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(functionIndex, CONTRIBUTION_COLUMN_INDEX)))
            ->text().toDouble()
        );

        for(int dimensionIndex = 0; dimensionIndex < dimensionsNumber; ++dimensionIndex)
        {
            means->last()->append
            (
                ((QLineEdit*)(
                    ((QTableWidget*)(ui->tableWidget_targetFunctions->cellWidget(functionIndex, MEAN_COLUMN_INDEX)))
                    ->cellWidget(dimensionIndex, 0)
                ))
                ->text().toDouble()
            );

            stDevs->last()->append
            (
                ((QLineEdit*)(
                    ((QTableWidget*)(ui->tableWidget_targetFunctions->cellWidget(functionIndex, STDEV_COLUMN_INDEX)))
                    ->cellWidget(dimensionIndex, 0)
                ))
                ->text().toDouble()
            );
        }

        elementalDistributions.append(new normalDistribution(seed, means->last(), stDevs->last()));
    }

    distribution* targetDistribution = new complexDistribution(seed, &elementalDistributions, &contributions);

    bool progressiveDistribution = ui->checkBox_dynamicDistribution->isChecked();

    dataParser *parser = new distributionDataParser();
    dataReader *reader;

    if(progressiveDistribution)
    {
        qreal progressionSize = ui->lineEdit_distributionProgression->text().toDouble();
        reader = new progressiveDistributionDataReader(targetDistribution, progressionSize);
    }
    else reader = new distributionDataReader(targetDistribution);

    reservoirSamplingAlgorithm *samplingAlgorithm;

    int stepsNumber = ui->lineEdit_iterationsNumber->text().toDouble();

    int samplingAlgorithmID = ui->comboBox_samplingAlgorithm->currentIndex();

    switch(samplingAlgorithmID)
    {
        case BIASED_RESERVOIR_SAMPLING_ALGORITHM:
            samplingAlgorithm = new biasedReservoirSamplingAlgorithm(reader, parser, sampleSize, stepsNumber);
        break;
        case BASIC_RESERVOIR_SAMPLING_ALGORITHM:
        default:
            samplingAlgorithm = new basicReservoirSamplingAlgorithm(reader, parser, sampleSize, stepsNumber);
        break;
    }

    samplingAlgorithm->fillReservoir(&objects);

    foreach(auto object, objects)
    {
        samples.append(&(static_cast<distributionDataSample*>(object)->values));
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

            foreach(qreal dimensionVal, *prototypePoint) domain->last()->append(dimensionVal);
        }
        else fillTestDomain(domain, prototypePoint);

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
    refreshTargetFunctionTable();
}

void MainWindow::refreshKernelsTable()
{
    // Get new number of rows
    int newNumberOfRows = ui->spinBox_dimensionsNumber->value();

    // If new number of rows is equal to current number of rows do nothing
    if(newNumberOfRows == ui->tableWidget_dimensionKernels->rowCount()) return;

    // Set new row count
    ui->tableWidget_dimensionKernels->setRowCount(newNumberOfRows);

    QLocale locale = QLocale::English;
    locale.setNumberOptions(QLocale::c().numberOptions());

    QDoubleValidator* smoothingParameterValidator = new QDoubleValidator(MIN_SMOOTHING_P, MAX_SMOOTHING_P, DECIMAL_NUMBERS, this);
    smoothingParameterValidator->setLocale(locale);
    smoothingParameterValidator->setNotation(QDoubleValidator::StandardNotation);

    for(int rowNumber = 0; rowNumber < newNumberOfRows; ++rowNumber) addKernelToTable(rowNumber, smoothingParameterValidator);
}

void MainWindow::addKernelToTable(int rowNumber, QDoubleValidator* smoothingParameterValidator)
{
    // Add combobox with kernels

    // TODO TR: Ensure that this doesn't result in memory leaks
    ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, KERNEL_COLUMN_INDEX, new QComboBox());

    ((QComboBox*)(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, KERNEL_COLUMN_INDEX)))->insertItems(0, kernelTypes);

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

void MainWindow::refreshTargetFunctionTable()
{
    int numberOfRows = ui->tableWidget_targetFunctions->rowCount(),
        dimensionsNumber = ui->spinBox_dimensionsNumber->value();

    // Ensure that rows number is at least 1
    if(numberOfRows == 0)
        numberOfRows = 1;

    // Set row count
    ui->tableWidget_targetFunctions->setRowCount(numberOfRows);

    QLocale locale = QLocale::English;
    locale.setNumberOptions(QLocale::c().numberOptions());

    QDoubleValidator* meanValidator = new QDoubleValidator(-10.0, 10.0, 3, this);
    meanValidator->setLocale(locale);
    meanValidator->setNotation(QDoubleValidator::StandardNotation);

    QDoubleValidator* stDevValidator = new QDoubleValidator(-5.0, 5.0, 3, this);
    stDevValidator->setLocale(locale);
    stDevValidator->setNotation(QDoubleValidator::StandardNotation);

    QDoubleValidator* contributionValidator = new QDoubleValidator(0.0, 100.0, 3, this);
    contributionValidator->setLocale(locale);
    contributionValidator->setNotation(QDoubleValidator::StandardNotation);

    QTableWidget *targetFunctionTablePointer = (QTableWidget*)ui->tableWidget_targetFunctions,
                 *meansTablePointer, *stDevsTablePointer;

    for(int rowIndex = 0; rowIndex < numberOfRows; ++rowIndex)
    {
        // TODO TR: Ensure that this doesn't result in memory leaks
        targetFunctionTablePointer->setCellWidget(rowIndex, MEAN_COLUMN_INDEX, new QTableWidget());

        meansTablePointer = ((QTableWidget*)ui->tableWidget_targetFunctions->cellWidget(rowIndex, MEAN_COLUMN_INDEX));
        meansTablePointer->setRowCount(dimensionsNumber);
        meansTablePointer->setColumnCount(1);
        meansTablePointer->horizontalHeader()->hide();

        // TODO TR: Ensure that this doesn't result in memory leaks
        targetFunctionTablePointer->setCellWidget(rowIndex, STDEV_COLUMN_INDEX, new QTableWidget());

        stDevsTablePointer = (QTableWidget*)ui->tableWidget_targetFunctions->cellWidget(rowIndex, STDEV_COLUMN_INDEX);
        stDevsTablePointer->setRowCount(dimensionsNumber);
        stDevsTablePointer->setColumnCount(1);
        stDevsTablePointer->horizontalHeader()->hide();

        for(int dimensionNumber = 0; dimensionNumber < dimensionsNumber; ++dimensionNumber)
        {
           meansTablePointer->setCellWidget(dimensionNumber, 0, new QLineEdit());
           ((QLineEdit*)(meansTablePointer->cellWidget(dimensionNumber, 0)))->setText("0.0");
           ((QLineEdit*)(meansTablePointer->cellWidget(dimensionNumber, 0)))->setValidator(meanValidator);

           stDevsTablePointer->setCellWidget(dimensionNumber, 0, new QLineEdit());
           ((QLineEdit*)(stDevsTablePointer->cellWidget(dimensionNumber, 0)))->setText("1.0");
           ((QLineEdit*)(stDevsTablePointer->cellWidget(dimensionNumber, 0)))->setValidator(stDevValidator);
        }

        // TODO TR: Ensure that this doesn't result in memory leaks
        targetFunctionTablePointer->setCellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX, new QLineEdit());
        ((QLineEdit*)(targetFunctionTablePointer->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->setMaxLength(6);
        ((QLineEdit*)(targetFunctionTablePointer->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->setValidator(contributionValidator);
        QObject::connect(((QLineEdit*)(targetFunctionTablePointer->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX))), SIGNAL(textEdited(QString)), this, SLOT(updateLastContribution()));
    }

    // Disable last contribution cell, as it's filled automatically
    ((QLineEdit*)(targetFunctionTablePointer->cellWidget(numberOfRows -1, CONTRIBUTION_COLUMN_INDEX)))->setEnabled(false);

    uniformContributions();
}

void MainWindow::uniformContributions()
{
    int numberOfRows = ui->tableWidget_targetFunctions->rowCount(), lastRowIndex = numberOfRows - 1;

    for(int rowIndex = 0; rowIndex < lastRowIndex; ++rowIndex)
        ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->setText(QString::number(100.0/numberOfRows));

    ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(lastRowIndex, CONTRIBUTION_COLUMN_INDEX)))->setText(QString::number(countLastContribution()));
}

qreal MainWindow::countLastContribution()
{
    qreal result = 100.0;

    int lastRowIndex = ui->tableWidget_targetFunctions->rowCount()-1;

    for(int rowIndex = 0; rowIndex < lastRowIndex; ++rowIndex)
        result -= ((QLineEdit*)(ui->tableWidget_targetFunctions->cellWidget(rowIndex, CONTRIBUTION_COLUMN_INDEX)))->text().toDouble();

    return result;
}

void MainWindow::updateLastContribution()
{
    int lastRowIndex = ui->tableWidget_targetFunctions->rowCount()-1;
    qreal lastContributionValue = countLastContribution();

    ((QLineEdit*)(ui
                  ->tableWidget_targetFunctions
                  ->cellWidget(lastRowIndex, CONTRIBUTION_COLUMN_INDEX)))
                  ->setText(QString::number(lastContributionValue));
}

void MainWindow::on_pushButton_countSmoothingParameters_clicked()
{
    QVector<QVector<qreal> *> means, stDevs;

    generateSamples(&means, &stDevs);

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

void MainWindow::on_pushButton_addTargetFunction_clicked()
{
    int newRowsNumber = ui->tableWidget_targetFunctions->rowCount() +1;

    // Set row count
    ui->tableWidget_targetFunctions->setRowCount(newRowsNumber);

    refreshTargetFunctionTable();
}

void MainWindow::on_pushButton_removeTargetFunction_clicked()
{
    int newRowsNumber = ui->tableWidget_targetFunctions->rowCount() -1;

    // Set row count
    ui->tableWidget_targetFunctions->setRowCount(newRowsNumber);

    refreshTargetFunctionTable();
}
