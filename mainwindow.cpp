#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <cfloat>
#include <QDebug>
#include <algorithm>
#include <chrono>
#include <QDateTime>
#include <Benchmarking/errorsCalculator.h>

#include "UI/plotLabel.h"
#include "UI/plotLabelIntDataPreparator.h"

#include "Functions/complexfunction.h"

#include "Reservoir_sampling/biasedReservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/basicReservoirSamplingAlgorithm.h"

#include "Reservoir_sampling/distributiondataparser.h"
#include "Reservoir_sampling/progressivedistributiondatareader.h"

#include "ClusterKernelWrappers/varianceBasedClusterKernel.h"
#include "ClusterKernelWrappers/univariateStreamElement.h"

#include "UI/QwtContourPlotUI.h"

ClusterKernel *CreateNewVarianceBasedClusterKernel(ClusterKernelStreamElement *stream_element){
  auto newClusterKernel = new VarianceBasedClusterKernel(stream_element);
  return newClusterKernel;
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
  kernel_types_ << "Normal" << "Triangle" << "Epanecznikow" << "Dull";

  ui->setupUi(this);

  // Adding contour plot
  contour_plot_ = new Plot(ui->widget_contour_plot);
  auto l = new QGridLayout(ui->widget_contour_plot);
  l->setAlignment(ui->widget_contour_plot, Qt::AlignCenter);
  l->addWidget(contour_plot_, 0, 0, 1, 1, Qt::AlignTop | Qt::AlignLeft);
  QSizePolicy p(QSizePolicy::Maximum, QSizePolicy::Maximum);
  contour_plot_->setSizePolicy(p);

  // Setting up the validators
  SetupValidators();
  SetupPlot();
  SetupKernelsTable();

  testNewFunctionalities();

  // Set font
  auto appFont = QFont("Courier New");
  appFont.setStyleHint(QFont::TypeWriter);
  setFont(appFont);
}

MainWindow::~MainWindow() {
  delete contour_plot_;
  delete ui;
}

void log(const QString &msg){
 qDebug() << QDateTime::currentDateTime().toString("h:m:s ap") << ": " << msg;
}

void MainWindow::testNewFunctionalities() {
  log("Start test.");
  log("Nothing to test.");
  log("Finish test.");
}

void MainWindow::SetupValidators() {
  QLocale locale = QLocale::English;
  locale.setNumberOptions(QLocale::c().numberOptions());

  std::unique_ptr<QIntValidator> naturalNumbersValidator(new QIntValidator(0, std::numeric_limits<int>::max(), this));
  std::unique_ptr<QIntValidator>
      positiveNaturalNumbersValidator(new QIntValidator(1, std::numeric_limits<int>::max(), this));
  std::unique_ptr<QDoubleValidator> xAxisValidator(new QDoubleValidator(kMinX, kMaxX, kDecimalNumbers, this));
  xAxisValidator->setLocale(locale);
  xAxisValidator->setNotation(QDoubleValidator::StandardNotation);

  std::unique_ptr<QDoubleValidator> yAxisValidator(new QDoubleValidator(kMinY, kMaxY, kDecimalNumbers, this));
  yAxisValidator->setLocale(locale);
  yAxisValidator->setNotation(QDoubleValidator::StandardNotation);

  ui->lineEdit_sampleSize->setValidator(positiveNaturalNumbersValidator.get());
  ui->lineEdit_seed->setValidator(naturalNumbersValidator.get());
  ui->lineEdit_milisecondsDelay->setValidator(naturalNumbersValidator.get());

  ui->lineEdit_minX->setValidator(xAxisValidator.get());
  ui->lineEdit_maxX->setValidator(xAxisValidator.get());

  ui->lineEdit_minY->setValidator(yAxisValidator.get());
  ui->lineEdit_maxY->setValidator(yAxisValidator.get());
}

void MainWindow::SetupPlot() {
  ui->widget_plot->xAxis->setRange(kDefaultMinX, kDefaultMaxX);
  ui->widget_plot->yAxis->setRange(kDefaultMinY, kDefaultMaxY);
}

void MainWindow::SetupKernelsTable() {
  ui->tableWidget_dimensionKernels
    ->horizontalHeader()->setStretchLastSection(true);

  RefreshKernelsTable();
  RefreshTargetFunctionTable();
}

void MainWindow::DrawPlots(DESDA *DESDAAlgorithm) {
  ClearPlot();
  ResizePlot();

  std::vector<std::vector<double>> drawable_domain = {}; // This is required for types :P
  for(auto value : drawable_domain_){
    drawable_domain.push_back({value});
  }

  // Generate plot of model function
  if(ui->checkBox_showEstimatedPlot->isChecked()) {
    QVector<qreal> modelDistributionY =
      QVector<qreal>::fromStdVector(
  GetFunctionsValueOnDomain(target_function_.get(), drawable_domain)
      );
    AddPlot(&modelDistributionY, model_plot_pen_);
  }

  // Generate less elements KDE plot (navy blue)
  if(ui->checkBox_showEstimationPlot->isChecked()) {
    auto less_elements_estimator_y = QVector<double>::fromStdVector(
        DESDAAlgorithm->getKDEValues(&drawable_domain));
    AddPlot(&less_elements_estimator_y, kde_plot_pen_);
  }

  // Generate weighted estimator plot (light blue)
  if(ui->checkBox_showWeightedEstimationPlot->isChecked()) {
    auto weighted_estimator_y = QVector<double>::fromStdVector(
        DESDAAlgorithm->getWeightedKDEValues(&drawable_domain));
    AddPlot(&weighted_estimator_y, weighted_plot_pen_);
  }

  // Generate full estimator plot (BLACK)
  if(ui->checbox_showFullEstimator->isChecked()) {
    auto windowed_estimator_y = QVector<double>::fromStdVector(
        DESDAAlgorithm->getWindowKDEValues(&drawable_domain));
    AddPlot(&windowed_estimator_y, windowed_plot_pen_);
  }

  // Generate plot for kernel prognosis derivative
  if(ui->checkBox_kernelPrognosedPlot->isChecked())
    AddPlot(&kernel_prognosis_derivative_values_, derivative_plot_pen_);

  // Generate plot for standardized prognosis derivative, assuming that
  // normal derivative was generated first
  if(ui->checkBox_standarizedDerivative->isChecked()) {
    QVector<double> standardizedDerivativeY = {};
    for(auto val : kernel_prognosis_derivative_values_) {
      standardizedDerivativeY.push_back(
          0.1 * val / DESDAAlgorithm->_maxAbsDerivativeValueInCurrentStep
                                       );
    }
    AddPlot(&standardizedDerivativeY, standardized_derivative_plot_pen_);
  }

  if(ui->checkBox_sigmoidallyEnhancedKDE->isChecked()) {
    auto sigmoidally_enhanced_plot_y = QVector<double>::fromStdVector(
        DESDAAlgorithm->getEnhancedKDEValues(&drawable_domain));
    AddPlot(&sigmoidally_enhanced_plot_y, desda_kde_plot_pen_);
  }

  if(ui->checkBox_showUnusualClusters->isChecked()) {
    atypical_elements_values_and_derivatives_ =
        DESDAAlgorithm->getAtypicalElementsValuesAndDerivatives();
    quantile_estimator_value_ = DESDAAlgorithm->_quantileEstimator;
    MarkUncommonClusters();
  }

  if(ui->checkBox_REESEKDE->isChecked()) {
    auto rare_elements_enhanced_plot_y = QVector<double>::fromStdVector(
        DESDAAlgorithm->getRareElementsEnhancedKDEValues(&drawable_domain));
    AddPlot(&rare_elements_enhanced_plot_y, desda_rare_elements_kde_plot_pen_);
  }
  // Draw plots
  ui->widget_plot->replot();
}

void MainWindow::DrawPlots(EnhancedClusterKernelAlgorithm *CKAlgorithm) {
  ClearPlot();
  ResizePlot();

  std::vector<std::vector<double>> drawable_domain = {}; // This is required for types :P
  for(auto value : drawable_domain_){
    drawable_domain.push_back({value});
  }

  // Generate plot of model function
  if(ui->checkBox_showEstimatedPlot->isChecked()) {
    QVector<qreal> modelDistributionY =
        QVector<qreal>::fromStdVector(
            GetFunctionsValueOnDomain(target_function_.get(), drawable_domain)
                                     );
    AddPlot(&modelDistributionY, model_plot_pen_);
  }

  // Generate less elements KDE plot (navy blue)
  if(ui->checkBox_showEstimationPlot->isChecked()) {
    auto less_elements_estimator_y = QVector<double>::fromStdVector(
        CKAlgorithm->GetKDEValuesOnDomain(drawable_domain));
    AddPlot(&less_elements_estimator_y, kde_plot_pen_);
  }
}

void MainWindow::AddPlot(const QVector<qreal> *Y, const QPen &pen) {
  ui->widget_plot->addGraph();
  ui->widget_plot->graph(ui->widget_plot->graphCount() - 1)->setData(drawable_domain_, *Y);
  ui->widget_plot->graph(ui->widget_plot->graphCount() - 1)->setPen(pen);
}

void MainWindow::ResizePlot() {
  // Resize plot
  qreal minX = ui->lineEdit_minX->text().toDouble(),
        maxX = ui->lineEdit_maxX->text().toDouble(), // Standard
  //maxX = 3 + ui->lineEdit_distributionProgression->text().toDouble() * 3000, // For progression
  //maxX = 3 + ui->lineEdit_distributionProgression->text().toDouble(), // /* should jump */
  minY = ui->lineEdit_minY->text().toDouble(),
      maxY = ui->lineEdit_maxY->text().toDouble();

  // Check if sizes are entered correctly
  if(minX < maxX) {
    // If so change
    ui->widget_plot->xAxis->setRange(minX, maxX);
  }
  else {
    // If not log it and correct
    log("Minimal x value cannot be lower than it's maximal value.");
    minX = kDefaultMinX;
    maxX = kDefaultMaxX;
  }

  if(minY < maxY) {
    // If so change
    ui->widget_plot->yAxis->setRange(minY, maxY);
  }
  else {
    // If not log it and correct
    log("Minimal y value cannot be lower than it's maximal value.");
    minY = kDefaultMinY;
    maxY = kDefaultMaxY;
  }

  ui->widget_plot->xAxis->setTickLabelFont(QFont(font().family(), 14));
  ui->widget_plot->yAxis->setTickLabelFont(QFont(font().family(), 14));

  QVector<double> ticks;
  QVector<QString> labels;
  // Switch to an array

  double i = minX;

  while(i < maxX + 1){
    ticks << i;
    labels << QString::number(i);
    i += 1;
  }

  QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
  textTicker->addTicks(ticks, labels);

  ui->widget_plot->xAxis->setTicker(textTicker);
}

void MainWindow::ClearPlot() {
  while(ui->widget_plot->graphCount() != 0)
    ui->widget_plot->removeGraph(0);

  for(auto a : lines_on_plot_) {
    ui->widget_plot->removeItem(a);
  }

  lines_on_plot_.clear();
}

unsigned long long MainWindow::MarkUncommonClusters() {
  for(auto x : atypical_elements_values_and_derivatives_) {
    // Only works for distribution data
    auto verticalLine = new QCPItemLine(ui->widget_plot);
    verticalLine->start->setCoords(x.first, 0);
    verticalLine->end->setCoords(x.first, -quantile_estimator_value_);
    if(x.second > 0)
      verticalLine->setPen(QPen(Qt::green));
    else
      verticalLine->setPen(QPen(Qt::red));
    lines_on_plot_.push_back(verticalLine);
  }

  return atypical_elements_values_and_derivatives_.size();
}

QString MainWindow::FormatNumberForDisplay(double number) {
  // According to PK the number should be displayed as #.######
  QString result = " ";

  if(number < 0) result = "";

  QStringList splitNumber = QString::number(number, 'f', 7).split(".");
  result += splitNumber[0];

  if(splitNumber.size() == 1) return result;

  result += ".";

  for(int i = 0; i < 6 && i < splitNumber[1].size(); ++i)
    result += splitNumber[1][i];

  return result;
}

void MainWindow::FillStandardDeviations(vector<std::shared_ptr<vector<double>>> *stDevs) {
  int dimensionsNumber = ui->spinBox_dimensionsNumber->value(),
      targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

  for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex) {
    stDevs->push_back(std::make_shared<vector<double>>());

    for(int dimensionIndex = 0; dimensionIndex < dimensionsNumber; ++dimensionIndex) {
      stDevs->back().get()->push_back
          (
              (dynamic_cast<QLineEdit *>(
                  (dynamic_cast<QTableWidget *>(ui->tableWidget_targetFunctions
                                                 ->cellWidget(functionIndex, static_cast<int>(TargetFunctionSettingsColumns::kStandardDeviationColumnIndex))))
                      ->cellWidget(dimensionIndex, 0)
              ))
                  ->text().toDouble()
          );
    }
  }
}

void MainWindow::FillMeans(vector<std::shared_ptr<vector<double>>> *means) {
  int dimensionsNumber = ui->spinBox_dimensionsNumber->value(),
      targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

  for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex) {
    means->push_back(std::make_shared<vector<double>>());

    for(int dimensionIndex = 0; dimensionIndex < dimensionsNumber; ++dimensionIndex) {
      means->back()->push_back
          (
              (dynamic_cast<QLineEdit *>(
                  (dynamic_cast<QTableWidget *>(ui->tableWidget_targetFunctions
                                                 ->cellWidget(functionIndex, static_cast<int>(TargetFunctionSettingsColumns::kMeanColumnIndex))))
                      ->cellWidget(dimensionIndex, 0)
              ))
                  ->text().toDouble()
          );
    }
  }
}

void MainWindow::FillDomain(QVector<std::shared_ptr<point>> *domain, std::shared_ptr<point> *prototypePoint) {
  // Check if domain is null pointer
  if(domain == nullptr) return;

  std::shared_ptr<point> pPoint;

  // Check if prototype is a null pointer
  if(prototypePoint == nullptr) {
    // If so make it a point pointer
    pPoint = std::make_shared<point>();
  }
  else {
    pPoint = *prototypePoint;
  }

  qreal val = ui->lineEdit_minX->text().toDouble();
  //qreal maxVal = 3 + ui->lineEdit_distributionProgression->text().toDouble() * 3000;  // Traveling case
  //qreal maxVal = 3 + ui->lineEdit_distributionProgression->text().toDouble(); /* should jump */
  qreal maxVal = ui->lineEdit_maxX->text().toDouble();

  //while(val <= ui->lineEdit_maxX->text().toDouble())
  while(val <= maxVal) {
    pPoint->push_back(val);

    if(pPoint->size() == (size_t) ui->spinBox_dimensionsNumber->value()) {
      domain->append(std::make_shared<point>());

      for(auto dimensionVal : *(pPoint)){
          domain->back()->push_back(dimensionVal);
      }
    }
    else {
      FillDomain(domain, prototypePoint);
    }

    pPoint->erase(pPoint->end() - 1);

    val += ui->lineEdit_domainDensity->text().toDouble();
  }
}

distribution *MainWindow::GenerateTargetDistribution(
    vector<std::shared_ptr<vector<double>>> *means,
    vector<std::shared_ptr<vector<double>>> *stDevs) {
  int seed = ui->lineEdit_seed->text().toInt();
  vector<double> contributions;
  vector<std::shared_ptr<distribution>> elementalDistributions;

  //int targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();
  int targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

  //double maxMean = ui->lineEdit_distributionProgression->text().toDouble() *  3000;
  double maxMean = ui->lineEdit_maxX->text().toDouble();

  for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex) {
    contributions.push_back
                     (
                         (dynamic_cast<QLineEdit *>(ui->tableWidget_targetFunctions
                                                     ->cellWidget(functionIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
                             ->text().toDouble()
                     );

    elementalDistributions.push_back(
        std::shared_ptr<distribution>(
            new normalDistribution(seed,
                                   (*means)[functionIndex].get(),
                                   (*stDevs)[functionIndex].get(),
                                   maxMean)));
  }

  return new complexDistribution(seed, &elementalDistributions, &contributions);
}

reservoirSamplingAlgorithm *MainWindow::GenerateReservoirSamplingAlgorithm(dataReader *reader,
                                                                           dataParser *parser) {
  int sampleSize = ui->lineEdit_sampleSize->text().toInt(),
      stepsNumber = ui->lineEdit_iterationsNumber->text().toInt(),
      samplingAlgorithmID = ui->comboBox_samplingAlgorithm->currentIndex();

  switch(samplingAlgorithmID) {
    case static_cast<int>(ReservoirSamplingAlgorithms::kBiasedReservoirSamplingAlgorithm):
      return new biasedReservoirSamplingAlgorithm(reader, parser, sampleSize, stepsNumber);
    case static_cast<int>(ReservoirSamplingAlgorithms::kBasicReservoirSamplingAlgorithm):
    default:
      return new basicReservoirSamplingAlgorithm(reader, parser, sampleSize, stepsNumber);
  }
}

kernelDensityEstimator *MainWindow::GenerateKernelDensityEstimator(
    int dimensionsNumber) {
  vector<int> kernelsIDs;
  vector<double> smoothingParameters;
  vector<std::string> carriersRestrictions;



  for(int rowNumber = 0; rowNumber < dimensionsNumber; ++rowNumber) {
    kernelsIDs.push_back(
        (dynamic_cast<QComboBox *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kKernelColumnIndex))))
            ->currentIndex());
    smoothingParameters.push_back((dynamic_cast<QLineEdit *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber,
                                                                                                         static_cast<int>(KernelSettingsColumns::kCarrierRestrictionColumnIndex))))
                                      ->text().toDouble());
    carriersRestrictions.push_back((dynamic_cast<QLineEdit *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber,
                                                                                                          static_cast<int>(KernelSettingsColumns::kCarrierRestrictionColumnIndex))))
                                       ->text().toStdString());
  }

  return new kernelDensityEstimator(
      &samples_,
      &smoothingParameters,
      &carriersRestrictions,
      PRODUCT,
      &kernelsIDs
                                   );
}

function *MainWindow::GenerateTargetFunction(
    vector<std::shared_ptr<vector<double>>> *means,
    vector<std::shared_ptr<vector<double>>> *stDevs) {
  vector<double> contributions;
  vector<std::shared_ptr<function>> elementalFunctions;
  int targetFunctionElementsNumber = ui->tableWidget_targetFunctions
                                       ->rowCount();

  // Check if contributions are set correctly. If they are, then last
  // contribution is >= 0;
  if(dynamic_cast<QLineEdit *>(ui->tableWidget_targetFunctions
                                ->cellWidget(targetFunctionElementsNumber - 1, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))
     )->text().toDouble() <= 0
      ) {
    // If not then uniform distributions and log error
    log("Contributions aren't set correctly. Uniforming contributions.");
    UniformContributions();
  }

  for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex) {

    contributions.push_back
                     (
                         (dynamic_cast<QLineEdit *>(ui->tableWidget_targetFunctions
                                                     ->cellWidget(functionIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
                             ->text().toDouble()
                     );

    elementalFunctions.push_back(
        std::shared_ptr<function>(new multivariateNormalProbabilityDensityFunction(means->at(functionIndex).get(),
                                                                                   stDevs->at(functionIndex).get())));
  }

  return new complexFunction(&contributions, &elementalFunctions);
}

int MainWindow::CanAnimationBePerformed(int dimensionsNumber) {
  if(dimensionsNumber == 1){
    return 1;
  }
  else{
    log("Dimensions number is not equal 1. Animation cannot be performed.");
  }
  return -2;
}

void MainWindow::on_spinBox_dimensionsNumber_editingFinished() {
  RefreshKernelsTable();
  RefreshTargetFunctionTable();
}

void MainWindow::RefreshKernelsTable() {
  // Get new number of rows
  int newNumberOfRows = ui->spinBox_dimensionsNumber->value();

  // If new number of rows is equal to current number of rows do nothing
  if(newNumberOfRows == ui->tableWidget_dimensionKernels->rowCount()) return;

  // Set new row count
  ui->tableWidget_dimensionKernels->setRowCount(newNumberOfRows);

  QLocale locale = QLocale::English;
  locale.setNumberOptions(QLocale::c().numberOptions());

  auto smoothingParameterValidator = new QDoubleValidator(kMinSmoothingParameter, kMaxSmoothingParameter, kDecimalNumbers, this);
  smoothingParameterValidator->setLocale(locale);
  smoothingParameterValidator->setNotation(QDoubleValidator::StandardNotation);

  for(int rowNumber = 0; rowNumber < newNumberOfRows; ++rowNumber)
    AddKernelToTable(rowNumber, smoothingParameterValidator);
}

void MainWindow::AddKernelToTable(int rowNumber,
                                  QDoubleValidator *smoothingParameterValidator) {
  // Add combobox with kernels
  ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kKernelColumnIndex), new QComboBox());

  (dynamic_cast<QComboBox *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kKernelColumnIndex))))
      ->insertItems(0, kernel_types_);

  // Add input box with validator for smoothing parameters
  ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kSmoothingParameterColumnIndex), new QLineEdit());

  (dynamic_cast<QLineEdit *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kSmoothingParameterColumnIndex))))
      ->setText("1.0");
  (dynamic_cast<QLineEdit *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kSmoothingParameterColumnIndex))))
      ->setValidator(smoothingParameterValidator);

  // Add input box for carrier restriction value
  ui->tableWidget_dimensionKernels->setCellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kCarrierRestrictionColumnIndex), new QLineEdit());

  (dynamic_cast<QLineEdit *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kCarrierRestrictionColumnIndex))))
      ->setText("None.");
}

void MainWindow::RefreshTargetFunctionTable() {
  int numberOfRows = ui->tableWidget_targetFunctions->rowCount(),
      dimensionsNumber = ui->spinBox_dimensionsNumber->value();

  // Ensure that rows number is at least 1
  if(numberOfRows == 0)
    numberOfRows = 1;

  // Set row count
  ui->tableWidget_targetFunctions->setRowCount(numberOfRows);

  QLocale locale = QLocale::English;
  locale.setNumberOptions(QLocale::c().numberOptions());

  auto meanValidator = new QDoubleValidator(-10.0, 10.0, 3, this);
  meanValidator->setLocale(locale);
  meanValidator->setNotation(QDoubleValidator::StandardNotation);

  auto stDevValidator = new QDoubleValidator(-5.0, 5.0, 3, this);
  stDevValidator->setLocale(locale);
  stDevValidator->setNotation(QDoubleValidator::StandardNotation);

  auto contributionValidator = new QDoubleValidator(0.0, 100.0, 3, this);
  contributionValidator->setLocale(locale);
  contributionValidator->setNotation(QDoubleValidator::StandardNotation);

  QTableWidget *targetFunctionTablePointer = static_cast<QTableWidget *>(ui->tableWidget_targetFunctions),
      *meansTablePointer, *stDevsTablePointer;

  for(int rowIndex = 0; rowIndex < numberOfRows; ++rowIndex) {
    // TODO TR: Ensure that this doesn't result in memory leaks
    targetFunctionTablePointer->setCellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kMeanColumnIndex), new QTableWidget());

    meansTablePointer =
        dynamic_cast<QTableWidget *>(ui->tableWidget_targetFunctions->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kMeanColumnIndex)));
    meansTablePointer->setRowCount(dimensionsNumber);
    meansTablePointer->setColumnCount(1);
    meansTablePointer->horizontalHeader()->hide();

    // TODO TR: Ensure that this doesn't result in memory leaks
    targetFunctionTablePointer->setCellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kStandardDeviationColumnIndex), new QTableWidget());

    stDevsTablePointer =
        dynamic_cast<QTableWidget *>(ui->tableWidget_targetFunctions->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kStandardDeviationColumnIndex)));
    stDevsTablePointer->setRowCount(dimensionsNumber);
    stDevsTablePointer->setColumnCount(1);
    stDevsTablePointer->horizontalHeader()->hide();

    for(int dimensionNumber = 0; dimensionNumber < dimensionsNumber; ++dimensionNumber) {
      meansTablePointer->setCellWidget(dimensionNumber, 0, new QLineEdit());
      (dynamic_cast<QLineEdit *>(meansTablePointer->cellWidget(dimensionNumber, 0)))->setText("0.0");
      (dynamic_cast<QLineEdit *>(meansTablePointer->cellWidget(dimensionNumber, 0)))->setValidator(meanValidator);

      stDevsTablePointer->setCellWidget(dimensionNumber, 0, new QLineEdit());
      (dynamic_cast<QLineEdit *>(stDevsTablePointer->cellWidget(dimensionNumber, 0)))->setText("1.0");
      (dynamic_cast<QLineEdit *>(stDevsTablePointer->cellWidget(dimensionNumber, 0)))->setValidator(stDevValidator);
    }

    // TODO TR: Ensure that this doesn't result in memory leaks
    targetFunctionTablePointer->setCellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex), new QLineEdit());
    (dynamic_cast<QLineEdit *>(targetFunctionTablePointer->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
        ->setMaxLength(6);
    (dynamic_cast<QLineEdit *>(targetFunctionTablePointer->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
        ->setValidator(contributionValidator);
    QObject::connect(
        (dynamic_cast<QLineEdit *>(targetFunctionTablePointer->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex)))),
        SIGNAL(textEdited(QString)), this, SLOT(UpdateLastContribution()));
  }

  // Disable last contribution cell, as it's filled automatically
  (dynamic_cast<QLineEdit *>(targetFunctionTablePointer->cellWidget(numberOfRows - 1, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
      ->setEnabled(false);

  UniformContributions();
}

void MainWindow::UniformContributions() {
  int numberOfRows = ui->tableWidget_targetFunctions->rowCount(),
      lastRowIndex = numberOfRows - 1;

  for(int rowIndex = 0; rowIndex < lastRowIndex; ++rowIndex) {
    dynamic_cast<QLineEdit *>(
        ui->tableWidget_targetFunctions
          ->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))
    )->setText(QString::number(100.0 / numberOfRows));
  }

  dynamic_cast<QLineEdit *>(
      ui->tableWidget_targetFunctions
        ->cellWidget(lastRowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))
  )->setText(QString::number(CountLastContribution()));
}

qreal MainWindow::CountLastContribution() {
  qreal result = 100.0;
  int lastRowIndex = ui->tableWidget_targetFunctions->rowCount() - 1;

  for(int rowIndex = 0; rowIndex < lastRowIndex; ++rowIndex)
    result -= (dynamic_cast<QLineEdit *>
    (ui->tableWidget_targetFunctions
       ->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex)))
    )->text().toDouble();

  return result;
}

void MainWindow::UpdateLastContribution() {
  int lastRowIndex = ui->tableWidget_targetFunctions->rowCount() - 1;
  qreal lastContributionValue = CountLastContribution();

  (dynamic_cast<QLineEdit *>(ui
      ->tableWidget_targetFunctions
      ->cellWidget(lastRowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex)))
  )->setText(QString::number(lastContributionValue));
}

void MainWindow::on_pushButton_addTargetFunction_clicked() {
  int newRowsNumber = ui->tableWidget_targetFunctions->rowCount() + 1;

  // Set row count
  ui->tableWidget_targetFunctions->setRowCount(newRowsNumber);

  RefreshTargetFunctionTable();
}

void MainWindow::on_pushButton_removeTargetFunction_clicked() {
  int newRowsNumber = ui->tableWidget_targetFunctions->rowCount() - 1;

  // Set row count
  ui->tableWidget_targetFunctions->setRowCount(newRowsNumber);

  RefreshTargetFunctionTable();
}

void MainWindow::on_pushButton_start_clicked() {
  //Run1DExperimentWithDESDA();
  Run1DExperimentWithClusterKernels();
}

void MainWindow::on_pushButton_clicked() {
  log("2D Experiment start.");

  screen_generation_frequency_ = 1;
  int seed = ui->lineEdit_seed->text().toInt();
  int m0 = ui->lineEdit_sampleSize->text().toInt();

  // Prepare image location.
  QString expNum = "1369 (2D)";
  this->setWindowTitle("Experiment #" + expNum);
  QString expDesc =
      "iw=" + QString::number(screen_generation_frequency_)
      + ", euclidean KPSS, v=tor dla 2D, me=1k, m0=4k, mMin=400, sz261";
  //QString driveDir = "\\\\beabourg\\private\\"; // WIT PCs
  QString driveDir = "Y:\\"; // WIT PCs after update
  //QString driveDir = "D:\\Test\\"; // Home
  //QString driveDir = "d:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\";
  QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " ("
                    + expDesc + ")\\";
  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  // Contour levels calculation.
  QList<double> contourLevels;
  double level = 0.05;
  while(level < 0.21){
    contourLevels += level;
    level += 0.05;
  }

  // Add clusters_ to the estimator
  means_ = {std::make_shared<std::vector<double >>()};
  means_.back()->push_back(0);
  means_.back()->push_back(0);

  standard_deviations_ = {std::make_shared<std::vector<double >>()};
  standard_deviations_.back()->push_back(1);
  standard_deviations_.back()->push_back(1);

  auto densityFunction =
      new multivariateNormalProbabilityDensityFunction(means_.back().get(), standard_deviations_.back().get());
  contour_plot_->addQwtPlotSpectrogram(new SpectrogramData2(densityFunction, -10.0), QPen(QColor(255, 0, 0)));

  // Create estimator object
  std::shared_ptr<kernelDensityEstimator>
      estimator(GenerateKernelDensityEstimator(2));

  estimator->_shouldConsiderWeights = true;

  std::vector<double> pt = {0, 0};

  contour_plot_->addQwtPlotSpectrogram(new SpectrogramData2(estimator.get(), -10.0), QPen(QColor(0, 255, 0)));

  // After adding plots set contours and stuff.
  contour_plot_->setContours(contourLevels);
  contour_plot_->showContour(true);
  contour_plot_->setAlpha(0);

  // Set limit on axes.
  contour_plot_->setAxesLimit(5);

  std::shared_ptr < distribution >
  targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));

  std::vector<double> meansForDistribution = {0.0, 0.0};
  std::vector<double> stDevsForDistribution = {1.0, 1.0};

  parser_.reset(new distributionDataParser(&attributes_data_));

  reader_.reset(
      new progressiveDistributionDataReader(targetDistribution.get(), 0,
                                            0,  /* Delay */
                                            new normalDistribution(0, &meansForDistribution, &stDevsForDistribution,
                                                                   55))
               );

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  reservoirSamplingAlgorithm *samplingAlgorithm =
      GenerateReservoirSamplingAlgorithm(reader_.get(), parser_.get());

  objects_.clear();
  clusters_ = &stored_medoids_;
  clusters_->clear();

  int pluginRank = 3;
  double newWeightB = 0.5;
  groupingThread gt(&stored_medoids_, parser_);

  derivative_estimator_.reset(GenerateKernelDensityEstimator(2));
  enhanced_kde_.reset(GenerateKernelDensityEstimator(2));

  DESDA DESDAAlgorithm(
      estimator,
      derivative_estimator_,
      enhanced_kde_,
      ui->lineEdit_weightModifier->text().toDouble(),
      samplingAlgorithm,
      clusters_,
      &stored_medoids_,
      ui->lineEdit_rarity->text().toDouble(),
      &gt, newWeightB, pluginRank
                      );

  // Start the test
  step_number_ = 0;

  time_t startTime, endTime;

  l1_n_ = 0;
  l2_n_ = 0;
  sup_n_ = 0;
  mod_n_ = 0;
  double actual_l1 = 0;
  double actual_l2 = 0;
  double actual_sup = 0;
  double actual_mod = 0;

  int errorCalculationsNumber = 0;
  double sum_l1 = 0, sum_l2 = 0, sum_sup = 0, sum_mod = 0;
  QwtContourPlotUI plotUi(&step_number_, screen_generation_frequency_, seed,
                          &DESDAAlgorithm, &l1_n_, &l2_n_, &sup_n_, &mod_n_,
                          &actual_l1, &actual_l2, &actual_sup, &actual_mod);
  plotUi.attach(contour_plot_);
  plotUi.updateTexts();
  //QVector<int> initialDrawingSteps = {1, 2, 3, 4, 5, 6, 7, 8, 9 , 10};
  QVector<int> initialDrawingSteps = {};

  std::vector<double> model_function_values = {};
  std::vector<double> estimator_values = {};
  double domain_area = 0;
  std::vector<std::vector<double>> error_domain = {};

  ErrorsCalculator errors_calculator(
    &model_function_values, &estimator_values, &error_domain, &domain_area
  );

  log("Experiment started.");
  for(step_number_ = 1; step_number_ < 13001; ++step_number_) {

    log("New step.");
    startTime = time(nullptr);

    log("Performing new step.");
    DESDAAlgorithm.performStep();
    log("Step performed.");

    // Drawing
    if(step_number_ % screen_generation_frequency_ == 0 ||
       initialDrawingSteps.contains(step_number_)) {

      log("Drawing started.");

      log("Estimator preparation.");
      DESDAAlgorithm.prepareEstimatorForContourPlotDrawing();
      log("Estimator preparation finished.");
      // Error calculation

      if(step_number_ >= 0) {
        log("Error calculation started.");
        ++errorCalculationsNumber;
        error_domain = Generate2DPlotErrorDomain(&DESDAAlgorithm);
        domain_area = Calculate2DDomainArea(error_domain);
        model_function_values =
            GetFunctionsValueOnDomain(densityFunction, error_domain);
        estimator_values =
            GetFunctionsValueOnDomain(estimator.get(), error_domain);

        actual_l1 = errors_calculator.CalculateL1Error();
        actual_l2 = errors_calculator.CalculateL2Error();
        actual_sup = errors_calculator.CalculateSupError();
        actual_mod = errors_calculator.CalculateModError();
        sum_l1 += actual_l1;
        sum_l2 += actual_l2;
        sum_sup += actual_sup;
        sum_mod = actual_mod;

        l1_n_ = sum_l1 / errorCalculationsNumber;
        l2_n_ = sum_l2 / errorCalculationsNumber;
        sup_n_ = sum_sup / errorCalculationsNumber;
        mod_n_ = sum_mod / errorCalculationsNumber;
        log("Error calculation finished.");
      }

      densityFunction->setMeans(*means_.back().get());

      log("Texts updates.");
      plotUi.updateTexts();

      log("Replotting.");
      contour_plot_->replot();

      log("Restoring weights.");
      DESDAAlgorithm.restoreClustersCWeights();

      endTime = time(nullptr);

      log("Processing.");

      QCoreApplication::processEvents();

      QString imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image name: " + imageName);
      log("Saved: " + QString(ui->widget_contour_plot_holder->grab().save(imageName)));
      log("Drawing finished.");
    }

    if(step_number_ == 10) {
      initialDrawingSteps.clear();  // To reduce comparisons for drawing.
    }

    endTime = time(nullptr);

    log("Step time: " + QString::number(endTime - startTime) + " s");
  }

  log("Done!");
}

void MainWindow::resizeEvent(QResizeEvent *event) {
  int offset = 10; // Offset in px, so that scale is in
  QMainWindow::resizeEvent(event);
  int newSize = std::min(ui->widget_contour_plot->height(),
                         ui->widget_contour_plot->width()) - offset;
  //ui->widget->resize(newSize, newSize);
  contour_plot_->resize(2 * newSize, newSize);
}

double MainWindow::Calculate2DDomainArea(const std::vector<std::vector<double>> &domain) {
  // Given how 2D domain is calculated, it's easy to see, that only first and
  // last point of the domain are of interest (last holding max values, and
  // first holding minima).
  auto xLen = domain[domain.size() - 1][0] - domain[0][0];
  auto yLen = domain[domain.size() - 1][1] - domain[0][1];

  return xLen * yLen;
}

std::vector<double> MainWindow::GetFunctionsValueOnDomain(function *func, const std::vector<std::vector<double>> &domain) {
  std::vector<double> values = {};

  for(auto pt : domain) {
    values.push_back(func->getValue(&pt));
  }

  return values;
}

std::vector<std::vector<double>> MainWindow::Generate2DPlotErrorDomain(DESDA *DESDAAlgorithm) {
  std::vector<point> domainValues = {};
  auto xDomainValues = DESDAAlgorithm->getErrorDomain(0);
  auto yDomainValues = DESDAAlgorithm->getErrorDomain(1);

  for(auto x : xDomainValues) {
    for(auto y : yDomainValues) {
      domainValues.push_back({x, y});
    }
  }

  return domainValues;
}

std::vector<std::vector<double>> MainWindow::Generate1DPlotErrorDomain(
    DESDA *DESDAAlgorithm) {
  std::vector<point> domainValues = {};
  auto xDomainValues = DESDAAlgorithm->getErrorDomain(0);

  for(auto x : xDomainValues) {
    domainValues.push_back({x});
  }

  return domainValues;
}

std::vector<std::vector<double>> MainWindow::Generate1DWindowedPlotErrorDomain(
    DESDA *DESDAAlgorithm) {
  std::vector<point> domainValues = {};
  auto xDomainValues = DESDAAlgorithm->getWindowedErrorDomain();

  for(auto x : xDomainValues) {
    domainValues.push_back({x});
  }

  return domainValues;
}

void MainWindow::Run1DExperimentWithDESDA() {
  int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

  if(!CanAnimationBePerformed(dimensionsNumber)) return;

  QString seedString = ui->lineEdit_seed->text();

  // Log that application started generating KDE
  // Standard seed was 5625.
  log("KDE animation started.");
  log("Seed: " + seedString);
  log("Sample size: " + ui->lineEdit_sampleSize->text());

  step_number_ = 0;

  srand(static_cast<unsigned int>(seedString.toInt()));

  FillMeans(&means_);
  FillStandardDeviations(&standard_deviations_);

  target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

  std::shared_ptr<kernelDensityEstimator>
      estimator(GenerateKernelDensityEstimator(dimensionsNumber));

  estimator->_shouldConsiderWeights = false;

  derivative_estimator_.reset(GenerateKernelDensityEstimator(dimensionsNumber));
  enhanced_kde_.reset(GenerateKernelDensityEstimator(dimensionsNumber));

  std::shared_ptr<distribution>
      targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));
  vector<double> alternativeDistributionMean = {0.0};
  vector<double> alternativeDistributionStDevs = {1.0};
  qreal progressionSize =
      ui->lineEdit_distributionProgression->text().toDouble();

  parser_.reset(new distributionDataParser(&attributes_data_));

  reader_.reset(
      new progressiveDistributionDataReader(targetDistribution.get(),
                                            progressionSize,
                                            0,  /* Delay */
                                            new normalDistribution(seedString.toInt(), &alternativeDistributionMean,
                                                                   &alternativeDistributionStDevs, 55))
               );

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  reservoirSamplingAlgorithm *algorithm =
      GenerateReservoirSamplingAlgorithm(reader_.get(), parser_.get());

  objects_.clear();

  int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();
  int medoidsNumber = 50;
  groupingThread gt(&stored_medoids_, parser_);

  gt.setAttributesData(&attributes_data_);

  log("Attributes data set.");

  int sampleSize = ui->lineEdit_sampleSize->text().toInt();
  gt.initialize(medoidsNumber, sampleSize);

  double newWeightB = 0.5;

  clusters_ = &stored_medoids_;

  weightedSilvermanSmoothingParameterCounter smoothingParamCounter(clusters_, 0);

  derivative_estimator_->_shouldConsiderWeights = false;

  int pluginRank = 3;
  DESDA DESDAAlgorithm(
      estimator,
      derivative_estimator_,
      enhanced_kde_,
      ui->lineEdit_weightModifier->text().toDouble(),
      algorithm,
      clusters_,
      &stored_medoids_,
      ui->lineEdit_rarity->text().toDouble(),
      &gt, newWeightB, pluginRank
                      );
  QString expNum = "1281";
  this->setWindowTitle("Experiment #" + expNum);
  QString expDesc = "reservoir, plugin " + QString::number(pluginRank) +
                    ", 0,2N(-5,1)0,4N(0,1)0,4N(5,1), v=0-1-0-1, m0=" + QString::number(DESDAAlgorithm._maxM) +
                    ", mMin=" + QString::number(DESDAAlgorithm._minM) +
                    //", mKPSS=" + QString::number(DESDAAlgorithm._kpssM) + // This can just be commented out.
                    ", sz475";
  screen_generation_frequency_ = 10;

  //QString driveDir = "\\\\beabourg\\private\\"; // WIT PCs
  QString driveDir = "D:\\Test\\"; // Home
  QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " ("
                    + expDesc + ")\\";

  ClearPlot();
  ResizePlot();

  // Initial screen should only contain exp number (as requested).
  plotLabel expNumLabel(ui->widget_plot, 0.02, 0.25, "Exp." + expNum);
  expNumLabel.setFont(QFont("Courier New", 250));

  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  QString imageName = dirPath + QString::number(0) + ".png";

  log("Image saved: " + QString(ui->widget_plot->savePng(imageName,0, 0, 1, -1)));
  expNumLabel.setText("");

  QVector<std::shared_ptr<plotLabel>> plotLabels = {};
  double horizontalOffset = 0.01, verticalOffset = 0.01, verticalStep = 0.03;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "i     = ", &step_number_,
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  verticalOffset += verticalStep;


  plotLabels.push_back(std::make_shared<plotLabel>(MainWindow::ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "iw    = "
                                                                                     + QString::number(
                                                                                         screen_generation_frequency_)));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "seed  = " + seedString));
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel KPSSTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                          "KPSS     = 0");
  verticalOffset += verticalStep;

  plotLabel sgmKPSSTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                             "sgmKPSS  = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel mKPSSTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                           "mKPSS = " + QString::number(DESDAAlgorithm._kpssM));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot, horizontalOffset, verticalOffset,
                                                   "m0    = " + ui->lineEdit_sampleSize->text()));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "mmin  = "
                                                                                     + QString::number(
                                                                                         DESDAAlgorithm._minM)));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "m     = ", &(DESDAAlgorithm._m),
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel betaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                          "beta0 = " + QString::number(DESDAAlgorithm._beta0));
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel rTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "r     =" + FormatNumberForDisplay(DESDAAlgorithm._r));
  verticalOffset += verticalStep;

  plotLabel qTextLabel(ui->widget_plot, horizontalOffset, verticalOffset,
                       "q     =" + FormatNumberForDisplay(
                           DESDAAlgorithm._quantileEstimator));

  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "rare  = ",
                                                   &(DESDAAlgorithm._rareElementsNumber),
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  //====================  SECOND COLUMN =================//

  horizontalOffset = 0.20;
  verticalOffset = 0.01 + 9 * verticalStep;

  //==================== ERRORS SUM =================//

  horizontalOffset = 0.87;
  verticalOffset = 0.01;

  plotLabel L1WTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_w = 0");
  verticalOffset += verticalStep;

  plotLabel L1MTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_m = 0");
  verticalOffset += verticalStep;

  plotLabel L1DTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_d = 0");
  verticalOffset += verticalStep;

  plotLabel L1PTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_p = 0");
  verticalOffset += verticalStep;

  plotLabel L1NTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_n = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel L2WTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_w = 0");
  verticalOffset += verticalStep;

  plotLabel L2MTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_m = 0");
  verticalOffset += verticalStep;

  plotLabel L2DTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_d = 0");
  verticalOffset += verticalStep;

  plotLabel L2PTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_p = 0");
  verticalOffset += verticalStep;

  plotLabel L2NTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_n = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel supWTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_w = 0");
  verticalOffset += verticalStep;

  plotLabel supMTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_m = 0");
  verticalOffset += verticalStep;

  plotLabel supDTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_d = 0");
  verticalOffset += verticalStep;

  plotLabel supPTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_p = 0");
  verticalOffset += verticalStep;

  plotLabel supNTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_n = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel modWTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_w = 0");
  verticalOffset += verticalStep;

  plotLabel modMTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_m = 0");
  verticalOffset += verticalStep;

  plotLabel modDTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_d = 0");
  verticalOffset += verticalStep;

  plotLabel modPTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_p = 0");
  verticalOffset += verticalStep;

  plotLabel modNTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_n = 0");
  verticalOffset += verticalStep;

  FillDomain(&domain_, nullptr);
  for(const auto& pt : domain_) drawable_domain_.push_back(pt->at(0));

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

  int numberOfErrorCalculations = 1;
  QVector<int> additionalScreensSteps = {};

  /*
  for(int i = 990; i < 1011; ++i){
      additionalScreensSteps.append(i);
  }
  */

  double error_domain_length = 0;
  double windowed_error_domain_length = 0;
  std::vector<std::vector<double>> error_domain = {};
  std::vector<std::vector<double>> windowed_error_domain = {};

  std::vector<double> windowed_model_values = {};
  std::vector<double> windowed_kde_values = {};

  std::vector<double> model_values = {};
  std::vector<double> less_elements_kde_values = {};
  std::vector<double> weighted_kde_values = {};
  std::vector<double> enhanced_kde_values = {};
  std::vector<double> rare_elements_kde_values = {};

  ErrorsCalculator windowed_errors_calculator(
      &windowed_model_values, &windowed_kde_values, &windowed_error_domain, &windowed_error_domain_length
                                             );
  ErrorsCalculator less_elements_kde_errors_calculator(
      &model_values, &less_elements_kde_values, &error_domain, &error_domain_length
                                                      );
  ErrorsCalculator weighted_kde_errors_calculator(
      &model_values, &weighted_kde_values, &error_domain, &error_domain_length
                                                 );
  ErrorsCalculator enhanced_kde_errors_calculator(
      &model_values, &enhanced_kde_values, &error_domain, &error_domain_length
                                                 );
  ErrorsCalculator rare_elements_kde_errors_calculator(
      &model_values, &rare_elements_kde_values, &error_domain, &error_domain_length
                                                      );

  for(step_number_ = 1; step_number_ < stepsNumber; ++step_number_) {
    clock_t executionStartTime = clock();

    DESDAAlgorithm.performStep();

    target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

    if(step_number_ % screen_generation_frequency_ == 0 || step_number_ < 10
       || additionalScreensSteps.contains(step_number_)) {
      log("Drawing in step number " + QString::number(step_number_) + ".");

      kernel_prognosis_derivative_values_ =
          DESDAAlgorithm.getKernelPrognosisDerivativeValues(&drawable_domain_);

      // Error calculations
      if(step_number_ >= 1) {

        log("Getting windowed domain.");
        //windowed_error_domain_ = DESDAAlgorithm.getWindowedErrorDomain();
        windowed_error_domain = Generate1DWindowedPlotErrorDomain(&DESDAAlgorithm);
        log("Getting non-windowed domain.");
        //error_domain_ = DESDAAlgorithm.getErrorDomain();
        error_domain = Generate1DPlotErrorDomain(&DESDAAlgorithm);

        log("Getting model plot on windowed.");
        windowed_model_values = GetFunctionsValueOnDomain(target_function_.get(), windowed_error_domain);
        log("Getting KDE plot on windowed.");
        windowed_kde_values = DESDAAlgorithm.getWindowKDEValues(&windowed_error_domain);

        log("Getting model plot.");
        model_values = GetFunctionsValueOnDomain(target_function_.get(), error_domain);
        log("Getting KDE plot on lesser elements.");
        //less_elements_estimator_error_y_ = DESDAAlgorithm.getKDEValues(&error_domain_);
        less_elements_kde_values = DESDAAlgorithm.getKDEValues(&error_domain);
        log("Getting weighted KDE plot.");
        //weighted_estimator_error_y_ = DESDAAlgorithm.getWeightedKDEValues(&error_domain_);
        weighted_kde_values = DESDAAlgorithm.getWeightedKDEValues(&error_domain);
        log("Getting sgm KDE plot.");
        //sigmoidally_enhanced_error_plot_y_ = DESDAAlgorithm.getEnhancedKDEValues(&error_domain_);
        enhanced_kde_values = DESDAAlgorithm.getEnhancedKDEValues(&error_domain);
        log("Getting rare KDE plot.");
        //rare_elements_enhanced_error_plot_Y = DESDAAlgorithm.getRareElementsEnhancedKDEValues(&error_domain_);
        rare_elements_kde_values = DESDAAlgorithm.getRareElementsEnhancedKDEValues(&error_domain);

        error_domain_length =
            error_domain[error_domain.size() - 1][0] - error_domain[0][0];
        windowed_error_domain_length =
            windowed_error_domain[windowed_error_domain.size() - 1][0] - windowed_error_domain[0][0];

        l1_w_ += windowed_errors_calculator.CalculateL1Error();
        l2_w_ += windowed_errors_calculator.CalculateL2Error();
        sup_w_ += windowed_errors_calculator.CalculateSupError();
        mod_w_ += windowed_errors_calculator.CalculateModError();

        l1_m_ += less_elements_kde_errors_calculator.CalculateL1Error();
        l2_m_ += less_elements_kde_errors_calculator.CalculateL2Error();
        sup_m_ += less_elements_kde_errors_calculator.CalculateSupError();
        mod_m_ += less_elements_kde_errors_calculator.CalculateModError();

        l1_d_ += weighted_kde_errors_calculator.CalculateL1Error();
        l2_d_ += weighted_kde_errors_calculator.CalculateL2Error();
        sup_d_ += weighted_kde_errors_calculator.CalculateSupError();
        mod_d_ += weighted_kde_errors_calculator.CalculateModError();

        l1_p_ += enhanced_kde_errors_calculator.CalculateL1Error();
        l2_p_ += enhanced_kde_errors_calculator.CalculateL2Error();
        sup_p_ += enhanced_kde_errors_calculator.CalculateSupError();
        mod_p_ += enhanced_kde_errors_calculator.CalculateModError();

        l1_n_ += rare_elements_kde_errors_calculator.CalculateL1Error();
        l2_n_ += rare_elements_kde_errors_calculator.CalculateL2Error();
        sup_n_ += rare_elements_kde_errors_calculator.CalculateSupError();
        mod_n_ += rare_elements_kde_errors_calculator.CalculateModError();

        ++numberOfErrorCalculations;
      }

      // ============ SUMS =========== //

      L1WTextLabel
          .setText("L1_w  =" + FormatNumberForDisplay(
              l1_w_ / numberOfErrorCalculations));
      L1MTextLabel
          .setText("L1_m  =" + FormatNumberForDisplay(
              l1_m_ / numberOfErrorCalculations));
      L1DTextLabel
          .setText("L1_d  =" + FormatNumberForDisplay(
              l1_d_ / numberOfErrorCalculations));
      L1PTextLabel
          .setText("L1_p  =" + FormatNumberForDisplay(
              l1_p_ / numberOfErrorCalculations));
      L1NTextLabel
          .setText("L1_n  =" + FormatNumberForDisplay(
              l1_n_ / numberOfErrorCalculations));
      L2WTextLabel
          .setText("L2_w  =" + FormatNumberForDisplay(
              l2_w_ / numberOfErrorCalculations));
      L2MTextLabel
          .setText("L2_m  =" + FormatNumberForDisplay(
              l2_m_ / numberOfErrorCalculations));
      L2DTextLabel
          .setText("L2_d  =" + FormatNumberForDisplay(
              l2_d_ / numberOfErrorCalculations));
      L2PTextLabel
          .setText("L2_p  =" + FormatNumberForDisplay(
              l2_p_ / numberOfErrorCalculations));
      L2NTextLabel
          .setText("L2_n  =" + FormatNumberForDisplay(
              l2_n_ / numberOfErrorCalculations));
      supWTextLabel
          .setText("sup_w =" + FormatNumberForDisplay(
              sup_w_ / numberOfErrorCalculations));
      supMTextLabel
          .setText("sup_m =" + FormatNumberForDisplay(
              sup_m_ / numberOfErrorCalculations));
      supDTextLabel
          .setText("sup_d =" + FormatNumberForDisplay(
              sup_d_ / numberOfErrorCalculations));
      supPTextLabel
          .setText("sup_p =" + FormatNumberForDisplay(
              sup_p_ / numberOfErrorCalculations));
      supNTextLabel
          .setText("sup_n =" + FormatNumberForDisplay(
              sup_n_ / numberOfErrorCalculations));
      modWTextLabel
          .setText("mod_w =" + FormatNumberForDisplay(
              mod_w_ / numberOfErrorCalculations));
      modMTextLabel
          .setText("mod_m =" + FormatNumberForDisplay(
              mod_m_ / numberOfErrorCalculations));
      modDTextLabel
          .setText("mod_d =" + FormatNumberForDisplay(
              mod_d_ / numberOfErrorCalculations));
      modPTextLabel
          .setText("mod_p =" + FormatNumberForDisplay(
              mod_p_ / numberOfErrorCalculations));
      modNTextLabel
          .setText("mod_n =" + FormatNumberForDisplay(
              mod_n_ / numberOfErrorCalculations));


      // ============= LEFT SIDE UPDATE ================ //

      qTextLabel.setText("q     =" + FormatNumberForDisplay(
          DESDAAlgorithm._quantileEstimator));

      KPSSTextLabel.setText("KPSS     = " + FormatNumberForDisplay(
          DESDAAlgorithm.getStationarityTestValue()));

      sgmKPSSTextLabel.setText("sgmKPSS  = " + FormatNumberForDisplay(
          DESDAAlgorithm._sgmKPSS));

      rTextLabel.setText("r     =" + FormatNumberForDisplay(DESDAAlgorithm._r));

      DrawPlots(&DESDAAlgorithm);

      betaTextLabel.setText("beta0 = " + QString::number(DESDAAlgorithm._beta0));

      for(const auto &label : plotLabels) label->updateText();

      ui->widget_plot->replot();
      QCoreApplication::processEvents();

      if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

      imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image saved: " + QString(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
    }
  }

  log("Animation finished.");
}

void MainWindow::Run1DExperimentWithClusterKernels() {
  int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

  if(!CanAnimationBePerformed(dimensionsNumber)) return;

  QString seedString = ui->lineEdit_seed->text();

  // Log that application started generating KDE
  // Standard seed was 5625.
  log("KDE animation with Cluster Kernels started.");
  log("Seed: " + seedString);
  log("Sample size: " + ui->lineEdit_sampleSize->text());

  int number_of_cluster_kernels = 1000;
  step_number_ = 0;

  srand(static_cast<unsigned int>(seedString.toInt()));

  // Creating target function.
  FillMeans(&means_);
  FillStandardDeviations(&standard_deviations_);
  target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

  std::shared_ptr<distribution>
      targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));
  vector<double> alternativeDistributionMean = {0.0};
  vector<double> alternativeDistributionStDevs = {1.0};
  qreal progressionSize =
      ui->lineEdit_distributionProgression->text().toDouble();

  parser_.reset(new distributionDataParser(&attributes_data_));

  reader_.reset(
      new progressiveDistributionDataReader(targetDistribution.get(),
                                            progressionSize,
                                            0,  /* Delay */
                                            new normalDistribution(seedString.toInt(), &alternativeDistributionMean,
                                                                   &alternativeDistributionStDevs, 55))
               );

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();

  log("Attributes data set.");

  int sampleSize = ui->lineEdit_sampleSize->text().toInt();

  QString expNum = "CKKDE_TEST_6";
  this->setWindowTitle("Experiment #" + expNum);
  QString expDesc = "Cluster Kernels, m = " + QString::number(number_of_cluster_kernels) + ", mean-var-resampling, weighted list-based algorithm, alpha=0.01";
  screen_generation_frequency_ = 10;

  //QString driveDir = "\\\\beabourg\\private\\"; // WIT PCs
  //QString driveDir = "D:\\Test\\"; // Home
  QString driveDir = "d:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\";
  QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " ("
                    + expDesc + ")\\";

  ClearPlot();
  ResizePlot();

  // Initial screen should only contain exp number (as requested).
  plotLabel expNumLabel(ui->widget_plot, 0.02, 0.25, "Exp." + expNum);
  expNumLabel.setFont(QFont("Courier New", 250));

  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  QString imageName = dirPath + QString::number(0) + ".png";

  log("Image saved: " + QString(ui->widget_plot->savePng(imageName,0, 0, 1, -1)));
  expNumLabel.setText("");

  QVector<std::shared_ptr<plotLabel>> plotLabels = {};
  double horizontalOffset = 0.01, verticalOffset = 0.01, verticalStep = 0.03;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "i     = ", &step_number_,
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  verticalOffset += verticalStep;


  plotLabels.push_back(std::make_shared<plotLabel>(MainWindow::ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "iw    = "
                                                                                     + QString::number(
                                                                                         screen_generation_frequency_)));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "seed  = " + seedString));
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;


  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot, horizontalOffset, verticalOffset,
                                                   "m0    = " + ui->lineEdit_sampleSize->text()));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "m     = ", &(number_of_cluster_kernels),
                                                   std::make_shared<plotLabelIntDataPreparator>()));

  //====================  SECOND COLUMN =================//

  horizontalOffset = 0.20;
  verticalOffset = 0.01 + 9 * verticalStep;

  //==================== ERRORS SUM =================//

  horizontalOffset = 0.87;
  verticalOffset = 0.01;

  plotLabel L1WTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1_w = 0");
  verticalOffset += verticalStep;

  plotLabel L2WTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2_w = 0");
  verticalOffset += verticalStep;

  plotLabel supWTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup_w = 0");
  verticalOffset += verticalStep;

  plotLabel modWTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod_w = 0");
  verticalOffset += verticalStep;

  FillDomain(&domain_, nullptr);
  for(const auto& pt : domain_) drawable_domain_.push_back(pt->at(0));

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

  int numberOfErrorCalculations = 1;
  QVector<int> additionalScreensSteps = {};

  /*
  for(int i = 990; i < 1011; ++i){
      additionalScreensSteps.append(i);
  }
  */

  double error_domain_length = 0;
  std::vector<std::vector<double>> error_domain = {};

  std::vector<double> model_values = {};
  std::vector<double> kde_values = {};

  ErrorsCalculator errors_calculator(
    &model_values, &kde_values, &error_domain, &error_domain_length
  );

  log("Crating CK Algorithm!");
  auto CKAlgorithm = EnhancedClusterKernelAlgorithm(number_of_cluster_kernels,
                                                    CreateNewVarianceBasedClusterKernel);

  for(step_number_ = 1; step_number_ < stepsNumber; ++step_number_) {
    clock_t executionStartTime = clock();

    Point stream_value = {};
    reader_->getNextRawDatum(&stream_value);
    UnivariateStreamElement element(stream_value);

    log("Performing step: " + QString::number(step_number_));
    CKAlgorithm.PerformStep(&element);
    log("Step performed.");

    target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

    if(step_number_ % screen_generation_frequency_ == 0 || step_number_ < 10
       || additionalScreensSteps.contains(step_number_)) {
      log("Drawing in step number " + QString::number(step_number_) + ".");

      // Error calculations
      if(step_number_ >= 1) {

        log("Getting error domain.");
        error_domain = CKAlgorithm.GetErrorDomain();

        log("Getting model plot on windowed.");
        model_values = GetFunctionsValueOnDomain(target_function_.get(), error_domain);
        log("Getting KDE plot on windowed.");
        kde_values = CKAlgorithm.GetKDEValuesOnDomain(error_domain);

        log("Getting model plot.");
        model_values = GetFunctionsValueOnDomain(target_function_.get(), error_domain);

        log("Calculating domain length.");

        error_domain_length =
            error_domain[error_domain.size() - 1][0] - error_domain[0][0];

        log("Calculating errors.");
        l1_w_ += errors_calculator.CalculateL1Error();
        l2_w_ += errors_calculator.CalculateL2Error();
        sup_w_ += errors_calculator.CalculateSupError();
        mod_w_ += errors_calculator.CalculateModError();

        ++numberOfErrorCalculations;
        log("Errors calculated.");
      }

      // ============ SUMS =========== //

      L1WTextLabel
          .setText("L1_w  =" + FormatNumberForDisplay(
              l1_w_ / numberOfErrorCalculations));
      L2WTextLabel
          .setText("L2_w  =" + FormatNumberForDisplay(
              l2_w_ / numberOfErrorCalculations));
      supWTextLabel
          .setText("sup_w =" + FormatNumberForDisplay(
              sup_w_ / numberOfErrorCalculations));
      modWTextLabel
          .setText("mod_w =" + FormatNumberForDisplay(
              mod_w_ / numberOfErrorCalculations));

      DrawPlots(&CKAlgorithm);

      for(const auto &label : plotLabels) label->updateText();

      ui->widget_plot->replot();
      QCoreApplication::processEvents();

      if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

      imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image saved: " + QString(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
    }
  }

  log("Animation finished.");
}



