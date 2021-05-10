#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <cfloat>
#include <QDebug>
#include <algorithm>
#include <chrono>
#include <QDateTime>
#include <Benchmarking/errorsCalculator.h>
#include <UI/plotLabelDoubleDataPreparator.h>

#include "kerDepCcWde.h"
#include "kerDepWindowedWde.h"

#include "UI/plotLabel.h"
#include "UI/plotLabelIntDataPreparator.h"

#include "Functions/complexfunction.h"

#include "Reservoir_sampling/biasedReservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/basicReservoirSamplingAlgorithm.h"

#include "Reservoir_sampling/distributiondataparser.h"
#include "Reservoir_sampling/progressivedistributiondatareader.h"
#include "Reservoir_sampling/textDataReader.h"

#include "ClusterKernelWrappers/varianceBasedClusterKernel.h"
#include "ClusterKernelWrappers/univariateStreamElement.h"

#include "LinearWDE.h"
#include "WeightedLinearWDE.h"
#include "WeightedThresholdedWDE.h"
#include "ThresholdingStrategies/hardThresholdingStrategy.h"
#include "ThresholdingStrategies/softThresholdingStrategy.h"

#include "SOMKEWrappers/somkeNormalKernel.h"
#include "SOMKEWrappers/MergingStrategies/somkeFixedMemoryMergingStrategy.h"
#include "SOMKEWrappers/MergingStrategies/somkeFixedThresholdMergingStrategy.h"

#include <QDateTime>
#include <QDate>
#include <QTime>

#include "UI/QwtContourPlotUI.h"

ClusterKernel *CreateNewVarianceBasedClusterKernel(ClusterKernelStreamElement *stream_element) {
  auto newClusterKernel = new VarianceBasedClusterKernel(stream_element);
  return newClusterKernel;
}

WaveletDensityEstimator *CreateWaveletDensityEstimatorFromBlock(const vector<double> &values_block) {
  auto wde = new LinearWDE();
  wde->UpdateWDEData(values_block);
  return wde;
}

WaveletDensityEstimator *CreateWeightedWaveletDensityEstimatorFromBlock(const vector<double> &values_block) {
  auto wde = new WeightedLinearWDE();
  wde->UpdateWDEData(values_block);
  return wde;
}

WaveletDensityEstimator *CreateWeightedThresholdedWaveletDensityEstimatorFromBlock(const vector<double> &values_block) {
  auto thresholding_strategy = ThresholdingStrategyPtr(new SoftThresholdingStrategy);
  auto wde = new WeightedThresholdedWDE(thresholding_strategy);
  wde->UpdateWDEData(values_block);
  return wde;
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

void log(const QString &msg) {
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
  for(auto value : drawable_domain_) {
    drawable_domain.push_back({value});
  }

  // Generate plot of model function
  if(ui->checkBox_showEstimatedPlot->isChecked()) {
    auto target_function_values = GetFunctionsValueOnDomain(target_function_.get(), drawable_domain);
    auto modelDistributionY = QVector<qreal>(target_function_values.begin(), target_function_values.end());
    AddPlot(&modelDistributionY, model_plot_pen_);
  }

  // Generate less elements KDE plot (navy blue)
  if(ui->checkBox_showEstimationPlot->isChecked()) {
    auto less_elements_estimator_values = DESDAAlgorithm->getKDEValues(&drawable_domain);
    auto less_elements_estimator_y = QVector<double>(less_elements_estimator_values.begin(),
                                                     less_elements_estimator_values.end());
    AddPlot(&less_elements_estimator_y, kde_plot_pen_);
  }

  // Generate weighted estimator plot (light blue)
  if(ui->checkBox_showWeightedEstimationPlot->isChecked()) {
    auto weighted_estimator_values = DESDAAlgorithm->getWeightedKDEValues(&drawable_domain);
    auto weighted_estimator_y = QVector<double>(weighted_estimator_values.begin(),
                                                weighted_estimator_values.end());
    AddPlot(&weighted_estimator_y, weighted_plot_pen_);
  }

  // Generate full estimator plot (BLACK)
  if(ui->checbox_showFullEstimator->isChecked()) {
    auto windowed_estimator_values = DESDAAlgorithm->getWindowKDEValues(&drawable_domain);
    auto windowed_estimator_y = QVector<double>(windowed_estimator_values.begin(),
                                                windowed_estimator_values.end());
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
    auto sigmoidally_enhanced_plot_values = DESDAAlgorithm->getEnhancedKDEValues(&drawable_domain);
    auto sigmoidally_enhanced_plot_y = QVector<double>(sigmoidally_enhanced_plot_values.begin(),
                                                       sigmoidally_enhanced_plot_values.end());
    AddPlot(&sigmoidally_enhanced_plot_y, desda_kde_plot_pen_);
  }

  if(ui->checkBox_showUnusualClusters->isChecked()) {
    atypical_elements_values_and_derivatives_ =
        DESDAAlgorithm->getAtypicalElementsValuesAndDerivatives();
    quantile_estimator_value_ = DESDAAlgorithm->_quantileEstimator;
    MarkUncommonClusters();
  }

  if(ui->checkBox_REESEKDE->isChecked()) {
    auto rare_elements_enhanced_plot_values = DESDAAlgorithm->getRareElementsEnhancedKDEValues(&drawable_domain);
    auto rare_elements_enhanced_plot_y = QVector<double>(rare_elements_enhanced_plot_values.begin(),
                                                         rare_elements_enhanced_plot_values.end());
    AddPlot(&rare_elements_enhanced_plot_y, desda_rare_elements_kde_plot_pen_);
  }
  // Draw plots
  ui->widget_plot->replot();
}

void MainWindow::DrawPlots(EnhancedClusterKernelAlgorithm *CKAlgorithm) {
  ClearPlot();
  ResizePlot();

  std::vector<std::vector<double>> drawable_domain = {}; // This is required for types :P
  for(auto value : drawable_domain_) {
    drawable_domain.push_back({value});
  }

  // Generate plot of model function
  if(ui->checkBox_showEstimatedPlot->isChecked()) {
    auto model_distribution_values = GetFunctionsValueOnDomain(target_function_.get(), drawable_domain);
    QVector<qreal> modelDistributionY = QVector<qreal>(model_distribution_values.begin(),
                                                       model_distribution_values.end());
    AddPlot(&modelDistributionY, model_plot_pen_);
  }

  // Generate less elements KDE plot (navy blue)
  if(ui->checkBox_showEstimationPlot->isChecked()) {
    auto less_elements_estimator_values = CKAlgorithm->GetKDEValuesOnDomain(drawable_domain);
    auto less_elements_estimator_y = QVector<double>(less_elements_estimator_values.begin(),
                                                     less_elements_estimator_values.end());
    AddPlot(&less_elements_estimator_y, kde_plot_pen_);
  }
}

void MainWindow::DrawPlots(KerDEP_CC_WDE *WDEAlgorithm) {
  ClearPlot();
  ResizePlot();

  std::vector<std::vector<double>> drawable_domain = {}; // This is required for types :P
  for(auto value : drawable_domain_) {
    drawable_domain.push_back({value});
  }

  // Generate plot of model function
  if(ui->checkBox_showEstimatedPlot->isChecked()) {
    auto model_distribution_values = GetFunctionsValueOnDomain(target_function_.get(), drawable_domain);
    QVector<qreal> modelDistributionY = QVector<qreal>(model_distribution_values.begin(),
                                                       model_distribution_values.end());
    AddPlot(&modelDistributionY, model_plot_pen_);
  }

  // Generate less elements KDE plot (navy blue)
  if(ui->checkBox_showEstimationPlot->isChecked()) {
    auto estimator_values = WDEAlgorithm->GetEstimatorValuesOnDomain(drawable_domain);
    auto estimator_y = QVector<double>(estimator_values.begin(), estimator_values.end());
    AddPlot(&estimator_y, kde_plot_pen_);
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

  while(i < maxX + 1) {
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
                                                  ->cellWidget(functionIndex,
                                                               static_cast<int>(TargetFunctionSettingsColumns::kStandardDeviationColumnIndex))))
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
                                                  ->cellWidget(functionIndex,
                                                               static_cast<int>(TargetFunctionSettingsColumns::kMeanColumnIndex))))
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

      for(auto dimensionVal : *(pPoint)) {
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
                                                      ->cellWidget(functionIndex,
                                                                   static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
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
        (dynamic_cast<QComboBox *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber,
                                                                                static_cast<int>(KernelSettingsColumns::kKernelColumnIndex))))
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
                                 ->cellWidget(targetFunctionElementsNumber - 1,
                                              static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))
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
                                                      ->cellWidget(functionIndex,
                                                                   static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
                             ->text().toDouble()
                     );

    elementalFunctions.push_back(
        std::shared_ptr<function>(new multivariateNormalProbabilityDensityFunction(means->at(functionIndex).get(),
                                                                                   stDevs->at(functionIndex).get())));
  }

  return new complexFunction(&contributions, &elementalFunctions);
}

int MainWindow::CanAnimationBePerformed(int dimensionsNumber) {
  if(dimensionsNumber == 1) {
    return 1;
  }
  else {
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

  auto smoothingParameterValidator =
      new QDoubleValidator(kMinSmoothingParameter, kMaxSmoothingParameter, kDecimalNumbers, this);
  smoothingParameterValidator->setLocale(locale);
  smoothingParameterValidator->setNotation(QDoubleValidator::StandardNotation);

  for(int rowNumber = 0; rowNumber < newNumberOfRows; ++rowNumber)
    AddKernelToTable(rowNumber, smoothingParameterValidator);
}

void MainWindow::AddKernelToTable(int rowNumber,
                                  QDoubleValidator *smoothingParameterValidator) {
  // Add combobox with kernels
  ui->tableWidget_dimensionKernels
    ->setCellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kKernelColumnIndex), new QComboBox());

  (dynamic_cast<QComboBox *>(ui->tableWidget_dimensionKernels
                               ->cellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kKernelColumnIndex))))
      ->insertItems(0, kernel_types_);

  // Add input box with validator for smoothing parameters
  ui->tableWidget_dimensionKernels
    ->setCellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kSmoothingParameterColumnIndex),
                    new QLineEdit());

  (dynamic_cast<QLineEdit *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber,
                                                                          static_cast<int>(KernelSettingsColumns::kSmoothingParameterColumnIndex))))
      ->setText("1.0");
  (dynamic_cast<QLineEdit *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber,
                                                                          static_cast<int>(KernelSettingsColumns::kSmoothingParameterColumnIndex))))
      ->setValidator(smoothingParameterValidator);

  // Add input box for carrier restriction value
  ui->tableWidget_dimensionKernels
    ->setCellWidget(rowNumber, static_cast<int>(KernelSettingsColumns::kCarrierRestrictionColumnIndex),
                    new QLineEdit());

  (dynamic_cast<QLineEdit *>(ui->tableWidget_dimensionKernels->cellWidget(rowNumber,
                                                                          static_cast<int>(KernelSettingsColumns::kCarrierRestrictionColumnIndex))))
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
    targetFunctionTablePointer
        ->setCellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kMeanColumnIndex),
                        new QTableWidget());

    meansTablePointer =
        dynamic_cast<QTableWidget *>(ui->tableWidget_targetFunctions->cellWidget(rowIndex,
                                                                                 static_cast<int>(TargetFunctionSettingsColumns::kMeanColumnIndex)));
    meansTablePointer->setRowCount(dimensionsNumber);
    meansTablePointer->setColumnCount(1);
    meansTablePointer->horizontalHeader()->hide();

    // TODO TR: Ensure that this doesn't result in memory leaks
    targetFunctionTablePointer
        ->setCellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kStandardDeviationColumnIndex),
                        new QTableWidget());

    stDevsTablePointer =
        dynamic_cast<QTableWidget *>(ui->tableWidget_targetFunctions->cellWidget(rowIndex,
                                                                                 static_cast<int>(TargetFunctionSettingsColumns::kStandardDeviationColumnIndex)));
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
    targetFunctionTablePointer
        ->setCellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex),
                        new QLineEdit());
    (dynamic_cast<QLineEdit *>(targetFunctionTablePointer
        ->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
        ->setMaxLength(6);
    (dynamic_cast<QLineEdit *>(targetFunctionTablePointer
        ->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
        ->setValidator(contributionValidator);
    QObject::connect(
        (dynamic_cast<QLineEdit *>(targetFunctionTablePointer
            ->cellWidget(rowIndex, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex)))),
        SIGNAL(textEdited(QString)), this, SLOT(UpdateLastContribution()));
  }

  // Disable last contribution cell, as it's filled automatically
  (dynamic_cast<QLineEdit *>(targetFunctionTablePointer
      ->cellWidget(numberOfRows - 1, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
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
  Run1DExperimentWithDESDA();
  //Run1DExperimentWithClusterKernels();
  //Run1DExperimentWithWDE();
  //Run1DExperimentWithSOMKE();
}

void MainWindow::on_pushButton_clicked() {
  log("2D Experiment start.");

  screen_generation_frequency_ = 10;
  int seed = ui->lineEdit_seed->text().toInt();
  int m0 = ui->lineEdit_sampleSize->text().toInt();

  // Contour levels calculation.
  QList<double> contourLevels;
  double level = 0.025;
  while(level < 0.21) {
    contourLevels += level;
    level += 0.025;
  }

  // Add clusters_ to the estimator
  means_ = {std::make_shared<std::vector<double >>()};
  means_.back()->push_back(0);
  means_.back()->push_back(0);

  standard_deviations_ = {std::make_shared<std::vector<double >>()};
  standard_deviations_.back()->push_back(1);
  standard_deviations_.back()->push_back(1);

  //auto densityFunction =
  //    new multivariateNormalProbabilityDensityFunction(means_.back().get(), standard_deviations_.back().get());
  //contour_plot_->addQwtPlotSpectrogram(new SpectrogramData2(densityFunction, -10.0), QPen(QColor(255, 0, 0)));

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

  //std::shared_ptr<distribution>
  //    targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));
  std::vector<double> meansForDistribution = {0.0, 0.0};
  std::vector<double> stDevsForDistribution = {1.0, 1.0};

  parser_.reset(new distributionDataParser(&attributes_data_));

  std::string data_path = "k:\\Coding\\Python\\KerDEP_Data_Preparator\\MetroInterstateTraffic\\result_2D.txt";
  //data_path = "y:\\Data\\metro2017_2D.txt";

  // reader_.reset(new TextDataReader("k:\\Coding\\Python\\KerDEP_Data_Preparator\\result.txt"));
  // reader_.reset(new TextDataReader("k:\\Coding\\Python\\KerDEP_Data_Preparator\\AirQuality\\result.txt"));
  reader_.reset(new TextDataReader(data_path));
  //reader_.reset(new TextDataReader("k:\\Coding\\Python\\KerDEP_Data_Preparator\\Cracow_Temp_2016\\result.txt"));

  /*
  reader_.reset(
      new progressiveDistributionDataReader(targetDistribution.get(), 0,
                                            0,
                                            new normalDistribution(0, &meansForDistribution, &stDevsForDistribution,
                                                                   55))
               );

  */

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
      ui->lineEdit_rarity->text().toDouble(), newWeightB, pluginRank
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

  QDate data_start_date(2016, 10, 1); // Metro
  QTime data_start_time(0, 0, 0);
  QDateTime data_date_time(data_start_date, data_start_time);

  QwtContourPlotUI plotUi(&step_number_, screen_generation_frequency_, seed,
                          &DESDAAlgorithm, &l1_n_, &l2_n_, &sup_n_, &mod_n_,
                          &actual_l1, &actual_l2, &actual_sup, &actual_mod,
                          &data_date_time);
  plotUi.attach(contour_plot_);
  plotUi.updateTexts();
  plotUi.SetErrorsPrinting(false);
  //QVector<int> initialDrawingSteps = {1, 2, 3, 4, 5, 6, 7, 8, 9 , 10};
  QVector<int> initialDrawingSteps = {1};
  std::vector<double> model_function_values = {};
  std::vector<double> estimator_values = {};
  double domain_area = 0;
  std::vector<std::vector<double>> error_domain = {};
  ErrorsCalculator errors_calculator(
      &model_function_values, &estimator_values, &error_domain, &domain_area
                                    );

  // Prepare image location.
  QString expNum = "1527-1 (2D)";
  this->setWindowTitle("Experiment #" + expNum);
  QString expDesc =
      "iw=" + QString::number(screen_generation_frequency_)
      + ", v=0, seed = " + QString::number(seed) +
      ", m0=" + QString::number(m0) + ", mMin=" + QString::number(DESDAAlgorithm._minM) + ", sz001";
  //QString driveDir = "\\\\beabourg\\private\\"; // WIT PCs
  //QString driveDir = "Y:\\"; // WIT PCs after update
  QString driveDir = "D:\\Test\\"; // Home
  //QString driveDir = "d:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\";
  //QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " (" + expDesc + ")\\";
  QString dirPath = driveDir + "Eksperyment " + expNum + " (" + expDesc + ")\\";
  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  log("Experiment started.");
  for(step_number_ = 1; step_number_ < 42000; ++step_number_) {

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

      /*
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
       */

      // densityFunction->setMeans(*means_.back().get());

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
      log("Saved: " + QString::number(ui->widget_contour_plot_holder->grab().save(imageName)));
      log("Drawing finished.");
    }

    if(step_number_ == 10) {
      initialDrawingSteps.clear();  // To reduce comparisons for drawing.
    }

    endTime = time(nullptr);

    data_date_time = data_date_time.addSecs(3600); // Add hour to the date

    log("Step time: " + QString::number(endTime - startTime) + " s");
  }

  log("Done!");
}

void MainWindow::resizeEvent(QResizeEvent *event) {
  int offset = 10; // Offset in px, so that scale is in
  QMainWindow::resizeEvent(event);
  int newSize = std::min(ui->widget_contour_plot->height(),
                         ui->widget_contour_plot->width()) - offset;
  //contour_plot_->resize(newSize, newSize);
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

std::vector<double> MainWindow::GetFunctionsValueOnDomain(function *func,
                                                          const std::vector<std::vector<double>> &domain) {
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

  /*
  reader_.reset(
      new progressiveDistributionDataReader(targetDistribution.get(),
                                            progressionSize,
                                            0,  // Delay
                                            new normalDistribution(seedString.toInt(), &alternativeDistributionMean,
                                                                   &alternativeDistributionStDevs, 55))
               );
  */

  //std::string data_path = "k:\\Coding\\Python\\KerDEP_Data_Preparator\\MetroInterstateTraffic\\result.txt";
  // std::string data_path = "k:\\Coding\\Python\\KerDEP_Data_Preparator\\AirQuality\\result.txt";
  // std::string data_path = "k:\\Coding\\Python\\KerDEP_Data_Preparator\\Cracow_Temp_2016\\result.txt";
  std::string data_path = "k:\\Coding\\Python\\KerDEP_Data_Preparator\\BikeSharing\\result.txt";
  data_path = "y:\\Data\\BikeSharingPrices.txt";

  // reader_.reset(new TextDataReader("k:\\Coding\\Python\\KerDEP_Data_Preparator\\AirQuality\\result.txt"));
  reader_.reset(new TextDataReader(data_path));
  //reader_.reset(new TextDataReader("k:\\Coding\\Python\\KerDEP_Data_Preparator\\Cracow_Temp_2016\\result.txt"));

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  reservoirSamplingAlgorithm *algorithm =
      GenerateReservoirSamplingAlgorithm(reader_.get(), parser_.get());

  objects_.clear();

  int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();

  log("Attributes data set.");

  int sampleSize = ui->lineEdit_sampleSize->text().toInt();

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
      ui->lineEdit_rarity->text().toDouble(), newWeightB, pluginRank
                      );

  QString expNum = "1523";
  this->setWindowTitle("Experiment #" + expNum);
  QString expDesc = "DESDA, Plugin" + QString::number(pluginRank) +
                    ", Cracow, m0=" + QString::number(DESDAAlgorithm._maxM) +
                    ", mMin=" + QString::number(DESDAAlgorithm._minM) +
                    ", home";

  screen_generation_frequency_ = 10;
  bool compute_errors = false;

  //QString driveDir = "D:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\"; // Home
  QString driveDir = "D:\\Test\\"; // Test
  //QString driveDir = "Y:\\"; // WIT PCs after update

  //QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " (" + expDesc + ")\\";
  //QString dirPath = driveDir + "Badania PK\\Eksperyment " + expNum + " (" + expDesc + ")\\";
  QString dirPath = driveDir + "Eksperyment " + expNum + " (" + expDesc + ")\\";

  ClearPlot();
  ResizePlot();

  // Initial screen should only contain exp number (as requested).
  plotLabel expNumLabel(ui->widget_plot, 0.02, 0.25, "Exp." + expNum);
  expNumLabel.setFont(QFont("Courier New", 250));

  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  QString imageName = dirPath + QString::number(0) + ".png";

  log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
  expNumLabel.setText("");

  plot_labels_ = {};
  label_horizontal_offset_ = label_vertical_offset_ = 0.01;

  // Exps with days
    // Bike Sharing Experiment
  QDate startDate(2011, 1, 1);
  QTime startTime(0, 0, 0);
    // Air Quality Italy Experiment
  //QDate startDate(2004, 3, 10);
  //QTime startTime(18, 0, 0);
   // Metro Minneapolis Experiment
  //QDate startDate(2016, 10, 1);
  //QTime startTime(0, 0, 0);

  QDateTime dateTime(startDate, startTime);

  plotLabel date_label(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_, "");
  label_vertical_offset_ += label_vertical_offset_step_;
  // END Exps with days

  AddIntLabelToPlot("i     = ", &step_number_);
  AddConstantLabelToPlot("iw    = " + QString::number(screen_generation_frequency_));
  AddConstantLabelToPlot("seed  = " + seedString);
  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabel KPSSTextLabel(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_,
                          "KPSS     = 0");
  label_vertical_offset_ += label_vertical_offset_step_;

  AddDoubleLabelToPlot("sgmKPSS  = ", &DESDAAlgorithm._sgmKPSS);
  label_vertical_offset_ += label_vertical_offset_step_;

  AddConstantLabelToPlot("mKPSS = " + QString::number(DESDAAlgorithm._kpssM));
  AddConstantLabelToPlot("m0    = " + ui->lineEdit_sampleSize->text());
  AddConstantLabelToPlot("mmin  = " + QString::number(DESDAAlgorithm._minM));
  AddIntLabelToPlot("m     = ", &(DESDAAlgorithm._m));
  label_vertical_offset_ += label_vertical_offset_step_;

  AddDoubleLabelToPlot("beta0 = ", &(DESDAAlgorithm._beta0));
  label_vertical_offset_ += label_vertical_offset_step_;

  AddDoubleLabelToPlot("r     = ", &(DESDAAlgorithm._r));
  AddDoubleLabelToPlot("q     = ", &(DESDAAlgorithm._quantileEstimator));
  AddIntLabelToPlot("rare  = ", &(DESDAAlgorithm._rareElementsNumber));
  AddIntLabelToPlot("trend = ", &(DESDAAlgorithm._trendsNumber));
  label_vertical_offset_ += 5 * label_vertical_offset_step_;

  AddColorsLegendToPlot();

  //====================  SECOND COLUMN =================//

  label_horizontal_offset_ = 0.20;
  label_vertical_offset_ = 0.01 + 9 * label_vertical_offset_step_;

  //==================== ERRORS SUM =================//

  label_horizontal_offset_ = 0.87;
  label_vertical_offset_ = 0.01;

  QVector<double> l1_errors_sums = {};
  QVector<double> l2_errors_sums = {};
  QVector<double> sup_errors_sums = {};
  QVector<double> mod_errors_sums = {};

  QVector<double_ptr> l1_errors = {};
  QVector<double_ptr> l2_errors = {};
  QVector<double_ptr> sup_errors = {};
  QVector<double_ptr> mod_errors = {};

  if(compute_errors) {

    QVector<QString> l1_labels = {"L1_w  = ", "L1_m  = ", "L1_d  = ", "L1_p  = ", "L1_n  = "};
    QVector<QString> l2_labels = {"L2_w  = ", "L2_m  = ", "L2_d  = ", "L2_p  = ", "L2_n  = "};
    QVector<QString> sup_labels = {"sup_w = ", "sup_m = ", "sup_d = ", "sup_p = ", "sup_n = "};
    QVector<QString> mod_labels = {"mod_w = ", "mod_m = ", "mod_d = ", "mod_p = ", "mod_n = "};

    for(size_t i = 0; i < l1_labels.size(); ++i){
      l1_errors.push_back(std::make_shared<double>(0));
      l1_errors_sums.push_back(0);
      l2_errors.push_back(std::make_shared<double>(0));
      l2_errors_sums.push_back(0);
      sup_errors.push_back(std::make_shared<double>(0));
      sup_errors_sums.push_back(0);
      mod_errors.push_back(std::make_shared<double>(0));
      mod_errors_sums.push_back(0);
    }

    AddDoubleLabelsToPlot(l1_labels, l1_errors);
    label_vertical_offset_ += label_vertical_offset_step_;

    AddDoubleLabelsToPlot(l2_labels, l2_errors);
    label_vertical_offset_ += label_vertical_offset_step_;

    AddDoubleLabelsToPlot(sup_labels, sup_errors);
    label_vertical_offset_ += label_vertical_offset_step_;

    AddDoubleLabelsToPlot(mod_labels, mod_errors);
    label_vertical_offset_ += label_vertical_offset_step_;
  }

  FillDomain(&domain_, nullptr);
  for(const auto &pt : domain_){ drawable_domain_.push_back(pt->at(0)); }

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

  int numberOfErrorCalculations = 1;
  QVector<int> additionalScreensSteps = {};

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

  QVector<ErrorsCalculator*> errors_calculators = {
      &windowed_errors_calculator, &less_elements_kde_errors_calculator, &weighted_kde_errors_calculator,
      &enhanced_kde_errors_calculator, &rare_elements_kde_errors_calculator
  };


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
      if(step_number_ >= 1 && compute_errors) {

        log("Getting windowed domain.");
        windowed_error_domain = Generate1DWindowedPlotErrorDomain(&DESDAAlgorithm);
        log("Getting non-windowed domain.");
        error_domain = Generate1DPlotErrorDomain(&DESDAAlgorithm);

        log("Getting model plot on windowed.");
        windowed_model_values = GetFunctionsValueOnDomain(target_function_.get(), windowed_error_domain);
        log("Getting KDE plot on windowed.");
        windowed_kde_values = DESDAAlgorithm.getWindowKDEValues(&windowed_error_domain);

        log("Getting model plot.");
        model_values = GetFunctionsValueOnDomain(target_function_.get(), error_domain);
        log("Getting KDE plot on lesser elements.");
        less_elements_kde_values = DESDAAlgorithm.getKDEValues(&error_domain);
        log("Getting weighted KDE plot.");
        weighted_kde_values = DESDAAlgorithm.getWeightedKDEValues(&error_domain);
        log("Getting sgm KDE plot.");
        enhanced_kde_values = DESDAAlgorithm.getEnhancedKDEValues(&error_domain);
        log("Getting rare KDE plot.");
        rare_elements_kde_values = DESDAAlgorithm.getRareElementsEnhancedKDEValues(&error_domain);

        error_domain_length = error_domain[error_domain.size() - 1][0] - error_domain[0][0];
        windowed_error_domain_length =
            windowed_error_domain[windowed_error_domain.size() - 1][0] - windowed_error_domain[0][0];

        AddL1ErrorsToSum(errors_calculators, l1_errors_sums);
        AddL2ErrorsToSum(errors_calculators, l2_errors_sums);
        AddSupErrorsToSum(errors_calculators, sup_errors_sums);
        AddModErrorsToSum(errors_calculators, mod_errors_sums);

        for(size_t i = 0; i < errors_calculators.size(); ++i){
          *l1_errors[i] = l1_errors_sums[i] / numberOfErrorCalculations;
          *l2_errors[i] = l2_errors_sums[i] / numberOfErrorCalculations;
          *sup_errors[i] = sup_errors_sums[i] / numberOfErrorCalculations;
          *mod_errors[i] = mod_errors_sums[i] / numberOfErrorCalculations;
        }

        ++numberOfErrorCalculations;
      }

      // ============= LEFT SIDE UPDATE ================ //

      KPSSTextLabel.setText("KPSS     = " + FormatNumberForDisplay(
          DESDAAlgorithm.getStationarityTestValue()));

      DrawPlots(&DESDAAlgorithm);

      for(const auto &label : plot_labels_){
        label->updateText();
      }

      date_label.setText(dateTime.toString());

      ui->widget_plot->replot();
      QCoreApplication::processEvents();

      if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

      imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
    }

    dateTime = dateTime.addSecs(3600); // Bike sharing
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

  int number_of_cluster_kernels = 100;
  step_number_ = 0;
  double h = 1;
  double sigma = 0;

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
  QString expNum = "1386 (CK)";
  this->setWindowTitle("Experiment #" + expNum);
  QString expDesc = "v=0, m = " + QString::number(number_of_cluster_kernels)
                    + ", mean-var-resampling, weighted list-based algorithm, weighted StDev, updated h coefficient, alpha=0.01";
  screen_generation_frequency_ = 10;

  //QString driveDir = "\\\\beabourg\\private\\"; // WIT PCs
  //QString driveDir = "D:\\Test\\"; // Home
  QString driveDir = "Y:\\"; // WIT PCs after update
  //QString driveDir = "d:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\";
  QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " ("
                    + expDesc + ")\\";

  ClearPlot();
  ResizePlot();

  // Initial screen should only contain exp number (as requested).
  plotLabel expNumLabel(ui->widget_plot, 0.02, 0.25, "Exp." + expNum);
  expNumLabel.setFont(QFont("Courier New", 250));

  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  QString imageName = dirPath + QString::number(0) + ".png";

  log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
  expNumLabel.setText("");

  // Setting up the labels
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

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "m     = ",
                                                   &(number_of_cluster_kernels),
                                                   std::make_shared<plotLabelIntDataPreparator>()));

  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "h     = ", &(h),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));
  verticalOffset += verticalStep;
  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "sigma = ", &(sigma),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));



  //====================  SECOND COLUMN =================//

  horizontalOffset = 0.20;
  verticalOffset = 0.01 + 9 * verticalStep;

  //====================== ERRORS SUM ===================//

  horizontalOffset = 0.87;
  verticalOffset = 0.01;

  plotLabel L1TextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1   = 0");
  verticalOffset += verticalStep;
  plotLabel L1aTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1a  = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel L2TextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2   = 0");
  verticalOffset += verticalStep;
  plotLabel L2aTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2a  = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel supTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup  = 0");
  verticalOffset += verticalStep;
  plotLabel supaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "supa = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel modTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod  = 0");
  verticalOffset += verticalStep;
  plotLabel modaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "moda = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  FillDomain(&domain_, nullptr);
  for(const auto &pt : domain_) drawable_domain_.push_back(pt->at(0));

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

  int numberOfErrorCalculations = 0;
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
  double l1_sum = 0;
  double l2_sum = 0;
  double sup_sum = 0;
  double mod_sum = 0;

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
        l1_w_ = errors_calculator.CalculateL1Error();
        l2_w_ = errors_calculator.CalculateL2Error();
        sup_w_ = errors_calculator.CalculateSupError();
        mod_w_ = errors_calculator.CalculateModError();
        l1_sum += l1_w_;
        l2_sum += l2_w_;
        sup_sum += sup_w_;
        mod_sum += mod_w_;

        ++numberOfErrorCalculations;
        log("Errors calculated.");
      }

      // ============ SUMS =========== //

      L1TextLabel
          .setText("L1   =" + FormatNumberForDisplay(
              l1_sum / numberOfErrorCalculations));
      L2TextLabel
          .setText("L2   =" + FormatNumberForDisplay(
              l2_sum / numberOfErrorCalculations));
      supTextLabel
          .setText("sup  =" + FormatNumberForDisplay(
              sup_sum / numberOfErrorCalculations));
      modTextLabel
          .setText("mod  =" + FormatNumberForDisplay(
              mod_sum / numberOfErrorCalculations));

      L1aTextLabel
          .setText("L1a  =" + FormatNumberForDisplay(l1_w_));
      L2aTextLabel
          .setText("L2a  =" + FormatNumberForDisplay(l2_w_));
      supaTextLabel
          .setText("supa =" + FormatNumberForDisplay(sup_w_));
      modaTextLabel
          .setText("moda =" + FormatNumberForDisplay(mod_w_));

      DrawPlots(&CKAlgorithm);

      h = CKAlgorithm.GetBandwidth();
      sigma = CKAlgorithm.GetStandardDeviation();

      for(const auto &label : plotLabels) label->updateText();

      ui->widget_plot->replot();
      QCoreApplication::processEvents();

      if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

      imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
    }
  }

  log("Animation finished.");
}

void MainWindow::Run1DExperimentWithWDE() {
  // TR TODO: This is basically the same as it is in Cluster Kernels... Initialization should
  // be made an separate function.
  int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

  if(!CanAnimationBePerformed(dimensionsNumber)) return;

  QString seedString = ui->lineEdit_seed->text();

  // Log that application started generating KDE
  // Standard seed was 5625.
  log("KDE animation with WDE started.");
  log("Seed: " + seedString);
  log("Sample size: " + ui->lineEdit_sampleSize->text());

  step_number_ = 0;
  double sigma = 0;

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
  double weight_modifier = 0.95; // omega
  unsigned int maximal_number_of_coefficients = 100; // M
  unsigned int current_coefficients_number = 0; // #coef
  int number_of_elements_per_block = 1000; // b

  QString expNum = "1491 TEST (Thresholded Weighted Window WDE)";
  //QString expNum = "THRESHOLDED_WDE_TEST_1";
  this->setWindowTitle("Experiment #" + expNum);
  QString expDesc = "v=tor klasyczny, soft threshold, b=" + QString::number(number_of_elements_per_block) +
                    ", omega=" + QString::number(weight_modifier) +
                    ", M=" + QString::number(maximal_number_of_coefficients);
  screen_generation_frequency_ = 10;

  //QString driveDir = "\\\\beabourg\\private\\"; // WIT PCs
  //QString driveDir = "D:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\Doktorat\\"; // Home
  //QString driveDir = "Y:\\"; // WIT PCs after update
  QString driveDir = "D:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\";
  QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " ("
                    + expDesc + ")\\";

  ClearPlot();
  ResizePlot();

  // Initial screen should only contain exp number (as requested).
  plotLabel expNumLabel(ui->widget_plot, 0.02, 0.25, "Exp." + expNum);
  expNumLabel.setFont(QFont("Courier New", 250));

  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  QString imageName = dirPath + QString::number(0) + ".png";

  log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
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

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "b     = ",
                                                   &(number_of_elements_per_block),
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "M     = ",
                                                   &(maximal_number_of_coefficients),
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "#coef = ",
                                                   &(current_coefficients_number),
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "omega = ", &(weight_modifier),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));


  //====================  SECOND COLUMN =================//

  horizontalOffset = 0.20;
  verticalOffset = 0.01 + 9 * verticalStep;

  //====================== ERRORS SUM ===================//

  horizontalOffset = 0.87;
  verticalOffset = 0.01;

  plotLabel L1TextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1   = 0");
  verticalOffset += verticalStep;
  plotLabel L1aTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1a  = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel L2TextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2   = 0");
  verticalOffset += verticalStep;
  plotLabel L2aTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2a  = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel supTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup  = 0");
  verticalOffset += verticalStep;
  plotLabel supaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "supa = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel modTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod  = 0");
  verticalOffset += verticalStep;
  plotLabel modaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "moda = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  FillDomain(&domain_, nullptr);
  for(const auto &pt : domain_) drawable_domain_.push_back(pt->at(0));

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

  int numberOfErrorCalculations = 0;
  QVector<int> additionalScreensSteps = {1};

  /*
  for(int i = 990; i < 1011; ++i){
      additionalScreensSteps.append(i);
  }
  */

  double error_domain_length = 0;
  std::vector<std::vector<double>> error_domain = {};
  std::vector<double> model_values = {};
  std::vector<double> wde_values = {};
  ErrorsCalculator errors_calculator(
      &model_values, &wde_values, &error_domain, &error_domain_length
                                    );
  Windowed_WDE WDE_Algorithm = Windowed_WDE(maximal_number_of_coefficients, weight_modifier,
                                            CreateWeightedThresholdedWaveletDensityEstimatorFromBlock,
      //CreateWeightedWaveletDensityEstimatorFromBlock,
                                            number_of_elements_per_block);
  double l1_sum = 0;
  double l2_sum = 0;
  double sup_sum = 0;
  double mod_sum = 0;

  for(step_number_ = 1; step_number_ < stepsNumber; ++step_number_) {
    clock_t executionStartTime = clock();
    Point stream_value = {};
    reader_->getNextRawDatum(&stream_value);

    log("Performing step: " + QString::number(step_number_));
    WDE_Algorithm.PerformStep(&stream_value);
    log("Step performed.");

    target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

    if(step_number_ % screen_generation_frequency_ == 0 || additionalScreensSteps.contains(step_number_)) {
      log("Drawing in step number " + QString::number(step_number_) + ".");

      // Error calculations
      if(step_number_ >= 1000) {

        log("Getting error domain.");
        error_domain = WDE_Algorithm.GetErrorDomain();

        log("Getting model plot on windowed.");
        model_values = GetFunctionsValueOnDomain(target_function_.get(), error_domain);
        log("Getting KDE plot on windowed.");
        wde_values = WDE_Algorithm.GetEstimatorValuesOnDomain(error_domain);

        log("Getting model plot.");
        model_values = GetFunctionsValueOnDomain(target_function_.get(), error_domain);

        log("Calculating domain length.");

        error_domain_length =
            error_domain[error_domain.size() - 1][0] - error_domain[0][0];

        log("Calculating errors.");
        l1_w_ = errors_calculator.CalculateL1Error();
        l2_w_ = errors_calculator.CalculateL2Error();
        sup_w_ = errors_calculator.CalculateSupError();
        mod_w_ = errors_calculator.CalculateModError();
        l1_sum += l1_w_;
        l2_sum += l2_w_;
        sup_sum += sup_w_;
        mod_sum += mod_w_;

        ++numberOfErrorCalculations;
        log("Errors calculated.");
      }

      // ============ SUMS =========== //

      L1TextLabel
          .setText("L1   =" + FormatNumberForDisplay(
              l1_sum / numberOfErrorCalculations));
      L2TextLabel
          .setText("L2   =" + FormatNumberForDisplay(
              l2_sum / numberOfErrorCalculations));
      supTextLabel
          .setText("sup  =" + FormatNumberForDisplay(
              sup_sum / numberOfErrorCalculations));
      modTextLabel
          .setText("mod  =" + FormatNumberForDisplay(
              mod_sum / numberOfErrorCalculations));

      L1aTextLabel
          .setText("L1a  =" + FormatNumberForDisplay(l1_w_));
      L2aTextLabel
          .setText("L2a  =" + FormatNumberForDisplay(l2_w_));
      supaTextLabel
          .setText("supa =" + FormatNumberForDisplay(sup_w_));
      modaTextLabel
          .setText("moda =" + FormatNumberForDisplay(mod_w_));

      current_coefficients_number = WDE_Algorithm.GetCurrentCoefficientsNumber();

      DrawPlots(&WDE_Algorithm);

      for(const auto &label : plotLabels) label->updateText();

      ui->widget_plot->replot();
      QCoreApplication::processEvents();

      if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

      imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
    }
  }

  log("Experiment finished!");
}

void MainWindow::Run1DExperimentWithSOMKE() {
  int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

  if(!CanAnimationBePerformed(dimensionsNumber)) return;

  QString seedString = ui->lineEdit_seed->text();

  // Log that application started generating KDE
  // Standard seed was 5625.
  log("KDE animation with WDE started.");
  log("Seed: " + seedString);
  log("Sample size: " + ui->lineEdit_sampleSize->text());

  step_number_ = 0;
  double sigma = 0;

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
  int neurons_number = 100;
  int epochs_number = 3000;
  int data_window_size = 500;
  int max_number_of_som_seq_entries = 1;
  double beta = 0;
  double alpha = 1.0;

  double sigma0 = 25.0;
  double tau1 = 1000 / log(sigma0);
  double tau2 = 1000.0;
  double eta0 = 3.0;

  QString expNum = "1497 (SOMKE)";
  this->setWindowTitle("Experiment #" + expNum);
  QString expDesc = "v=tor klasyczny, fixed threshold, original training"
                    ", max_entries=" + QString::number(max_number_of_som_seq_entries) +
                    ", neurons_num=" + QString::number(neurons_number) +
                    ", window_size=" + QString::number(data_window_size) +
                    ", epochs_num= " + QString::number(epochs_number);
  screen_generation_frequency_ = 10;

  //QString driveDir = "\\\\beabourg\\private\\"; // WIT PCs
  //QString driveDir = "D:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\Doktorat\\"; // Home
  //QString driveDir = "Y:\\"; // WIT PCs after update
  QString driveDir = "D:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\";
  QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " ("
                    + expDesc + ")\\";

  ClearPlot();
  ResizePlot();

  // Initial screen should only contain exp number (as requested).
  plotLabel expNumLabel(ui->widget_plot, 0.02, 0.25, "Exp." + expNum);
  expNumLabel.setFont(QFont("Courier New", 250));

  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  QString imageName = dirPath + QString::number(0) + ".png";

  log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
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

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "Neurons  = ", &(neurons_number),
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  /*
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "Seq_max  = ",
                                                   &(max_number_of_som_seq_entries),
                                                   std::make_shared<plotLabelIntDataPreparator>()));
 */

  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "win_size = ", &(data_window_size),
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "epochs   = ", &(epochs_number),
                                                   std::make_shared<plotLabelIntDataPreparator>()));

  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "beta     = ", &(beta),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "alpha    = ", &(alpha),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "sigma0   = ", &(sigma0),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "eta0     = ", &(eta0),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "tau1     = ", &(tau1),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  verticalOffset += verticalStep;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "tau2     = ", &(tau2),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));


  //====================  SECOND COLUMN =================//

  horizontalOffset = 0.20;
  verticalOffset = 0.01 + 9 * verticalStep;

  //====================== ERRORS SUM ===================//

  horizontalOffset = 0.87;
  verticalOffset = 0.01;

  plotLabel L1TextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1   = 0");
  verticalOffset += verticalStep;
  plotLabel L1aTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L1a  = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel L2TextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2   = 0");
  verticalOffset += verticalStep;
  plotLabel L2aTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2a  = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel supTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "sup  = 0");
  verticalOffset += verticalStep;
  plotLabel supaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "supa = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  plotLabel modTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "mod  = 0");
  verticalOffset += verticalStep;
  plotLabel modaTextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "moda = 0");
  verticalOffset += verticalStep;
  verticalOffset += verticalStep;

  FillDomain(&domain_, nullptr);
  for(const auto &pt : domain_) drawable_domain_.push_back(pt->at(0));

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

  int numberOfErrorCalculations = 0;
  QVector<int> additionalScreensSteps = {1};

  /*
  for(int i = 990; i < 1011; ++i){
      additionalScreensSteps.append(i);
  }
  */

  double error_domain_length = 0;
  std::vector<std::vector<double>> error_domain = {};
  std::vector<double> model_values = {};
  std::vector<double> somke_values = {};
  ErrorsCalculator errors_calculator(
      &model_values, &somke_values, &error_domain, &error_domain_length
                                    );
  KernelPtr kernel(new SOMKENormalKernel());
  // MergingStrategyPtr merging_strategy(new SOMKEFixedMemoryMergingStrategy(max_number_of_som_seq_entries, beta));
  MergingStrategyPtr merging_strategy(new SOMKEFixedThresholdMergingStrategy(alpha, beta));
  SOMKEAlgorithm somke_algorithm(kernel, merging_strategy, neurons_number, epochs_number, data_window_size);
  double l1_sum = 0;
  double l2_sum = 0;
  double sup_sum = 0;
  double mod_sum = 0;

  for(step_number_ = 1; step_number_ < stepsNumber; ++step_number_) {
    clock_t executionStartTime = clock();
    Point stream_value = {};
    reader_->getNextRawDatum(&stream_value);

    log("Performing step: " + QString::number(step_number_));
    somke_algorithm.PerformStep(stream_value);
    log("Step performed.");

    target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

    if(step_number_ % screen_generation_frequency_ == 0 || additionalScreensSteps.contains(step_number_)) {
      log("Drawing in step number " + QString::number(step_number_) + ".");

      // Error calculations
      if(step_number_ >= 1000) {

        log("Getting error domain.");
        error_domain = somke_algorithm.divergence_domain_;

        log("Getting model plot on windowed.");
        model_values = GetFunctionsValueOnDomain(target_function_.get(), error_domain);
        log("Getting KDE plot on windowed.");
        somke_values = {};
        for(auto pt : error_domain) {
          somke_values.push_back(somke_algorithm.GetValue(pt));
        }

        log("Getting model plot.");
        model_values = GetFunctionsValueOnDomain(target_function_.get(), error_domain);

        log("Calculating domain length.");

        error_domain_length =
            error_domain[error_domain.size() - 1][0] - error_domain[0][0];

        log("Calculating errors.");
        l1_w_ = errors_calculator.CalculateL1Error();
        l2_w_ = errors_calculator.CalculateL2Error();
        sup_w_ = errors_calculator.CalculateSupError();
        mod_w_ = errors_calculator.CalculateModError();
        l1_sum += l1_w_;
        l2_sum += l2_w_;
        sup_sum += sup_w_;
        mod_sum += mod_w_;

        ++numberOfErrorCalculations;
        log("Errors calculated.");
      }

      // ============ SUMS =========== //

      L1TextLabel
          .setText("L1   =" + FormatNumberForDisplay(
              l1_sum / numberOfErrorCalculations));
      L2TextLabel
          .setText("L2   =" + FormatNumberForDisplay(
              l2_sum / numberOfErrorCalculations));
      supTextLabel
          .setText("sup  =" + FormatNumberForDisplay(
              sup_sum / numberOfErrorCalculations));
      modTextLabel
          .setText("mod  =" + FormatNumberForDisplay(
              mod_sum / numberOfErrorCalculations));

      L1aTextLabel
          .setText("L1a  =" + FormatNumberForDisplay(l1_w_));
      L2aTextLabel
          .setText("L2a  =" + FormatNumberForDisplay(l2_w_));
      supaTextLabel
          .setText("supa =" + FormatNumberForDisplay(sup_w_));
      modaTextLabel
          .setText("moda =" + FormatNumberForDisplay(mod_w_));

      DrawPlots(&somke_algorithm);

      for(const auto &label : plotLabels) label->updateText();

      ui->widget_plot->replot();
      QCoreApplication::processEvents();

      if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

      imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
    }
  }

  log("Experiment finished!");
}

void MainWindow::DrawPlots(SOMKEAlgorithm *somke_algorithm) {
  ClearPlot();
  ResizePlot();

  std::vector<std::vector<double>> drawable_domain = {}; // This is required for types :P
  for(auto value : drawable_domain_) {
    drawable_domain.push_back({value});
  }

  // Generate plot of model function
  if(ui->checkBox_showEstimatedPlot->isChecked()) {
    auto model_distribution_values = GetFunctionsValueOnDomain(target_function_.get(), drawable_domain);
    QVector<qreal> modelDistributionY = QVector<qreal>(model_distribution_values.begin(),
                                                       model_distribution_values.end());
    AddPlot(&modelDistributionY, model_plot_pen_);
  }

  // Generate less elements KDE plot (navy blue)
  if(ui->checkBox_showEstimationPlot->isChecked()) {
    vector<double> estimator_values = {};
    for(auto pt : drawable_domain) {
      estimator_values.push_back(somke_algorithm->GetValue(pt));
    }
    auto estimator_y = QVector<double>(estimator_values.begin(), estimator_values.end());
    AddPlot(&estimator_y, kde_plot_pen_);
  }
}

void MainWindow::AddDoubleLabelsToPlot(const QVector<QString> &labels, const QVector<double_ptr> &values) {
  // TODO TR: I assume that labels.size() == values_references.size()
  for(size_t i = 0; i < labels.size(); ++i){
    AddDoubleLabelToPlot(labels[i], values[i].get());
  }
}

void MainWindow::AddDoubleLabelToPlot(const QString &label, double *value) {
  plot_labels_.push_back(std::make_shared<plotLabel>(ui->widget_plot, label_horizontal_offset_,
                                                     label_vertical_offset_, label, value,
                                                     std::make_shared<plotLabelDoubleDataPreparator>()));
  label_vertical_offset_ += label_vertical_offset_step_;
}

void MainWindow::AddIntLabelToPlot(const QString &label, int *value) {
  plot_labels_.push_back(std::make_shared<plotLabel>(ui->widget_plot, label_horizontal_offset_,
                                                     label_vertical_offset_, label, value,
                                                     std::make_shared<plotLabelIntDataPreparator>()));
  label_vertical_offset_ += label_vertical_offset_step_;
}

void MainWindow::AddConstantLabelToPlot(const QString &label) {
  plot_labels_.push_back(std::make_shared<plotLabel>(ui->widget_plot, label_horizontal_offset_,
                                                     label_vertical_offset_, label));
  label_vertical_offset_ += label_vertical_offset_step_;
}

void MainWindow::AddL1ErrorsToSum(QVector<ErrorsCalculator*> &errors_calculators, QVector<double> &errors_sums) {
  for(size_t i = 0; i < errors_calculators.size(); ++i){
    errors_sums[i] += errors_calculators[i]->CalculateL1Error();
  }
}

void MainWindow::AddL2ErrorsToSum(QVector<ErrorsCalculator*> &errors_calculators, QVector<double> &errors_sums) {
  for(size_t i = 0; i < errors_calculators.size(); ++i){
    errors_sums[i] += errors_calculators[i]->CalculateL1Error();
  }
}

void MainWindow::AddSupErrorsToSum(QVector<ErrorsCalculator*> &errors_calculators, QVector<double> &errors_sums) {
  for(size_t i = 0; i < errors_calculators.size(); ++i){
    errors_sums[i] += errors_calculators[i]->CalculateSupError();
  }
}

void MainWindow::AddModErrorsToSum(QVector<ErrorsCalculator*> &errors_calculators, QVector<double> &errors_sums) {
  for(size_t i = 0; i < errors_calculators.size(); ++i){
    errors_sums[i] += errors_calculators[i]->CalculateModError();
  }
}

void MainWindow::AddColorsLegendToPlot() {
  if(ui->checkBox_showEstimatedPlot->isChecked()) {
    AddConstantLabelToPlot("Model");
    plot_labels_.back()->SetColor(model_plot_pen_.color());
  }
  if(ui->checbox_showFullEstimator->isChecked()){
    AddConstantLabelToPlot("Window");
    plot_labels_.back()->SetColor(windowed_plot_pen_.color());
  }

  if(ui->checkBox_showEstimationPlot->isChecked()){
    AddConstantLabelToPlot("Less elements");
    plot_labels_.back()->SetColor(kde_plot_pen_.color());
  }

  if(ui->checkBox_showWeightedEstimationPlot->isChecked()){
    AddConstantLabelToPlot("Weights");
    plot_labels_.back()->SetColor(weighted_plot_pen_.color());
  }

  if(ui->checkBox_kernelPrognosedPlot->isChecked()){
    AddConstantLabelToPlot("Derivative");
    plot_labels_.back()->SetColor(derivative_plot_pen_.color());
  }

  if(ui->checkBox_standarizedDerivative->isChecked()){
    AddConstantLabelToPlot("Std. Derivative");
    plot_labels_.back()->SetColor(standardized_derivative_plot_pen_.color());
  }

  if(ui->checkBox_sigmoidallyEnhancedKDE->isChecked()){
    AddConstantLabelToPlot("Prediction");
    plot_labels_.back()->SetColor(desda_kde_plot_pen_.color());
  }

  if(ui->checkBox_REESEKDE->isChecked()){
    AddConstantLabelToPlot("Rare Elements");
    plot_labels_.back()->SetColor(desda_rare_elements_kde_plot_pen_.color());
  }
}
