#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <cfloat>
#include <QDebug>
#include <algorithm>
#include <chrono>
#include <QDateTime>
#include <Benchmarking/errorsCalculator.h>
#include <UI/plotLabelDoubleDataPreparator.h>
#include "KDE/pluginsmoothingparametercounter.h"

#include "kerDepCcWde.h"
#include "kerDepWindowedWde.h"

#include "UI/plotLabel.h"
#include "UI/plotLabelIntDataPreparator.h"

#include "Functions/complexfunction.h"
#include "Distributions/alternatingSplittingDistribution.h"

#include "Reservoir_sampling/biasedReservoirSamplingAlgorithm.h"
#include "Reservoir_sampling/basicReservoirSamplingAlgorithm.h"

#include "Reservoir_sampling/distributiondataparser.h"
#include "Reservoir_sampling/progressivedistributiondatareader.h"
#include "Reservoir_sampling/sinusoidalDistributionDataReader.h"
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
#include <QtGlobal>

#include "UI/QwtContourPlotUI.h"


#include <qwt_plot.h>
#include <qwt_symbol.h>
#include <qwt_legend.h>


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

  FillErrorIndicesColors();

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

void MainWindow::FillErrorIndicesColors() {
  error_indices_colors.push_back(windowed_plot_pen_.color());
  error_indices_colors.push_back(kde_plot_pen_.color());
  error_indices_colors.push_back(weighted_plot_pen_.color());
  error_indices_colors.push_back(desda_kde_plot_pen_.color());
  error_indices_colors.push_back(desda_rare_elements_kde_plot_pen_.color());
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

  int replace_constant = 0;

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

  // Generate m=m0 estimator plot
  if(ui->checbox_showFullEstimator->isChecked()) {
    auto windowed_estimator_values = DESDAAlgorithm->getWindowKDEValues(&drawable_domain);
    auto windowed_estimator_y = QVector<double>(windowed_estimator_values.begin(),
                                                windowed_estimator_values.end());

    for(auto i = 0; i < drawable_domain_.size(); ++i){
      drawable_domain_[i] += replace_constant;
    }

    AddPlot(&windowed_estimator_y, windowed_plot_pen_);
  }

  // Generate variable m estimator plot
  if(ui->checkBox_showEstimationPlot->isChecked()) {
    auto less_elements_estimator_values = DESDAAlgorithm->getKDEValues(&drawable_domain);
    auto less_elements_estimator_y = QVector<double>(less_elements_estimator_values.begin(),
                                                     less_elements_estimator_values.end());

    for(auto i = 0; i < drawable_domain_.size(); ++i){
      drawable_domain_[i] += replace_constant;
    }

    AddPlot(&less_elements_estimator_y, kde_plot_pen_);
  }

  // Generate weights estimator plot
  if(ui->checkBox_showWeightedEstimationPlot->isChecked()) {
    auto weighted_estimator_values = DESDAAlgorithm->getWeightedKDEValues(&drawable_domain);
    auto weighted_estimator_y = QVector<double>(weighted_estimator_values.begin(),
                                                weighted_estimator_values.end());

    for(auto i = 0; i < drawable_domain_.size(); ++i){
      drawable_domain_[i] += replace_constant;
    }

    AddPlot(&weighted_estimator_y, weighted_plot_pen_);
  }

  // Generate plot for prognosis estimator
  if(ui->checkBox_sigmoidallyEnhancedKDE->isChecked()) {
    auto prognosis_enhanced_plot_values = DESDAAlgorithm->getEnhancedKDEValues(&drawable_domain);
    auto prognosis_enhanced_plot_y = QVector<double>(prognosis_enhanced_plot_values.begin(),
                                                     prognosis_enhanced_plot_values.end());

    for(auto i = 0; i < drawable_domain_.size(); ++i){
      drawable_domain_[i] += replace_constant;
    }

    AddPlot(&prognosis_enhanced_plot_y, desda_kde_plot_pen_);
  }

  // Generate plot for atypical estimator
  if(ui->checkBox_REESEKDE->isChecked()) {
    auto rare_elements_enhanced_plot_values = DESDAAlgorithm->getRareElementsEnhancedKDEValues(&drawable_domain);
    auto rare_elements_enhanced_plot_y = QVector<double>(rare_elements_enhanced_plot_values.begin(),
                                                         rare_elements_enhanced_plot_values.end());

    for(auto i = 0; i < drawable_domain_.size(); ++i){
      drawable_domain_[i] += replace_constant;
    }

    AddPlot(&rare_elements_enhanced_plot_y, desda_rare_elements_kde_plot_pen_);
  }

  // Generate prognosis derivative plot
  if(ui->checkBox_kernelPrognosedPlot->isChecked()) {
    AddPlot(&kernel_prognosis_derivative_values_, derivative_plot_pen_);
  }

  // Generate plot for standardized prognosis derivative, assuming that normal derivative was generated first
  if(ui->checkBox_standarizedDerivative->isChecked()) {
    QVector<double> standardizedDerivativeY = {};

    double normalization_factor = qAbs(kernel_prognosis_derivative_values_[0]);

    for(auto val : kernel_prognosis_derivative_values_) {
      if(qAbs(val) > normalization_factor){
        normalization_factor = qAbs(val);
      }
    }

    for(auto val : kernel_prognosis_derivative_values_) {
      standardizedDerivativeY.push_back( 0.1 * val / normalization_factor);
    }
    AddPlot(&standardizedDerivativeY, standardized_derivative_plot_pen_);
  }

  if(ui->checkBox_showUnusualClusters->isChecked()) {
    MarkUncommonClusters(DESDAAlgorithm);
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

  int i = minX;
  int i_increment = 1;

  if(maxX - minX >= 100){
    i_increment = 10;
  } else {
    i_increment = 5;
  }

  for(i = minX; i <= maxX + 1; i += 1){
    ticks << i;
    if(i % i_increment == 0){
      labels << QString::number(i);
    }
    else{
      labels << "";
    }
  }

  QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
  textTicker->addTicks(ticks, labels);

  ui->widget_plot->xAxis->setTicker(textTicker);
  ui->widget_plot->axisRect()->setMinimumMargins(QMargins(0,0,30,0));
  //ui->widget_plot->xAxis->setTickLabelRotation(90);
}

void MainWindow::ClearPlot() {
  while(ui->widget_plot->graphCount() != 0)
    ui->widget_plot->removeGraph(0);

  for(auto a : lines_on_plot_) {
    ui->widget_plot->removeItem(a);
  }

  lines_on_plot_.clear();
}

unsigned long long MainWindow::MarkUncommonClusters(DESDA *DESDAAlgorithm) {

  atypical_elements_points_and_derivatives_ =
      DESDAAlgorithm->getAtypicalElementsValuesAndDerivatives();
  quantile_estimator_value_ = DESDAAlgorithm->_quantileEstimator;

  // Last element
  if(screen_generation_frequency_ == 1) {
    auto clusters = DESDAAlgorithm->getClustersForEstimator();
    auto verticalLine = new QCPItemLine(ui->widget_plot);
    verticalLine->start->setCoords(std::stod(clusters[0]->getRepresentative()->attributesValues["Val0"]), 0);
    verticalLine->end->setCoords(std::stod(clusters[0]->getRepresentative()->attributesValues["Val0"]), -0.02);
    verticalLine->setPen(QPen(Qt::black));
    lines_on_plot_.push_back(verticalLine);
  }

  // Atypical
  for(auto x : atypical_elements_points_and_derivatives_) {
    // Only works for distribution data
    auto verticalLine = new QCPItemLine(ui->widget_plot);
    verticalLine->start->setCoords(x.first[0], 0);
    verticalLine->end->setCoords(x.first[0], -0.01);
    if(x.second > 0)
      verticalLine->setPen(QPen(Qt::green));
    else
      verticalLine->setPen(QPen(Qt::red));
    lines_on_plot_.push_back(verticalLine);
  }

  return atypical_elements_points_and_derivatives_.size();
}

void MainWindow::MarkUncommonClusters2D(DESDA *DESDAAlgorithm, std::deque<QwtPlotCurve> *uncommon_clusters_markers){
  atypical_elements_points_and_derivatives_ = DESDAAlgorithm->getAtypicalElementsValuesAndDerivatives();
  quantile_estimator_value_ = DESDAAlgorithm->_quantileEstimator;

  QPolygonF new_trends;
  QPolygonF vanishing_trends;

  // Mark last cluster
  if(screen_generation_frequency_ == 1) {
    QPolygonF last_cluster_point;
    auto last_cluster = DESDAAlgorithm->getClustersForEstimator()[0];

    last_cluster_point << QPointF(std::stod(last_cluster->getRepresentative()->attributesValues["Val0"]),
                                  std::stod(last_cluster->getRepresentative()->attributesValues["Val1"]));

    uncommon_clusters_markers->at(2).setSamples(last_cluster_point);
  }

  for(auto x : atypical_elements_points_and_derivatives_) {
    // Only works for distribution data
    if(x.second > 0) {
      new_trends << QPointF(x.first[0], x.first[1]);
    } else {
      vanishing_trends << QPointF(x.first[0], x.first[1]);
    }
  }

  // On 0 we've got red +, on 1 we've got green +, on 2 we've got black X (for debug).
  uncommon_clusters_markers->at(0).setSamples(vanishing_trends);
  uncommon_clusters_markers->at(1).setSamples(new_trends);
}

QString MainWindow::FormatNumberForDisplay(double number) {
  // According to PK the number should be displayed as #.###
  QString result = " ";

  if(number < 0) result = "";

  QStringList splitNumber = QString::number(number, 'f', 7).split(".");
  result += splitNumber[0];

  if(splitNumber.size() == 1) return result;

  result += ".";

  for(int i = 0; i < 3 && i < splitNumber[1].size(); ++i)
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

  means->clear();

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
    vector<std::shared_ptr<vector<double>>> *stDevs,
    double additionalMultiplier) {

  int seed = ui->lineEdit_seed->text().toInt();
  vector<double> contributions;
  vector<std::shared_ptr<distribution>> elementalDistributions;

  int targetFunctionElementsNumber = ui->tableWidget_targetFunctions->rowCount();

  for(int functionIndex = 0; functionIndex < targetFunctionElementsNumber; ++functionIndex) {
    contributions.push_back
                     (
                         (dynamic_cast<QLineEdit *>(ui->tableWidget_targetFunctions
                                                      ->cellWidget(functionIndex,
                                                                   static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))))
                             ->text().toDouble()
                     );
    qDebug() << "Setting seed " << seed << ".";
    elementalDistributions.push_back(
        std::shared_ptr<distribution>(
            new normalDistribution(seed,
                                   (*means)[functionIndex].get(),
                                   (*stDevs)[functionIndex].get())));
    seed += 1;
  }

  return new alternatingSplittingDistribution(seed, &elementalDistributions, &contributions, additionalMultiplier);
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
    int dimensionsNumber, const bool &radial) {
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
      &kernelsIDs,
      radial
                                   );
}

function *MainWindow::GenerateTargetFunction(
    vector<std::shared_ptr<vector<double>>> *means,
    vector<std::shared_ptr<vector<double>>> *stDevs,
    double additionalMultiplier) {
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

  log("Start pushed!");
  // Delay so that
  QTime dieTime= QTime::currentTime().addSecs(0);
  while (QTime::currentTime() < dieTime) {
    QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
  }

  //RunAccuracyExperiment();

  // Set target function
  //SetBimodalTargetFunction();
  // SetTrimodalTargetFunction();

  // Set number of iterations
  this->ui->lineEdit_iterationsNumber->setText("10000");
  pcName = "sz";
  screen_generation_frequency_ = 1000;
  error_computation_frequency_ = 10;
  int n_seeds = 20;
  stream_number = 0;

  for(int seed = n_seeds - 19; seed < n_seeds + 1; ++seed){
  //for(int seed = n_seeds; seed > 0; --seed){  // Reversed loop for other experiments.
    ui->lineEdit_seed->setText(QString::number(seed)); // Default seed.
    ui->label_dataStream->setText("y:\\data\\stream_" + QString::number(stream_number) + "\\stream_" + QString::number(stream_number) + "_" + QString::number(seed) +  ".csv");
    //ui->label_dataStream->setText("k:\\Coding\\Python\\Poligon\\Articles\\IBS_PhD\\streams\\stream_13\\""\\stream_" + QString::number(stream_number) + "_" + QString::number(seed) +  ".csv");
    //Run1DExperimentWithDESDA();
    //Run1DExperimentWithClusterKernels();
    Run1DExperimentWithWDE();
    //Run1DExperimentWithSOMKE();
  }

  //run_3d_experiment();
}

void MainWindow::on_pushButton_clicked() {

  log("Start pushed!");
  // Delay so that
  QTime dieTime= QTime::currentTime().addSecs(0);
  while (QTime::currentTime() < dieTime) {
    QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
  }

  log("2D Experiment start.");
  bool radial = false;

  // Symbols
  std::deque<QwtSymbol> uncommon_clusters_symbols;
  uncommon_clusters_symbols.emplace_back(QwtSymbol::Cross, QBrush(Qt::black), QPen(Qt::red, 2), QSize(12, 12));
  uncommon_clusters_symbols.emplace_back(QwtSymbol::Cross, QBrush(Qt::black), QPen(Qt::green, 2), QSize(12, 12));
  uncommon_clusters_symbols.emplace_back(QwtSymbol::Cross, QBrush(Qt::black), QPen(Qt::black, 2), QSize(24, 24));

  std::deque<QwtPlotCurve> uncommon_clusters_markers;
  for(size_t i = 0; i < uncommon_clusters_symbols.size(); ++i){
    uncommon_clusters_markers.emplace_back("");
    uncommon_clusters_markers[i].setStyle(QwtPlotCurve::NoCurve);
    uncommon_clusters_markers[i].setSymbol( &(uncommon_clusters_symbols[i]) );
    uncommon_clusters_markers[i].attach(contour_plot_);
  }

  //curve.setSamples( points );
  //curve.attach( contour_plot_ );

  // Initially these vectors were used in errors computation only. We now also use them for the spectrogram.
  QVector<double> error_xs = {};
  QVector<double> error_ys = {};
  std::vector<double> model_function_values = {};
  std::vector<double> estimator_values = {};
  std::vector<std::vector<double>> error_domain = {};

  screen_generation_frequency_ = 10;
  int seed = ui->lineEdit_seed->text().toInt();

  // Contour levels calculation.
  QList<double> contourLevels;
  double level = 0.0003125;
  double level_multipliers = 2;
  while(level < 0.09) {
    contourLevels += level;
    level *= level_multipliers;
  }
  contourLevels += 0.155;

  // Automatic dimension update
  if(ui->spinBox_dimensionsNumber->value() == 1){
    ui->spinBox_dimensionsNumber->setValue(2);
    ui->spinBox_dimensionsNumber->editingFinished();
  }

  // Add clusters_ to the estimator
  means_ = {std::make_shared<std::vector<double >>()};
  means_.back()->push_back(0);
  means_.back()->push_back(0);

  standard_deviations_ = {std::make_shared<std::vector<double >>()};
  standard_deviations_.back()->push_back(1);
  standard_deviations_.back()->push_back(1);

  auto densityFunction =
      new multivariateNormalProbabilityDensityFunction(means_.back().get(), standard_deviations_.back().get(), 0.7);

  // Create estimator object
  std::shared_ptr<kernelDensityEstimator>
      estimator(GenerateKernelDensityEstimator(2, radial));

  estimator->_shouldConsiderWeights = true;

  std::vector<double> pt = {0, 0};
  contour_plot_->ShowColorMap(false);
  //*
  contour_plot_->addQwtPlotSpectrogram(new SpectrogramData2(estimator.get(), &error_xs, &error_ys, &estimator_values),
                                       QPen(QColor(0, 0, 0)));
  //*/

  std::shared_ptr<distribution> targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));
  std::vector<double> meansForDistribution = {0.0, 0.0};
  std::vector<double> stDevsForDistribution = {1.0, 1.0};

  parser_.reset(new distributionDataParser(&attributes_data_));

  QString expNum = "R84 (2D); Product; ; Correlation 0.7";
  QString pc_id = "sz285";
  int drawing_start_step = 0;
  int errors_calculation_start_step = 0;

  /*
  QString expDesc = "Rio 2014 Temp-Hum, start=" + QString::number(drawing_start_step) + ", " + pc_id; QString experiment_description = "Rio de Janeiro; 2014; temperature-humidity"; QDate data_start_date(2013, 10, 1); std::string data_path = "y:\\Data\\rio_2014_temp_humidity.csv";
  //QString expDesc = "id=1, Cracow 2020 Temp-Hum, start="+QString::number(drawing_start_step) + ", " + pc_id; QString experiment_description = "Cracow; 2020; temperature-humidity"; QDate data_start_date(2019, 10, 1); std::string data_path = "y:\\Data\\cracow_2020_temp_humidity.csv";
  //QString expDesc = "TEST, " + pc_id; QString experiment_description = "TEST"; QDate data_start_date(2019, 10, 1); std::string data_path = "y:\\Data\\2d_trend_v3.csv";

  QTime data_start_time(0, 0, 0);
  QDateTime data_date_time(data_start_date, data_start_time);

  ui->lineEdit_iterationsNumber->setText("15000");
  reader_.reset(new TextDataReader(data_path, 2));

  bool should_compute_errors = false;

  // Set limits on axes.
  contour_plot_->setAxesLimit(0); // This function doesn't work as the arguments suggest.

  //*/

  //*
  // p2 = 0.75p1 lub p2=0
  bool should_compute_errors = true;
  QString p2 = "0.5";

  // Prepare the reader
  reader_.reset(new progressiveDistributionDataReader(targetDistribution.get(), 0,0, new normalDistribution(0, &meansForDistribution, &stDevsForDistribution), p2.toDouble()));

  // Only to remove problems initialize the date
  QTime data_start_time(0, 0, 0); QDate data_start_date(2019, 10, 1); QDateTime data_date_time(data_start_date, data_start_time);

  if(p2.size() > 1){
    p2=p2+"p_1";
  }

  if(p2 == "1"){
    p2 = "p_1";
  }

  if(p2 == "0"){
    p2 = "0";
  }

  // Multiple instructions in one line, for simplicity
  QString experiment_description = "assumed data stream; skewed normal; 2D; p_2 = " + p2; QString expDesc = "assumed data stream 2D, p_2 = " + p2 + ", " + pc_id;

  // Set limits on axes.
  contour_plot_->setAxesLimit(20); // This function doesn't work as the arguments suggest.

  //*
  contour_plot_->addQwtPlotSpectrogram(new SpectrogramData2(densityFunction, &error_xs, &error_ys, &model_function_values), QPen(QColor(255, 0, 0)));
  //*/

  // After adding plots set contours and stuff.
  contour_plot_->setContours(contourLevels);
  contour_plot_->showContour(true);
  contour_plot_->setAlpha(0); // Set alpha to 0 for lack of coloring

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  reservoirSamplingAlgorithm *samplingAlgorithm =
      GenerateReservoirSamplingAlgorithm(reader_.get(), parser_.get());

  objects_.clear();
  clusters_ = &stored_medoids_;
  clusters_->clear();

  int pluginRank = 3;
  groupingThread gt(&stored_medoids_, parser_);

  derivative_estimator_.reset(GenerateKernelDensityEstimator(2, radial));
  enhanced_kde_.reset(GenerateKernelDensityEstimator(2, radial));

  DESDA DESDAAlgorithm(
      estimator,
      derivative_estimator_,
      enhanced_kde_,
      samplingAlgorithm,
      clusters_,
      ui->lineEdit_rarity->text().toDouble(), pluginRank
                      );

  // Start the test
  step_number_ = 0;

  time_t startTime, endTime;

  double l1_n_ = 0;
  double l2_n_ = 0;
  double sup_n_ = 0;
  double mod_n_ = 0;
  double actual_l1 = 0;
  double actual_l2 = 0;
  double actual_sup = 0;
  double actual_mod = 0;
  int errorCalculationsNumber = 0;
  double sum_l2 = 0;

  QwtContourPlotUI plotUi(&step_number_, screen_generation_frequency_, seed,
                          &DESDAAlgorithm, &l1_n_, &l2_n_, &sup_n_, &mod_n_,
                          &actual_l1, &actual_l2, &actual_sup, &actual_mod,
                          &data_date_time, &level_multipliers, experiment_description);
  plotUi.attach(contour_plot_);
  plotUi.updateTexts();
  plotUi.SetErrorsPrinting(should_compute_errors);
  QVector<int> initialDrawingSteps = {};

  for(size_t i = 500; i <= drawing_start_step; i += 500){
    initialDrawingSteps.push_back(i);
  }

  double domain_area = 0;

  ErrorsCalculator errors_calculator(&model_function_values, &estimator_values, &error_domain, &domain_area);

  // Prepare image location.
  this->setWindowTitle("Experiment #" + expNum);
  QString driveDir = "Y:\\"; // WIT PCs after update
  //QString driveDir = "D:\\Test\\"; // Home
  //QString driveDir = "d:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\";
  QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " (" + expDesc + ")\\";
  //QString dirPath = driveDir + "Eksperyment " + expNum + " (" + expDesc + ")\\";
  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  int steps_number = ui->lineEdit_iterationsNumber->text().toInt();

  log("Experiment started.");
  for(step_number_ = 1; step_number_ <= steps_number; ++step_number_) {

    log("New step.");
    startTime = time(nullptr);

    log("Performing new step.");
    DESDAAlgorithm.performStep();
    log("Step performed.");

    bool compute_errors = step_number_ >= errors_calculation_start_step && should_compute_errors;
    bool draw_plot = (step_number_ % screen_generation_frequency_ == 0 && step_number_ >= drawing_start_step) || initialDrawingSteps.contains(step_number_);

    if(ui->checkBox_showUnusualClusters->isChecked() && draw_plot){
      log("Marking rare elements.");
      MarkUncommonClusters2D(&DESDAAlgorithm, &uncommon_clusters_markers);
    }

    log("Estimator preparation.");
    DESDAAlgorithm.prepareEstimatorForContourPlotDrawing();
    log("Estimator preparation finished.");

    // NOTE: We use error domain for spectrogram generation! That's why we compute the domain and values outside the if.
    if(compute_errors || draw_plot) {
      log("Computing domains.");
      error_xs = DESDAAlgorithm.getErrorDomain(0);
      error_ys = DESDAAlgorithm.getErrorDomain(1);
      error_domain = Generate2DPlotErrorDomain(error_xs, error_ys);
      log("Computing values of domains.");
      model_function_values = GetFunctionsValueOnDomain(densityFunction, error_domain);
      estimator_values = GetFunctionsValueOnDomain(estimator.get(), error_domain);
    }

    // Error calculation
    //*
    if(compute_errors) {
      log("Error calculation started.");
      ++errorCalculationsNumber;

      domain_area = Calculate2DDomainArea(error_domain);

      //actual_l1 = errors_calculator.CalculateL1Error();
      actual_l2 = errors_calculator.CalculateL2Error();
      //actual_sup = errors_calculator.CalculateSupError();
      //actual_mod = errors_calculator.CalculateModError();
      //sum_l1 += actual_l1;
      sum_l2 += actual_l2;
      //sum_sup += actual_sup;
      //sum_mod = actual_mod;

      //l1_n_ = sum_l1 / errorCalculationsNumber;
      l2_n_ = sum_l2 / errorCalculationsNumber;
      //sup_n_ = sum_sup / errorCalculationsNumber;
      //mod_n_ = sum_mod / errorCalculationsNumber;
      log("Error calculation finished.");
    }
    //*/
    if(should_compute_errors) {
      densityFunction->setMeans(*means_.back().get());
    }

    // Drawing
    if(draw_plot){

      log("Drawing started.");

      log("Texts updates.");
      plotUi.updateTexts();

      log("Replotting.");
      contour_plot_->replot();

      endTime = time(nullptr);

      log("Processing.");

      QCoreApplication::processEvents();

      QString imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image name: " + imageName);
      log("Saved: " + QString::number(ui->widget_contour_plot_holder->grab().save(imageName)));
      log("Drawing finished.");

    }

    log("Restoring weights.");
    DESDAAlgorithm.restoreClustersCWeights();

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

std::vector<std::vector<double>> MainWindow::Generate2DPlotErrorDomain(const QVector<double> &xDomainValues,
                                                                       const QVector<double> &yDomainValues) {
  std::vector<point> domainValues = {};

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


void MainWindow::SetBimodalTargetFunction(){
    // 0.6 N(0, 1) + 0.4 N(5, 1)
    ui->pushButton_addTargetFunction->click();
    // log("Klikłem");

    auto meansTablePointer =
        dynamic_cast<QTableWidget *>(ui->tableWidget_targetFunctions->cellWidget(0,
                                                                                 static_cast<int>(TargetFunctionSettingsColumns::kMeanColumnIndex)));
    (dynamic_cast<QLineEdit *>(meansTablePointer->cellWidget(0, 0)))->setText("5.0");

    //log("Ustaiwam drugie.");

    dynamic_cast<QLineEdit *>(
        ui->tableWidget_targetFunctions
            ->cellWidget(0, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))
        )->setText(QString::number(40));

    UpdateLastContribution();
    // log("Działa");

    // Update the plot bounds
    ui->lineEdit_maxX->setText("53");
    ui->lineEdit_minX->setText("-5");
}

void MainWindow::SetTrimodalTargetFunction(){
    // 0.4 N(0, 1) + 0.3 N(5, 1) + 0.3 N(-5, 1)
    ui->pushButton_addTargetFunction->click();
    ui->pushButton_addTargetFunction->click();

    //log("Klikłem");

    auto meansTablePointer =
        dynamic_cast<QTableWidget *>(ui->tableWidget_targetFunctions->cellWidget(0,
                                                                                 static_cast<int>(TargetFunctionSettingsColumns::kMeanColumnIndex)));
    (dynamic_cast<QLineEdit *>(meansTablePointer->cellWidget(0, 0)))->setText("5.0");

    meansTablePointer =
        dynamic_cast<QTableWidget *>(ui->tableWidget_targetFunctions->cellWidget(1,
                                                                                 static_cast<int>(TargetFunctionSettingsColumns::kMeanColumnIndex)));

    (dynamic_cast<QLineEdit *>(meansTablePointer->cellWidget(0, 0)))->setText("-5.0");

    //log("Ustaiwam drugie.");

    dynamic_cast<QLineEdit *>(
        ui->tableWidget_targetFunctions
            ->cellWidget(0, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))
        )->setText(QString::number(30));

    dynamic_cast<QLineEdit *>(
        ui->tableWidget_targetFunctions
            ->cellWidget(1, static_cast<int>(TargetFunctionSettingsColumns::kContributionColumnIndex))
        )->setText(QString::number(30));

    UpdateLastContribution();
    // log("Działa");

    // Update the plot bounds
    ui->lineEdit_maxX->setText("53");
    ui->lineEdit_minX->setText("-10");
}

void MainWindow::Run1DExperimentWithDESDA() {

  ui->widget_plot->clearGraphs();
  stored_medoids_.clear();

  int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

  if(!CanAnimationBePerformed(dimensionsNumber)) return;

  QString seedString = ui->lineEdit_seed->text();

  // Log that application started generating KDE
  // Standard seed was 5625.
  log("KDE animation started.");
  log("Seed: " + seedString);
  log("Sample size: " + ui->lineEdit_sampleSize->text());
  log("Experiment started.");  

  srand(static_cast<unsigned int>(seedString.toInt()));

  FillMeans(&means_);
  FillStandardDeviations(&standard_deviations_);

  target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

  std::shared_ptr<kernelDensityEstimator>
      estimator(GenerateKernelDensityEstimator(dimensionsNumber));

  estimator->_shouldConsiderWeights = false;

  derivative_estimator_.reset(GenerateKernelDensityEstimator(dimensionsNumber));
  enhanced_kde_.reset(GenerateKernelDensityEstimator(dimensionsNumber));

  double evenDistributionsMultiplier = 1;  // 1 to standardowa "ścieżka zdrowia / życia"

  std::shared_ptr<distribution>
      targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_, evenDistributionsMultiplier));

  parser_.reset(new distributionDataParser(&attributes_data_));
  std::vector<double> target;

  bool should_compute_errors = true;
  //double p2 = 0.1;

  QString streamDesc = "assumed trimodal";
  QString expDesc = "id=" + QString::number(screen_generation_frequency_) + ", "+streamDesc+" data stream, seed=" + seedString;
  QString plot_description = streamDesc + " data stream";
  ui->checkBox_showEstimatedPlot->setChecked(true);

  int drawing_start_step = 0;

  QString expNum = "T" + QString::number(stream_number) + "_" + seedString;

  expDesc += ", " + pcName;

  if(this->ui->label_dataStream->text().toStdString() == "Not selected."){
    log("Data stream not selected.");
    return;
  }

  std::string data_path = this->ui->label_dataStream->text().toStdString();
  reader_.reset(new TextDataReader(data_path));

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  reservoirSamplingAlgorithm *algorithm =
      GenerateReservoirSamplingAlgorithm(reader_.get(), parser_.get());

  objects_.clear();

  int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();

  log("Attributes data set.");

  clusters_ = &stored_medoids_;

  derivative_estimator_->_shouldConsiderWeights = false;

  int pluginRank = 2;
  DESDA DESDAAlgorithm(
      estimator,
      derivative_estimator_,
      enhanced_kde_,
      algorithm,
      clusters_,
      ui->lineEdit_rarity->text().toDouble(), pluginRank
                      );

  this->setWindowTitle("Experiment #" + expNum);

  ClearPlot();
  ResizePlot();

  set_paths(expNum, expDesc);

  QString imageName = "";

  // Set labels
  plot_labels_ = {};
  label_horizontal_offset_ = 0.02;
  label_vertical_offset_ = 0.01;

  plotLabel desc_label(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_,
                       plot_description);
  label_vertical_offset_ += label_vertical_offset_step_;
  label_vertical_offset_ += label_vertical_offset_step_;

  QVector<plotLabel> date_labels = {};
  //*
  if(!should_compute_errors) {
    date_labels.push_back(plotLabel(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_, ""));
    label_vertical_offset_ += label_vertical_offset_step_;
  }
  //*/

  AddIntLabelToPlot("t          = ", &step_number_);
  // AddConstantLabelToPlot("iw    = " + QString::number(screen_generation_frequency_));
  // AddConstantLabelToPlot("seed  = " + seedString);
  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabel KPSSTextLabel(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_, "KPSS         = 0");
  label_vertical_offset_ += label_vertical_offset_step_;

  AddDoubleLabelToPlot("sgmKPSS    =", &DESDAAlgorithm._sgmKPSS);
  label_vertical_offset_ += label_vertical_offset_step_;

  //AddConstantLabelToPlot("mKPSS = " + QString::number(DESDAAlgorithm._kpssM));
  AddIntLabelToPlot("m          = ", &(DESDAAlgorithm._m));
  //AddConstantLabelToPlot("m_min      = " + QString::number(DESDAAlgorithm._minM));
  //AddConstantLabelToPlot("m_0        = " + ui->lineEdit_sampleSize->text());
  label_vertical_offset_ += label_vertical_offset_step_;

  //AddDoubleLabelToPlot("beta0 = ", &(DESDAAlgorithm._beta0));
  //label_vertical_offset_ += label_vertical_offset_step_;

  AddDoubleLabelToPlot("r          =", &(DESDAAlgorithm._r));
  //AddDoubleLabelToPlot("q          =", &(DESDAAlgorithm._quantileEstimator));
  AddIntLabelToPlot("#atypical  = ", &(DESDAAlgorithm._rareElementsNumber));
  label_vertical_offset_ += label_vertical_offset_step_;
  //AddIntLabelToPlot("trend = ", &(DESDAAlgorithm._trendsNumber));
  //*
  QString signal_exclamation_points = "                    ";
  plotLabel signal_exclamation_point_label(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_,
                                           signal_exclamation_points);

  //AddDoubleLabelToPlot("TS         =" , &(DESDAAlgorithm.statistics_[0]));
  plotLabel statisticsTextLabel(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_, "TS           = 0");
  label_vertical_offset_ += label_vertical_offset_step_;

  //label_vertical_offset_ += label_vertical_offset_step_;
  //label_vertical_offset_ += 3 * label_vertical_offset_step_;

  AddColorsLegendToPlot();

  //====================  SECOND COLUMN =================//

  label_horizontal_offset_ = 0.25;
  label_vertical_offset_ = 0.01 + 9 * label_vertical_offset_step_;

  //==================== ERRORS SUM =================//

  label_horizontal_offset_ = 0.875;
  label_vertical_offset_ = 0.01;

  clear_errors();

  if(should_compute_errors) {

    //QVector<QString> l1_labels = {"L1_w  = ", "L1_m  = ", "L1_d  = ", "L1_p  = ", "L1_n  = "};
    QVector<QString> l2_labels = {"L2 =", "L2 =", "L2 =", "L2 =", "L2 ="};
    //QVector<QString> sup_labels = {"sup_w = ", "sup_m = ", "sup_d = ", "sup_p = ", "sup_n = "};
    //QVector<QString> mod_labels = {"mod_w = ", "mod_m = ", "mod_d = ", "mod_p = ", "mod_n = "};

    for(size_t i = 0; i < l2_labels.size(); ++i){
      //l1_errors.push_back(std::make_shared<double>(0));
      //l1_errors_sums.push_back(0);
      l2_errors.push_back(std::make_shared<double>(0));
      l2_errors_sums.push_back(0);
      //sup_errors.push_back(std::make_shared<double>(0));
      //sup_errors_sums.push_back(0);
      //mod_errors.push_back(std::make_shared<double>(0));
      //mod_errors_sums.push_back(0);
    }

    //AddErrorLabelsToPlot(l1_labels, l1_errors);
    //label_vertical_offset_ += label_vertical_offset_step_;

    AddErrorLabelsToPlot(l2_labels, l2_errors);
    label_vertical_offset_ += label_vertical_offset_step_;

    //AddErrorLabelsToPlot(sup_labels, sup_errors);
    //label_vertical_offset_ += label_vertical_offset_step_;

    //AddErrorLabelsToPlot(mod_labels, mod_errors);
    //label_vertical_offset_ += label_vertical_offset_step_;
  }

  FillDomain(&domain_, nullptr);
  for(const auto &pt : domain_){ drawable_domain_.push_back(pt->at(0)); }

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

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

  errors_calculators = {
      &windowed_errors_calculator, &less_elements_kde_errors_calculator, &weighted_kde_errors_calculator,
      &enhanced_kde_errors_calculator, &rare_elements_kde_errors_calculator
  };

  std::ifstream inFile(avg_l2_errors_file_path);
  numberOfErrorCalculations = std::count(std::istreambuf_iterator<char>(inFile),
             std::istreambuf_iterator<char>(), '\n');

  drawing_start_step = numberOfErrorCalculations * error_computation_frequency_;

  if(drawing_start_step > 0 && should_compute_errors) {
    // Open the file with average l2 errors and load them into the list.

    qDebug() << "Loading errors.";

    std::ifstream l2_errors_sums_file(l2_errors_sum_file_path);
    if(l2_errors_sums_file.is_open()) {
      l2_errors_sums.clear();
      std::string line;
      while(std::getline(l2_errors_sums_file, line)) {
        l2_errors_sums.push_back(std::stod(line));
      }
      l2_errors_sums_file.close();
    }
    else {
      log("Unable to open avg_l2_errors.txt file.");
    }
  }

  // Check if the final figure is already in the folder
  QString final_figure_path = dirPath + QString::number(stepsNumber) + ".png";
  if(drawing_start_step > 0 && should_compute_errors) {
      if(std::ifstream(final_figure_path.toStdString())) {
      log("Final figure already exists. Skipping.");
      return;
      }
  }

  ui->checkBox_kernelPrognosedPlot->toggle(); ui->checbox_showFullEstimator->toggle();

  for(step_number_ = 1; step_number_ <= stepsNumber; ++step_number_) {

    update_theoretical_density();

    DESDAAlgorithm.performStep();

    target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

    // Error calculations
    //*
    if(step_number_ > drawing_start_step && should_compute_errors && step_number_ % error_computation_frequency_ == 0) {

      ++numberOfErrorCalculations;

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

      // Save averaged l2 errors to the file
      std::ofstream in;
      in.open(avg_l2_errors_file_path, std::ios_base::app);
      for(size_t i = 0; i < errors_calculators.size(); ++i) {
        in << *l2_errors[i] << ";";
      }
      in << "\n";
      in.close();

      // Save the errors sum to the file
      std::ofstream in2;
      in2.open(l2_errors_sum_file_path);
      for(size_t i = 0; i < errors_calculators.size(); ++i) {
        in2 << l2_errors_sums[i] << "\n";
      }
      in2.close();
    }
    //*/

    // Drawing
    if(drawing_start_step < step_number_ && ( step_number_ % screen_generation_frequency_ == 0)) {
      log("Drawing in step number " + QString::number(step_number_) + ".");

      QVector<std::vector<double>> points_to_compute_derivative_on = {};

      for(double val : drawable_domain_){
        points_to_compute_derivative_on.push_back(std::vector<double>());
        points_to_compute_derivative_on.back().push_back(val);
      }

      kernel_prognosis_derivative_values_ =
          DESDAAlgorithm.getKernelPrognosisDerivativeValues(&points_to_compute_derivative_on);

      // ============= LEFT SIDE UPDATE ================ //

      KPSSTextLabel.setText("KPSS       =" + FormatNumberForDisplay(
          DESDAAlgorithm.getStationarityTestValue()));

      if(step_number_ >= drawing_start_step) {
        DrawPlots(&DESDAAlgorithm);
      }

      for(const auto &label : plot_labels_) {
        label->updateText();
      }

      signal_exclamation_points = "                   ";

      statisticsTextLabel.setText("TS         =" + FormatNumberForDisplay(
          DESDAAlgorithm.statistics_[0]));

      if(fabs(DESDAAlgorithm.statistics_[0]) >= 0.3){
        signal_exclamation_points += "!";
      } else if(fabs(DESDAAlgorithm.statistics_[0]) >= 0.2){
        signal_exclamation_points += "?";
      }

      signal_exclamation_point_label.setText(signal_exclamation_points);

      ui->widget_plot->replot();
      QCoreApplication::processEvents();
      imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
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

  int number_of_cluster_kernels = 100;
  double h = 1;
  double sigma = 0;

  // Creating target function.
  FillMeans(&means_);
  FillStandardDeviations(&standard_deviations_);
  target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

  std::shared_ptr<distribution>
      targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));

  parser_.reset(new distributionDataParser(&attributes_data_));

  QString expNum = "CK" + QString::number(stream_number) + "_" + seedString;
  this->setWindowTitle("Experiment #" + expNum);

  QString streamDesc = "assumed trimodal";
  QString expDesc = "id=" + QString::number(screen_generation_frequency_) + ", "+streamDesc+" data stream, seed=" + seedString;
  expDesc += ", " + pcName;

  if(this->ui->label_dataStream->text().toStdString() == "Not selected."){
      log("Data stream not selected.");
      return;
  }

  std::string data_path = this->ui->label_dataStream->text().toStdString();
  reader_.reset(new TextDataReader(data_path));

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();

  log("Attributes data set.");

  ClearPlot();
  ResizePlot();

  set_paths(expNum, expDesc);

  clear_errors();
  l2_errors.push_back(std::make_shared<double>(0));
  l2_errors_sums.push_back(0);

  bool should_compute_errors = true;

  // Setting up the labels
  QVector<std::shared_ptr<plotLabel>> plotLabels = {};
  double horizontalOffset = 0.01, verticalOffset = 0.01, verticalStep = 0.03;

  plotLabel desc_label(ui->widget_plot, horizontalOffset, verticalOffset,
                       streamDesc + " data stream; 1D");

  verticalOffset += 2 * verticalStep;

  QVector<plotLabel> date_labels = {};
  //*
  if(!should_compute_errors) {
    date_labels.push_back(plotLabel(ui->widget_plot, horizontalOffset, verticalOffset, ""));
    verticalOffset += verticalStep;
  }
  //*/

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   horizontalOffset, verticalOffset, "t     = ", &step_number_,
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  verticalOffset += 2 * verticalStep;

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

  label_vertical_offset_ = 0.75;
  label_horizontal_offset_ = horizontalOffset;

  AddConstantLabelToPlot("theoretical");
  plot_labels_.back()->SetColor(model_plot_pen_.color());

  AddConstantLabelToPlot("CK Estimator");
  plot_labels_.back()->SetColor(kde_plot_pen_.color());


  //====================  SECOND COLUMN =================//

  horizontalOffset = 0.20;
  verticalOffset = 0.01 + 9 * verticalStep;

  //====================== ERRORS SUM ===================//

  horizontalOffset = 0.85;
  verticalOffset = 0.75 + verticalStep;

  plotLabel L2TextLabel(ui->widget_plot, horizontalOffset, verticalOffset, "L2   = 0");
  L2TextLabel.SetColor(kde_plot_pen_.color());

  FillDomain(&domain_, nullptr);
  for(const auto &pt : domain_){
    drawable_domain_.push_back(pt->at(0));
  }

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

  double error_domain_length = 0;
  std::vector<std::vector<double>> error_domain = {};
  std::vector<double> model_values = {};
  std::vector<double> kde_values = {};

  ErrorsCalculator errors_calculator(
      &model_values, &kde_values, &error_domain, &error_domain_length
                                    );

  errors_calculators.append(&errors_calculator);

  std::ifstream inFile(avg_l2_errors_file_path);
  numberOfErrorCalculations = std::count(std::istreambuf_iterator<char>(inFile),
                                         std::istreambuf_iterator<char>(), '\n');

  int drawing_start_step = 0;
  drawing_start_step = numberOfErrorCalculations * error_computation_frequency_;

  if(drawing_start_step > 0 && should_compute_errors) {
      // Open the file with average l2 errors and load them into the list.

      qDebug() << "Loading errors.";

      std::ifstream l2_errors_sums_file(l2_errors_sum_file_path);
      if(l2_errors_sums_file.is_open()) {
          l2_errors_sums.clear();
          std::string line;
          while(std::getline(l2_errors_sums_file, line)) {
              l2_errors_sums.push_back(std::stod(line));
          }
          l2_errors_sums_file.close();
      }
      else {
          log("Unable to open avg_l2_errors.txt file.");
      }
  }

  log("Crating CK Algorithm!");
  auto CKAlgorithm = EnhancedClusterKernelAlgorithm(number_of_cluster_kernels,
                                                    CreateNewVarianceBasedClusterKernel);


  for(step_number_ = 1; step_number_ <= stepsNumber; ++step_number_) {

    update_theoretical_density();

    Point stream_value = {};
    reader_->getNextRawDatum(&stream_value);
    UnivariateStreamElement element(stream_value);

    log("Performing step: " + QString::number(step_number_));
    CKAlgorithm.PerformStep(&element);
    log("Step performed.");

    target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

    // Error calculations
    if(step_number_ > drawing_start_step && should_compute_errors && step_number_ % error_computation_frequency_ == 0) {
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
      ++numberOfErrorCalculations;
      compute_errors();

      log("Errors calculated.");
      log("Saving errors.");

      // Save averaged l2 errors to the file
      std::ofstream in;
      in.open(avg_l2_errors_file_path, std::ios_base::app);
      for(size_t i = 0; i < errors_calculators.size(); ++i) {
          in << *l2_errors[i] << ";";
      }
      in << "\n";
      in.close();

      // Save the errors sum to the file
      std::ofstream in2;
      in2.open(l2_errors_sum_file_path);
      for(size_t i = 0; i < errors_calculators.size(); ++i) {
          in2 << l2_errors_sums[i] << "\n";
      }
      in2.close();

      log("Errors saved.");
    }

    // Drawing
    if(drawing_start_step < step_number_ && ( step_number_ % screen_generation_frequency_ == 0)) {
      log("Drawing in step number " + QString::number(step_number_) + ".");

      L2TextLabel.setText("L2   =" + FormatNumberForDisplay(*l2_errors[0]));

      DrawPlots(&CKAlgorithm);

      h = CKAlgorithm.GetBandwidth();
      sigma = CKAlgorithm.GetStandardDeviation();

      for(const auto &label : plotLabels){
          label->updateText();
      }

      ui->widget_plot->replot();
      QCoreApplication::processEvents();

      QString imageName = dirPath + QString::number(step_number_) + ".png";
      log("Image saved: " + QString::number(ui->widget_plot->savePng(imageName, 0, 0, 1, -1)));
    }

  }

  log("Animation finished.");
}

void MainWindow::Run1DExperimentWithWDE() {
  // TR TODO: This is basically the same as it is in Cluster Kernels... Initialization should
  //          be made an separate function.
  int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

  if(!CanAnimationBePerformed(dimensionsNumber)) return;

  QString seedString = ui->lineEdit_seed->text();

  // Log that application started generating KDE
  // Standard seed was 5625.
  log("KDE animation with WDE started.");
  log("Seed: " + seedString);
  log("Sample size: " + ui->lineEdit_sampleSize->text());

  double weight_modifier = 0.99; // omega
  unsigned int maximal_number_of_coefficients = 100; // M
  unsigned int current_coefficients_number = 0; // #coef
  int number_of_elements_per_block = 500; // b

  srand(static_cast<unsigned int>(seedString.toInt()));

  // Creating target function.
  FillMeans(&means_);
  FillStandardDeviations(&standard_deviations_);
  target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

  std::shared_ptr<distribution>
      targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));

  parser_.reset(new distributionDataParser(&attributes_data_));

  QString expNum = "WDE" + QString::number(stream_number) + "_" + seedString;
  this->setWindowTitle("Experiment #" + expNum);

  QString streamDesc = "assumed";
  QString expDesc = "id=" + QString::number(screen_generation_frequency_) + ", "+streamDesc+" data stream, seed=" + seedString;
  expDesc += ", " + pcName;

  if(this->ui->label_dataStream->text().toStdString() == "Not selected."){
      log("Data stream not selected.");
      return;
  }

  std::string data_path = this->ui->label_dataStream->text().toStdString();
  reader_.reset(new TextDataReader(data_path));

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();

  log("Attributes data set.");

  ClearPlot();
  ResizePlot();

  set_paths(expNum, expDesc);

  clear_errors();
  l2_errors.push_back(std::make_shared<double>(0));
  l2_errors_sums.push_back(0);

  bool should_compute_errors = true;


  // Setting up the labels
  QVector<std::shared_ptr<plotLabel>> plotLabels = {};
  label_vertical_offset_ = label_horizontal_offset_ = 0.01;

  plotLabel desc_label(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_,
                       expDesc);

  label_vertical_offset_ += 2 * label_vertical_offset_step_;

  QVector<plotLabel> date_labels = {};
  //*
  if(!should_compute_errors) {
    date_labels.push_back(plotLabel(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_, ""));
    label_vertical_offset_ += label_vertical_offset_step_;
  }
  //*/

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "t     = ", &step_number_,
                                                   std::make_shared<plotLabelIntDataPreparator>()));

  label_vertical_offset_ += 2 * label_vertical_offset_step_;


  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "b     = ",
                                                   &(number_of_elements_per_block),
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "M     = ",
                                                   &(maximal_number_of_coefficients),
                                                   std::make_shared<plotLabelIntDataPreparator>()));

  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "omega = 0.99"));

  label_vertical_offset_ += 2 * label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "#coef = ",
                                                   &(current_coefficients_number),
                                                   std::make_shared<plotLabelIntDataPreparator>()));

  label_vertical_offset_ = 0.75;

  AddConstantLabelToPlot("theoretical");
  plot_labels_.back()->SetColor(model_plot_pen_.color());

  AddConstantLabelToPlot("WDE Estimator");
  plot_labels_.back()->SetColor(kde_plot_pen_.color());


  //====================  SECOND COLUMN =================//

  label_horizontal_offset_ = 0.20;
  label_vertical_offset_ = 0.01 + 9 * label_vertical_offset_step_;

  //====================== ERRORS SUM ===================//

  label_horizontal_offset_ = 0.85;
  label_vertical_offset_ = 0.78;

  plotLabel L2TextLabel(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_, "L2   = 0");
  L2TextLabel.SetColor(kde_plot_pen_.color());

  if(!should_compute_errors){
    L2TextLabel.setText("");
  }

  FillDomain(&domain_, nullptr);
  for(const auto &pt : domain_){
      drawable_domain_.push_back(pt->at(0));
  }

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

  numberOfErrorCalculations = 0;

  double error_domain_length = 0;
  std::vector<std::vector<double>> error_domain = {};
  std::vector<double> model_values = {};
  std::vector<double> wde_values = {};

  ErrorsCalculator errors_calculator(
      &model_values, &wde_values, &error_domain, &error_domain_length
                                    );

  errors_calculators = {&errors_calculator};

  std::ifstream inFile(avg_l2_errors_file_path);
  numberOfErrorCalculations = std::count(std::istreambuf_iterator<char>(inFile),
                                         std::istreambuf_iterator<char>(), '\n');

  int drawing_start_step = 0;
  drawing_start_step = numberOfErrorCalculations * error_computation_frequency_;

  if(drawing_start_step > 0 && should_compute_errors) {
      // Open the file with average l2 errors and load them into the list.

      qDebug() << "Loading errors.";

      std::ifstream l2_errors_sums_file(l2_errors_sum_file_path);
      if(l2_errors_sums_file.is_open()) {
          l2_errors_sums.clear();
          std::string line;
          while(std::getline(l2_errors_sums_file, line)) {
              l2_errors_sums.push_back(std::stod(line));
          }
          l2_errors_sums_file.close();
      }
      else {
          log("Unable to open avg_l2_errors.txt file.");
      }
  }

  Windowed_WDE WDE_Algorithm = Windowed_WDE(maximal_number_of_coefficients, weight_modifier,
                                            CreateWeightedThresholdedWaveletDensityEstimatorFromBlock,
                                            //CreateWeightedWaveletDensityEstimatorFromBlock,
                                            number_of_elements_per_block);

  for(step_number_ = 1; step_number_ <= stepsNumber; ++step_number_) {

    update_theoretical_density();

    Point stream_value = {};
    reader_->getNextRawDatum(&stream_value);

    log("Performing step: " + QString::number(step_number_));
    WDE_Algorithm.PerformStep(&stream_value);
    log("Step performed.");

    target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

    // Error calculations
    if(step_number_ > drawing_start_step && should_compute_errors && step_number_ % error_computation_frequency_ == 0) {

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
      ++numberOfErrorCalculations;
      compute_errors();

      log("Errors calculated.");
      log("Saving errors.");

      // Save averaged l2 errors to the file
      std::ofstream in;
      in.open(avg_l2_errors_file_path, std::ios_base::app);
      for(size_t i = 0; i < errors_calculators.size(); ++i) {
          in << *l2_errors[i] << ";";
      }
      in << "\n";
      in.close();

      // Save the errors sum to the file
      std::ofstream in2;
      in2.open(l2_errors_sum_file_path);
      for(size_t i = 0; i < errors_calculators.size(); ++i) {
          in2 << l2_errors_sums[i] << "\n";
      }
      in2.close();

      log("Errors saved.");
    }

    // Drawing
    if(drawing_start_step < step_number_ && ( step_number_ % screen_generation_frequency_ == 0)) {
      log("Drawing in step number " + QString::number(step_number_) + ".");

      L2TextLabel.setText("L2   =" + FormatNumberForDisplay(*l2_errors[0]));

      if(!should_compute_errors){
        L2TextLabel.setText("");
      }

      current_coefficients_number = WDE_Algorithm.GetCurrentCoefficientsNumber();

      DrawPlots(&WDE_Algorithm);

      for(const auto &label : plotLabels){
        label->updateText();
      }

      ui->widget_plot->replot();
      QCoreApplication::processEvents();

      QString imageName = dirPath + QString::number(step_number_) + ".png";
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
  log("KDE animation with SOMKE started.");
  log("Seed: " + seedString);
  log("Sample size: " + ui->lineEdit_sampleSize->text());

  step_number_ = 0;

  int neurons_number = 100;
  int epochs_number = 3000;
  int data_window_size = 500;
  double beta = 0;
  double alpha = 1.0;

  double sigma0 = 25.0;
  double tau1 = 1000 / log(sigma0);
  double tau2 = 1000.0;
  double eta0 = 3.0;

  srand(static_cast<unsigned int>(seedString.toInt()));

  // Creating target function.
  FillMeans(&means_);
  FillStandardDeviations(&standard_deviations_);
  target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

  std::shared_ptr<distribution>
      targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));

  parser_.reset(new distributionDataParser(&attributes_data_));

  QString expNum = "SOMKE" + QString::number(stream_number) + "_" + seedString;
  this->setWindowTitle("Experiment #" + expNum);

  QString streamDesc = "assumed trimodal";
  QString expDesc = "id=" + QString::number(screen_generation_frequency_) + ", "+streamDesc+" data stream, seed=" + seedString;
  expDesc += ", " + pcName;

  if(this->ui->label_dataStream->text().toStdString() == "Not selected."){
      log("Data stream not selected.");
      return;
  }

  std::string data_path = this->ui->label_dataStream->text().toStdString();
  reader_.reset(new TextDataReader(data_path));

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  int stepsNumber = ui->lineEdit_iterationsNumber->text().toInt();

  log("Attributes data set.");

  ClearPlot();
  ResizePlot();

  set_paths(expNum, expDesc);

  clear_errors();
  l2_errors.push_back(std::make_shared<double>(0));
  l2_errors_sums.push_back(0);

  bool should_compute_errors = true;

  // Setting up the labels
  QVector<std::shared_ptr<plotLabel>> plotLabels = {};
  label_horizontal_offset_ = label_vertical_offset_ = 0.01;

  plotLabel desc_label(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_,
                       expDesc);

  label_vertical_offset_ += 2 * label_vertical_offset_step_;

  QVector<plotLabel> date_labels = {};
  //*
  if(!should_compute_errors) {
    date_labels.push_back(plotLabel(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_, ""));
    label_vertical_offset_ += label_vertical_offset_step_;
  }
  //*/

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "t     = ", &step_number_,
                                                   std::make_shared<plotLabelIntDataPreparator>()));
  label_vertical_offset_ += 2 * label_vertical_offset_step_;


  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "#neurons = ", &(neurons_number),
                                                   std::make_shared<plotLabelIntDataPreparator>()));

  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "#epochs  = ", &(epochs_number),
                                                   std::make_shared<plotLabelIntDataPreparator>()));

  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "b        = ", &(data_window_size),
                                                   std::make_shared<plotLabelIntDataPreparator>()));

  label_vertical_offset_ += 2 * label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "beta     = ", &(beta),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "alpha    =", &(alpha),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "sigma0   =", &(sigma0),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "eta0     =", &(eta0),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "tau1     =", &(tau1),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  label_vertical_offset_ += label_vertical_offset_step_;

  plotLabels.push_back(std::make_shared<plotLabel>(ui->widget_plot,
                                                   label_horizontal_offset_, label_vertical_offset_, "tau2     =", &(tau2),
                                                   std::make_shared<plotLabelDoubleDataPreparator>()));

  label_vertical_offset_ = 0.75;

  AddConstantLabelToPlot("theoretical");
  plot_labels_.back()->SetColor(model_plot_pen_.color());

  AddConstantLabelToPlot("SOMKE Estimator");
  plot_labels_.back()->SetColor(kde_plot_pen_.color());


  //====================  SECOND COLUMN =================//

  label_horizontal_offset_ = 0.20;
  label_vertical_offset_ = 0.01 + 9 * label_vertical_offset_step_;

  //====================== ERRORS SUM ===================//

  label_horizontal_offset_ = 0.87;
  label_vertical_offset_ = 0.75 + label_vertical_offset_step_;

  plotLabel L2TextLabel(ui->widget_plot, label_horizontal_offset_, label_vertical_offset_, "L2   = 0");
  L2TextLabel.SetColor(kde_plot_pen_.color());

  FillDomain(&domain_, nullptr);
  for(const auto &pt : domain_){
      drawable_domain_.push_back(pt->at(0));
  }

  ui->widget_plot->replot();
  QCoreApplication::processEvents();

  double error_domain_length = 0;
  std::vector<std::vector<double>> error_domain = {};
  std::vector<double> model_values = {};
  std::vector<double> somke_values = {};

  ErrorsCalculator errors_calculator(
      &model_values, &somke_values, &error_domain, &error_domain_length
                                    );

  errors_calculators = {&errors_calculator};

  std::ifstream inFile(avg_l2_errors_file_path);
  numberOfErrorCalculations = std::count(std::istreambuf_iterator<char>(inFile),
                                         std::istreambuf_iterator<char>(), '\n');

  int drawing_start_step = 0;
  drawing_start_step = numberOfErrorCalculations * error_computation_frequency_;

  if(drawing_start_step > 0 && should_compute_errors) {
      // Open the file with average l2 errors and load them into the list.

      qDebug() << "Loading errors.";

      std::ifstream l2_errors_sums_file(l2_errors_sum_file_path);
      if(l2_errors_sums_file.is_open()) {
          l2_errors_sums.clear();
          std::string line;
          while(std::getline(l2_errors_sums_file, line)) {
              l2_errors_sums.push_back(std::stod(line));
          }
          l2_errors_sums_file.close();
      }
      else {
          log("Unable to open avg_l2_errors.txt file.");
      }
  }

  KernelPtr kernel(new SOMKENormalKernel());
  // MergingStrategyPtr merging_strategy(new SOMKEFixedMemoryMergingStrategy(max_number_of_som_seq_entries, beta));
  MergingStrategyPtr merging_strategy(new SOMKEFixedThresholdMergingStrategy(alpha, beta));
  SOMKEAlgorithm somke_algorithm(kernel, merging_strategy, neurons_number, epochs_number, data_window_size);

  for(step_number_ = 1; step_number_ <= stepsNumber; ++step_number_) {

    update_theoretical_density();

    Point stream_value = {};
    reader_->getNextRawDatum(&stream_value);

    log("Performing step: " + QString::number(step_number_));
    somke_algorithm.PerformStep(stream_value);
    log("Step performed.");

    target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

    // Error calculations
    if(step_number_ > drawing_start_step && should_compute_errors && step_number_ % error_computation_frequency_ == 0) {

      log("Getting error domain.");
      error_domain = somke_algorithm.divergence_domain_;

      // Failsafe, if divergence domain was not yet computed
      if(error_domain.empty()){

        double minX = ui->lineEdit_minX->text().toDouble();
        double maxX = ui->lineEdit_maxX->text().toDouble();
        int error_points = 1000;

        double error_domain_step = (maxX - minX) / error_points;

        error_domain = {};
        double domain_value = minX;

        while(error_domain.size() < error_points + 1){
            error_domain.push_back({domain_value});
            domain_value += error_domain_step;
        }
      }

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
      ++numberOfErrorCalculations;
      compute_errors();

      log("Saving errors.");

      // Save averaged l2 errors to the file
      std::ofstream in;
      in.open(avg_l2_errors_file_path, std::ios_base::app);
      for(size_t i = 0; i < errors_calculators.size(); ++i) {
          in << *l2_errors[i] << ";";
      }
      in << "\n";
      in.close();

      // Save the errors sum to the file
      std::ofstream in2;
      in2.open(l2_errors_sum_file_path);
      for(size_t i = 0; i < errors_calculators.size(); ++i) {
          in2 << l2_errors_sums[i] << "\n";
      }
      in2.close();

      log("Errors saved.");
    }

    // Drawing
    if(drawing_start_step < step_number_ && ( step_number_ % screen_generation_frequency_ == 0)) {
      log("Drawing in step number " + QString::number(step_number_) + ".");

      L2TextLabel.setText("L2   =" + FormatNumberForDisplay(*l2_errors[0]));

      if(!should_compute_errors){
        L2TextLabel.setText("");
      }

      DrawPlots(&somke_algorithm);

      for(const auto &label : plotLabels){
          label->updateText();
      }

      ui->widget_plot->replot();
      QCoreApplication::processEvents();

      QString imageName = dirPath + QString::number(step_number_) + ".png";
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

void MainWindow::AddErrorLabelsToPlot(const QVector<QString> &labels, const QVector<double_ptr> &values) {
  // TODO TR: I assume that labels.size() == values_references.size()

  label_vertical_offset_ = 0.75 + label_vertical_offset_step_;

  for(size_t i = 0; i < labels.size(); ++i){
    AddErrorLabelToPlot(labels[i], values[i].get());
    plot_labels_.back()->SetColor(error_indices_colors[i]);
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

void MainWindow::AddColorsLegendToPlot() {

  label_vertical_offset_ = 0.75;

  if(ui->checkBox_showEstimatedPlot->isChecked()) {
    AddConstantLabelToPlot("theoretical");
    plot_labels_.back()->SetColor(model_plot_pen_.color());
  }
  if(ui->checbox_showFullEstimator->isChecked()){
    AddConstantLabelToPlot("m=m0");
    plot_labels_.back()->SetColor(windowed_plot_pen_.color());
  }

  if(ui->checkBox_showEstimationPlot->isChecked()){
    AddConstantLabelToPlot("m variable (Sec. 3.1)");
    plot_labels_.back()->SetColor(kde_plot_pen_.color());
  }

  if(ui->checkBox_showWeightedEstimationPlot->isChecked()){
    AddConstantLabelToPlot("weights    (Sec. 3.2)");
    plot_labels_.back()->SetColor(weighted_plot_pen_.color());
  }

  if(ui->checkBox_sigmoidallyEnhancedKDE->isChecked()){
    AddConstantLabelToPlot("prediction (Sec. 3.3)");
    plot_labels_.back()->SetColor(desda_kde_plot_pen_.color());
  }

  if(ui->checkBox_REESEKDE->isChecked()){
    AddConstantLabelToPlot("atypical   (Sec. 3.4)");
    plot_labels_.back()->SetColor(desda_rare_elements_kde_plot_pen_.color());
  }

  if(ui->checkBox_kernelPrognosedPlot->isChecked()){
    AddConstantLabelToPlot("derivative");
    plot_labels_.back()->SetColor(derivative_plot_pen_.color());
  }

  if(ui->checkBox_standarizedDerivative){
    AddConstantLabelToPlot("derivative standardized");
    plot_labels_.back()->SetColor(standardized_derivative_plot_pen_.color());
  }

}

void MainWindow::AddErrorLabelToPlot(const QString &label, double *value) {
  plot_labels_.push_back(std::make_shared<plotLabel>(ui->widget_plot, label_horizontal_offset_,
                                                     label_vertical_offset_, label, value,
                                                     std::make_shared<plotLabelDoubleDataPreparator>(6)));

  label_vertical_offset_ += label_vertical_offset_step_;
}

void MainWindow::RunAccuracyExperiment() {

  auto max_seeds = 1000;
  //vector<double> samples_numbers = {5, 10, 100, 500, 1000}; // Test
  vector<double> samples_numbers = {50, 100, 1000, 5000, 10000};

  int dimensionsNumber = ui->tableWidget_dimensionKernels->rowCount();

  FillMeans(&means_);
  FillStandardDeviations(&standard_deviations_);

  vector<double> alternativeDistributionMean = {0.0};
  vector<double> alternativeDistributionStDevs = {1.0};

  vector<double> l1_sums = {0, 0, 0, 0, 0};
  vector<double> l2_sums = {0, 0, 0, 0, 0};
  vector<double> sup_sums = {0, 0, 0, 0, 0};

  // Seed = 0 is the same as Seed = 0
  for(auto seed = 1; seed <= max_seeds; ++seed){

    int samples_number_index = 0;

    srand(seed);
    ui->lineEdit_seed->setText(QString::number(seed));

    target_function_.reset(GenerateTargetFunction(&means_, &standard_deviations_));

    std::shared_ptr<kernelDensityEstimator>
        estimator(GenerateKernelDensityEstimator(dimensionsNumber));
    estimator->_shouldConsiderWeights = false;

    std::shared_ptr<distribution>
        targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));

    parser_.reset(new distributionDataParser(&attributes_data_));

    reader_.reset(
        new progressiveDistributionDataReader(targetDistribution.get(),
                                              0,
                                              0,  // Delay
                                              new normalDistribution(seed, &alternativeDistributionMean,
                                                                     &alternativeDistributionStDevs))
                 );

    reader_->gatherAttributesData(&attributes_data_);
    parser_->setAttributesOrder(reader_->getAttributesOrder());

    reservoirSamplingAlgorithm *algorithm =
        GenerateReservoirSamplingAlgorithm(reader_.get(), parser_.get());

    objects_.clear();
    stored_medoids_.clear();

    double min_val = NAN;
    double max_val = NAN;

    double mean = 0;
    QVector<double> samples = {};

    for(double samples_number : samples_numbers) {

      while(stored_medoids_.size() < samples_number) {
        auto sample_number = stored_medoids_.size() + 1;
        algorithm->performSingleStep(&objects_, sample_number);

        std::shared_ptr<cluster> newCluster =
            std::shared_ptr<cluster>(new cluster(sample_number, objects_.back()));
        newCluster->setTimestamp(sample_number);

        stored_medoids_.push_back(newCluster);

        double val = stod(objects_[sample_number - 1]->attributesValues["Val0"]);

        if(isnan(min_val)){
          min_val = max_val = val;
        }

        if(val > max_val) {
          max_val = val;
        }

        if(val < min_val) {
          min_val = val;
        }

        samples.push_back(val);
        mean += val;
      }

      double std = 0;
      mean /= samples_number;

      for(auto val : samples) {
        std = (val - mean) * (val - mean);
      }

      std /= samples_number;
      std = sqrt(std);

      pluginSmoothingParameterCounter counter(&samples, 3);
      auto h = counter.countSmoothingParameterValue();

      estimator->setSmoothingParameters({h});
      estimator->setClusters(stored_medoids_);

      double domain_length = max_val - min_val;
      std::vector<std::vector<double>> error_domain = {};
      double stepSize = domain_length / (1000);

      for(auto val = min_val; val <= max_val; val += stepSize) {
        error_domain.push_back({val});
      }

      std::vector<double> model_values = GetFunctionsValueOnDomain(target_function_.get(), error_domain);
      std::vector<double> kde_values = {};

      for(auto pt : error_domain) {
        kde_values.push_back(estimator->getValue(&pt));
      }

      std::vector<double> x = {0.0};
      qDebug() << "\tEstimator check: KDE(0) = " << estimator->getValue(&x);
      qDebug() << "\t L1: " << l1_sums;
      qDebug() << "\t L1: " << l2_sums;
      qDebug() << "\t L1: " << sup_sums;

      ErrorsCalculator errors_calculator(&model_values, &kde_values, &error_domain, &domain_length);

      l1_sums[samples_number_index] += errors_calculator.CalculateL1Error();
      l2_sums[samples_number_index] += errors_calculator.CalculateL2Error();
      sup_sums[samples_number_index] += errors_calculator.CalculateSupError();

      ++samples_number_index;
    }

    qDebug() << "Saving to file.";

    std::ofstream outfile;

    outfile.open("results.txt");

    outfile << "Seed " << seed << "\n\tl1_sums = [";

    for(auto val : l1_sums){
      outfile << val << ", ";
    }

    outfile << "]\n\tl2_sums = [";

    for(auto val : l2_sums){
      outfile << val << ", ";
    }

    outfile << "]\n\tSup = [";

    for(auto val : sup_sums){
      outfile << val << ", ";
    }

    outfile << "]";

    outfile.close();
  }

  std::ofstream outfile;

  // Average errors
  for(auto i = 0; i < samples_numbers.size(); ++i){
    l1_sums[i] /= max_seeds;
    l2_sums[i] /= max_seeds;
    sup_sums[i] /= max_seeds;
  }

  outfile.open("results.txt"); // append instead of overwrite

  outfile << "Results\n\tl1_sums = [";

  for(auto val : l1_sums){
    outfile << val << ", ";
  }

  outfile << "]\n\tl2_sums = [";

  for(auto val : l2_sums){
    outfile << val << ", ";
  }

  outfile << "]\n\tSup = [";

  for(auto val : sup_sums){
    outfile << val << ", ";
  }

  outfile << "]";

  outfile.close();

  log("Accuracy experiment finished.");
}

void MainWindow::run_3d_experiment() {

  log("Start pushed!");
  // Delay so that
  QTime dieTime= QTime::currentTime().addSecs(0);
  while (QTime::currentTime() < dieTime) {
    QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
  }

  log("3D Experiment start.");
  bool radial = false;

  //curve.setSamples( points );
  //curve.attach( contour_plot_ );

  // Initially these vectors were used in errors computation only. We now also use them for the spectrogram.
  QVector<double> error_xs = {};
  QVector<double> error_ys = {};
  QVector<double> error_zs = {};
  std::vector<double> model_function_values = {};
  std::vector<double> estimator_values = {};
  std::vector<std::vector<double>> error_domain = {};

  screen_generation_frequency_ = 10;

  // Automatic dimension update
  if(ui->spinBox_dimensionsNumber->value() == 1){
    ui->spinBox_dimensionsNumber->setValue(3);
    ui->spinBox_dimensionsNumber->editingFinished();
  }

  // Add clusters_ to the estimator
  means_ = {std::make_shared<std::vector<double >>()};
  means_.back()->push_back(0);
  means_.back()->push_back(0);
  means_.back()->push_back(0);

  standard_deviations_ = {std::make_shared<std::vector<double >>()};
  standard_deviations_.back()->push_back(1);
  standard_deviations_.back()->push_back(1);
  standard_deviations_.back()->push_back(1);

  auto densityFunction =
      new multivariateNormalProbabilityDensityFunction(means_.back().get(), standard_deviations_.back().get(), 0);

  // Create estimator object
  std::shared_ptr<kernelDensityEstimator>
      estimator(GenerateKernelDensityEstimator(3, radial));

  estimator->_shouldConsiderWeights = true;

  std::shared_ptr<distribution> targetDistribution(GenerateTargetDistribution(&means_, &standard_deviations_));
  std::vector<double> meansForDistribution = {0, 0, 0};
  std::vector<double> stDevsForDistribution = {1, 1, 1};

  parser_.reset(new distributionDataParser(&attributes_data_));

  QString expNum = "R104 (3D); Product; x move only";
  QString pc_id = "sz246";
  int drawing_start_step = 1750;
  int errors_calculation_start_step = drawing_start_step;

  bool should_compute_errors = true;

  // Prepare the reader
  reader_.reset(new progressiveDistributionDataReader(targetDistribution.get(), 0,0, new normalDistribution(0, &meansForDistribution, &stDevsForDistribution)));

  // Only to remove problems initialize the date
  QTime data_start_time(0, 0, 0); QDate data_start_date(2019, 10, 1); QDateTime data_date_time(data_start_date, data_start_time);

   // Multiple instructions in one line, for simplicity
  QString experiment_description = "assumed data stream; 3D"; QString expDesc = "assumed data stream 3D, " + pc_id;

  reader_->gatherAttributesData(&attributes_data_);
  parser_->setAttributesOrder(reader_->getAttributesOrder());

  reservoirSamplingAlgorithm *samplingAlgorithm =
      GenerateReservoirSamplingAlgorithm(reader_.get(), parser_.get());

  objects_.clear();
  clusters_ = &stored_medoids_;
  clusters_->clear();

  int pluginRank = 3;
  groupingThread gt(&stored_medoids_, parser_);

  derivative_estimator_.reset(GenerateKernelDensityEstimator(2, radial));
  enhanced_kde_.reset(GenerateKernelDensityEstimator(2, radial));

  DESDA DESDAAlgorithm(
      estimator,
      derivative_estimator_,
      enhanced_kde_,
      samplingAlgorithm,
      clusters_,
      ui->lineEdit_rarity->text().toDouble(), pluginRank
  );

  // Start the test
  step_number_ = 0;

  time_t startTime, endTime;

  double actual_l2 = 0;
  int errorCalculationsNumber = 0;
  double sum_l2 = 0.0730053 * errors_calculation_start_step / 10;

  double domain_area = 0;

  ErrorsCalculator errors_calculator(&model_function_values, &estimator_values, &error_domain, &domain_area);

  // Prepare image location.
  this->setWindowTitle("Experiment #" + expNum);
  QString driveDir = "Y:\\"; // WIT PCs after update
  //QString driveDir = "D:\\Test\\"; // Home
  //QString driveDir = "d:\\OneDrive - Instytut Badań Systemowych Polskiej Akademii Nauk\\";
  QString dirPath = driveDir + "TR Badania\\Eksperyment " + expNum + " (" + expDesc + ")\\";
  //QString dirPath = driveDir + "Eksperyment " + expNum + " (" + expDesc + ")\\";
  if(!QDir(dirPath).exists()) QDir().mkdir(dirPath);

  int steps_number = ui->lineEdit_iterationsNumber->text().toInt();

  log("Experiment started.");
  for(step_number_ = 1; step_number_ <= steps_number; ++step_number_) {

    log("New step.");
    startTime = time(nullptr);

    log("Performing new step.");
    DESDAAlgorithm.performStep();
    log("Step performed.");

    bool compute_errors = (step_number_ > errors_calculation_start_step) && should_compute_errors && (step_number_ % screen_generation_frequency_ == 0);
    bool draw_plot = (step_number_ % screen_generation_frequency_ == 0 && step_number_ >= drawing_start_step);

    log("Estimator preparation.");
    DESDAAlgorithm.prepareEstimatorForContourPlotDrawing();
    log("Estimator preparation finished.");

    // NOTE: We use error domain for spectrogram generation! That's why we compute the domain and values outside the if.
    if(compute_errors || draw_plot) {
      log("Computing domains.");
      log(compute_errors);

      error_xs = DESDAAlgorithm.getErrorDomain(0);
      error_ys = DESDAAlgorithm.getErrorDomain(1);
      error_zs = DESDAAlgorithm.getErrorDomain(2);

      // 3D error domain generation
      error_domain.clear();

      for(auto x : error_xs){
        for(auto y: error_ys){
          for(auto z: error_zs){
            error_domain.push_back({x, y, z});
          }
        }
      }

      log("Computing values of domains.");
      model_function_values = GetFunctionsValueOnDomain(densityFunction, error_domain);
      estimator_values = GetFunctionsValueOnDomain(estimator.get(), error_domain);
      log("Values computation finished.");


      log("Save data to file");
      std::ofstream file;
      std::string file_path = (dirPath + QString::number(step_number_) + ".csv").toStdString();

      std::cout << "Saving to file: " << file_path << "\n";

      file.open(file_path);

      std::string lines = "";

      for(int i = 0; i < model_function_values.size(); ++i){
        lines +=
            std::to_string(error_domain[i][0]) + ";" +
            std::to_string(error_domain[i][1]) + ";" +
            std::to_string(error_domain[i][2]) + ";" +
            std::to_string(model_function_values[i]) + ";" +
            std::to_string(estimator_values[i]) + "\n";
      }

      file << lines;
      file.close();

      log("Values saved.");
    }

    // Error calculation
    //*
    if(compute_errors) {
      log("Error calculation started.");
      //++errorCalculationsNumber;
      errorCalculationsNumber = step_number_ / screen_generation_frequency_;

      domain_area = Calculate2DDomainArea(error_domain);

      auto zLen = error_domain[error_domain.size() - 1][2] - error_domain[0][2];

      domain_area *= zLen;


      actual_l2 = errors_calculator.CalculateL2Error();
      sum_l2 += actual_l2;
      double l2_n_ = sum_l2 / errorCalculationsNumber;

      log("Save data to file");

      std::ofstream file;
      file.open((dirPath + "errors" + ".csv").toStdString(), std::ios_base::app);
      file << QString::number(l2_n_).toStdString() << "\n";
      file.close();

      log("Error calculation finished.");
    }
    //*/
    if(should_compute_errors) {
      densityFunction->setMeans(*means_.back().get());
    }

    log("Restoring weights.");
    DESDAAlgorithm.restoreClustersCWeights();

    endTime = time(nullptr);

    data_date_time = data_date_time.addSecs(3600); // Add hour to the date

    log("Step time: " + QString::number(endTime - startTime) + " s");
  }

  log("Done!");
}

void MainWindow::on_toolButton_findDataStream_clicked()
{
    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::AnyFile);
    dialog.setViewMode(QFileDialog::Detail);

    QStringList fileNames;
    if (dialog.exec()) {
        fileNames = dialog.selectedFiles();
    }

    if(fileNames.isEmpty()){
        return;
    }

    QString fileName = fileNames[0];

    this->ui->label_dataStream->setText(fileName);
}

void MainWindow::update_theoretical_density(){
    //if(step_number_ >= 1000 && step_number_ <= 4000){
    //  means_[0]->at(0) = sin(2 * 2 * 3.14 * (step_number_ - 1000) / 3000);
    //}

    //*
    if(step_number_ <= 2000){

        //log("Update!");
        //log(QString::number(means_[0]->size()));
        //log(QString::number(means_.size()));

        for(int i = 0; i < means_.size(); ++i){
            means_[i]->at(0) += 0.001;
        }
    }
    if(step_number_ > 2000 && step_number_ <= 6000){
        for(int i = 0; i < means_.size(); ++i){
            means_[i]->at(0) += 0.01;
        }
    }
    if(step_number_ == 8000){
        for(int i = 0; i < means_.size(); ++i){
            means_[i]->at(0) += 1;
        }
    }
    //*/
}

void MainWindow::set_paths(QString expNum, QString expDesc){
    //drive = "D:\\";
    drive = "Y:\\"; // WIT PCs after update

    dirPath = drive + "TR Badania\\Eksperyment " + expNum + " (" + expDesc + ")\\";
    //QString dirPath = driveDir + "Badania PK\\Eksperyment " + expNum + " (" + expDesc + ")\\";
    //QString dirPath = driveDir + "Eksperyment " + expNum + " (" + expDesc + ")\\";

    l2_errors_sum_file_path = QString(dirPath + "sum_l2_errors.txt").toStdString();
    avg_l2_errors_file_path = QString(dirPath + "avg_l2_errors.txt").toStdString();

    if(!QDir(dirPath).exists()){
        QDir().mkdir(dirPath);
        log(dirPath);
        log(QDir(dirPath).exists());
    }
}

void MainWindow::clear_errors(){
    l1_errors.clear();
    l1_errors_sums.clear();
    l2_errors.clear();
    l2_errors_sums.clear();
    sup_errors.clear();
    sup_errors_sums.clear();
    mod_errors.clear();
    mod_errors_sums.clear();

    errors_calculators.clear();
}

void MainWindow::compute_errors(){
    for(size_t i = 0; i < errors_calculators.size(); ++i){
        //l1_errors_sums[i] += errors_calculators[i]->CalculateL1Error();
        l2_errors_sums[i] += errors_calculators[i]->CalculateL2Error();
        *l2_errors[i] = l2_errors_sums[i] / numberOfErrorCalculations;
        //sup_errors_sums[i] += errors_calculators[i]->CalculateSupError();
        //mod_errors_sums[i] += errors_calculators[i]->CalculateModError();
    }
}
