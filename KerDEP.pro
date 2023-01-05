#-------------------------------------------------
#
# Project created by QtCreator 2017-02-21T16:19:57
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET      =   KerDEP
TEMPLATE    =   app
CONFIG      +=  static
CONFIG      +=  qwt

QMAKE_CXXFLAGS += -std=c++17


if(exists(k:/Libs/)){
    INCLUDEPATH += k:/Libs/Qwt-6.1.5/include/
    INCLUDEPATH += k:/Libs/Qwt-6.1.5/lib/
    INCLUDEPATH += k:/Libs/Qwt-6.1.5/
    INCLUDEPATH += k:/Libs/boost_1_75_0/
    INCLUDEPATH += k:/Libs/knnl/include/
    LIBS += -L "k:/Libs/Qwt-6.1.5/lib/" -lqwt
}

if(exists(y:/Code/)){
    INCLUDEPATH += y:/Qwt-6.1.5/include/
    INCLUDEPATH += y:/Qwt-6.1.5/lib/
    INCLUDEPATH += y:/Qwt-6.1.5/
    INCLUDEPATH += y:/boost_1_75_0/
    INCLUDEPATH += y:/knnl/include/
    LIBS += -L "y:/Qwt-6.1.5/lib/" -lqwt
}


INCLUDEPATH += $$PWD/ClusterKernelWrappers/

INCLUDEPATH += $$PWD/ClusterKernelsKDE/include/ClusterKernelsKDE/
INCLUDEPATH += $$PWD/ClusterKernelsKDE/src/

INCLUDEPATH += $$PWD/Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream
INCLUDEPATH += $$PWD/Compressed_Cumulative_WDE_Over_Stream/src

INCLUDEPATH += $$PWD/Compressed_Cumulative_WDE_Wrappers/

INCLUDEPATH += $$PWD/SOMKE/include/SOMKE/
INCLUDEPATH += $$PWD/SOMKE/src/


SOURCES     +=  main.cpp\
                Benchmarking/errorsCalculator.cpp \
                ClusterKernelWrappers/enhancedClusterKernelAlgorithm.cpp \
                ClusterKernelWrappers/epanecznikowKernelRealValuedFunction.cpp \
                ClusterKernelWrappers/univariateStreamElement.cpp \
                ClusterKernelWrappers/varianceBasedClusterKernel.cpp \
                ClusterKernelsKDE/src/ClusterKernelsAlgorithm.cpp \
                ClusterKernelsKDE/src/UnivariateListBasedClusterKernelAlgorithm.cpp \
                ClusterKernelsKDE/src/WeightedUnivariateListBasedClusterKernelAlgorithm.cpp \
                Compressed_Cumulative_WDE_Over_Stream/src/CompressedCumulativeWaveletDensityEstimator.cpp \
                Compressed_Cumulative_WDE_Over_Stream/src/TranslatedDilatedScalingFunction.cpp \
                Compressed_Cumulative_WDE_Over_Stream/src/TranslatedDilatedWaveletFunction.cpp \
                Compressed_Cumulative_WDE_Wrappers/LinearWDE.cpp \
                Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/hardThresholdingStrategy.cpp \
                Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/softThresholdingStrategy.cpp \
                Compressed_Cumulative_WDE_Wrappers/WeightedThresholdedWDE.cpp \
                Compressed_Cumulative_WDE_Wrappers/kerDepCcWde.cpp \
                Compressed_Cumulative_WDE_Wrappers/kerDepWindowedWde.cpp \
                Compressed_Cumulative_WDE_Wrappers/math_helpers.cpp \
                Compressed_Cumulative_WDE_Wrappers/weightedLinearWde.cpp \
                Reservoir_sampling/textDataReader.cpp \
                SOMKE/src/SOMKEAlgorithm.cpp \
                SOMKEWrappers/MergingStrategies/somkeFixedMemoryMergingStrategy.cpp \
                SOMKEWrappers/MergingStrategies/somkeFixedThresholdMergingStrategy.cpp \
                SOMKEWrappers/somkeNormalKernel.cpp \
                UI/QwtContourPlotUI.cpp \
                UI/plot.cpp \
                UI/plotLabelDoubleDataPreparator.cpp \
                UI/plotLabelIntDataPreparator.cpp \
                mainwindow.cpp \
                QCustomPlot/qcustomplot.cpp \
                KDE/kerneldensityestimator.cpp \
                Distributions/normaldistribution.cpp \
                Functions/Kernels/dullkernel.cpp \
                Functions/Kernels/normalkernel.cpp \
                Functions/Kernels/trianglekernel.cpp \
                Functions/Kernels/epanecznikowkernel.cpp \
                KDE/pluginsmoothingparametercounter.cpp \
                Functions/multivariatenormalprobabilitydensityfunction.cpp \
                Libraries/matrixoperationslibrary.cpp \
                Functions/complexfunction.cpp \
                Distributions/complexdistribution.cpp \
                Reservoir_sampling/biasedReservoirSamplingAlgorithm.cpp \
                Reservoir_sampling/distributiondataparser.cpp \
                Reservoir_sampling/basicReservoirSamplingAlgorithm.cpp \
                Reservoir_sampling/progressivedistributiondatareader.cpp \
                  Reservoir_sampling/sinusoidalDistributionDataReader.cpp \
                Reservoir_sampling/distributionDataSample.cpp \
                KDE/weightedSilvermanSmoothingParameterCounter.cpp \
                groupingThread/groupingThread.cpp \
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/categorical/smdCategoricalAttributesDistanceMeasure.cpp \
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/numerical/gowersNumericalAttributesDistanceMeasure.cpp \
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/numerical/smdNumericalAttributesDistanceMeasure.cpp \
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/averageLinkClusterDistanceMeasure.cpp \
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/centroidLinkClusterDistanceMeasure.cpp \
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/completeLinkClusterDistanceMeasure.cpp \
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/singleLinkClusterDistanceMeasure.cpp \
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.cpp \
                groupingThread/kMedoidsAlgorithm/categoricalAttributeData.cpp \
                groupingThread/kMedoidsAlgorithm/customObjectsDistanceMeasure.cpp \
                groupingThread/kMedoidsAlgorithm/kMedoidsAlgorithm.cpp \
                groupingThread/kMedoidsAlgorithm/numericalAttributeData.cpp \
                groupingThread/medoidStoringAlgorithm/medoidStoringAlgorithm.cpp \
                groupingThread/kMeansAlgorithm.cpp \
                DESDA.cpp \
                StationarityTests/kpssstationaritytest.cpp \
                UI/plotLabel.cpp

HEADERS     +=  mainwindow.h \
                Benchmarking/errorsCalculator.h \
                ClusterKernelWrappers/enhancedClusterKernelAlgorithm.h \
                ClusterKernelWrappers/epanecznikowKernelRealValuedFunction.h \
                ClusterKernelWrappers/univariateStreamElement.h \
                ClusterKernelWrappers/varianceBasedClusterKernel.h \
                ClusterKernelsKDE/include/ClusterKernelsKDE/ClusterKernel.h \
                ClusterKernelsKDE/include/ClusterKernelsKDE/ClusterKernelStreamElement.h \
                ClusterKernelsKDE/include/ClusterKernelsKDE/ClusterKernelsAlgorithm.h \
                ClusterKernelsKDE/include/ClusterKernelsKDE/RealValuedFunction.h \
                ClusterKernelsKDE/include/ClusterKernelsKDE/UnivariateListBasedClusterKernelAlgorithm.h \
                ClusterKernelsKDE/include/ClusterKernelsKDE/WeightedUnivariateListBasedClusterKernelAlgorithm.h \
                Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/CompressedCumulativeWaveletDensityEstimator.h \
                Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/TranslatedDilatedScalingFunction.h \
                Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/TranslatedDilatedWaveletFunction.h \
                Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/WaveletDensityEstimator.h \
                Compressed_Cumulative_WDE_Wrappers/LinearWDE.h \
                Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/ThresholdingStrategyInterface.h \
                Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/hardThresholdingStrategy.h \
                Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/softThresholdingStrategy.h \
                Compressed_Cumulative_WDE_Wrappers/WeightedThresholdedWDE.h \
                Compressed_Cumulative_WDE_Wrappers/kerDepCcWde.h \
                Compressed_Cumulative_WDE_Wrappers/kerDepWindowedWde.h \
                Compressed_Cumulative_WDE_Wrappers/math_helpers.h \
                Compressed_Cumulative_WDE_Wrappers/weightedLinearWde.h \
                QCustomPlot/qcustomplot.h \
                Functions/function.h \
                KDE/kerneldensityestimator.h \
                Distributions/distributions.h \
                Distributions/distribution.h \
                Distributions/normaldistribution.h \
                Functions/Kernels/kernel.h \
                Functions/Kernels/dullkernel.h \
                Functions/Kernels/normalkernel.h \
                Functions/Kernels/trianglekernel.h \
                Functions/Kernels/epanecznikowkernel.h \
                Functions/Kernels/kernels.h \
                KDE/pluginsmoothingparametercounter.h \
                Functions/multivariatenormalprobabilitydensityfunction.h \
                Libraries/matrixoperationslibrary.h \
                Functions/complexfunction.h \
                Distributions/complexdistribution.h \
                Reservoir_sampling/biasedReservoirSamplingAlgorithm.h \
                Reservoir_sampling/dataParser.h \
                Reservoir_sampling/dataReader.h \
                Reservoir_sampling/reservoirSamplingAlgorithm.h \
                Reservoir_sampling/distributiondataparser.h \
                Reservoir_sampling/basicReservoirSamplingAlgorithm.h \
                Reservoir_sampling/progressivedistributiondatareader.h \
                Reservoir_sampling/sinusoidalDistributionDataReader.h \
                Reservoir_sampling/distributionDataSample.h \
                KDE/smoothingParameterCounter.h \
                KDE/weightedSilvermanSmoothingParameterCounter.h \
                Reservoir_sampling/sample.h \
                Reservoir_sampling/textDataReader.h \
                SOMKE/include/SOMKE/Kernel.h \
                SOMKE/include/SOMKE/SOMKEAlgorithm.h \
                SOMKE/include/SOMKE/SOMKEMergingStrategy.h \
                SOMKE/include/SOMKE/SOMSequenceEntry.h \
                SOMKE/include/SOMKE/wtm_localized_training_algorithm.h \
                SOMKEWrappers/MergingStrategies/somkeFixedMemoryMergingStrategy.h \
                SOMKEWrappers/MergingStrategies/somkeFixedThresholdMergingStrategy.h \
                SOMKEWrappers/somkeNormalKernel.h \
                UI/QwtContourPlotUI.h \
                UI/i_plotLabelDataPreparator.h \
                UI/plot.h \
                UI/plotLabelDoubleDataPreparator.h \
                UI/plotLabelIntDataPreparator.h \
                groupingThread/groupingThread.h \
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/categorical/smdCategoricalAttributesDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/numerical/gowersNumericalAttributesDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/numerical/smdNumericalAttributesDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/averageLinkClusterDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/centroidLinkClusterDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/completeLinkClusterDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/singleLinkClusterDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/dataParsers/textDataParser.h \
                groupingThread/kMedoidsAlgorithm/dataReaders/textDataReader.h \
                groupingThread/kMedoidsAlgorithm/dataSamples/textDataSample.h \
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h \
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/distanceBasedGroupingAlgorithm.h \
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/groupingAlgorithm.h \
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/sample.h \
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/similarityBasedGroupingAlgorithm.h \
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/summarizedCluster.h \
                groupingThread/kMedoidsAlgorithm/attributeData.h \
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/categoricalAttributeData.h \
                groupingThread/kMedoidsAlgorithm/clustersDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/customObjectsDistanceMeasure.h \
                groupingThread/kMedoidsAlgorithm/dataParser.h \
                groupingThread/kMedoidsAlgorithm/dataReader.h \
                groupingThread/kMedoidsAlgorithm/kMedoidsAlgorithm.h \
                groupingThread/kMedoidsAlgorithm/numericalAttributeData.h \
                groupingThread/kMedoidsAlgorithm/objectsDistanceMeasure.h \
                groupingThread/medoidStoringAlgorithm/medoidStoringAlgorithm.h \
                groupingThread/kMeansAlgorithm.h \
                DESDA.h \
                StationarityTests/kpssstationaritytest.h \
                StationarityTests/i_stationaritytest.h \
                UI/plotLabel.h

FORMS       +=  mainwindow.ui

DISTFILES   +=  .gitignore \
  .gitmodules \
  .idea/.gitignore \
  .idea/KerDEP.iml \
  .idea/dataSources.local.xml \
  .idea/misc.xml \
  .idea/modules.xml \
  .idea/vcs.xml \
  .idea/workspace.xml \
  CMakeLists.txt \
  ClusterKernelsKDE/.gitignore \
  ClusterKernelsKDE/CMakeLists.txt \
  Compressed_Cumulative_WDE_Over_Stream/.gitignore \
  Compressed_Cumulative_WDE_Over_Stream/CMakeLists.txt \
  Compressed_Cumulative_WDE_Over_Stream/FindBoost.cmake \
  FindQwt.cmake \
  SOMKE/.gitignore \
  SOMKE/CMakeLists.txt \
  groupingThread/kMedoidsAlgorithm/.gitignore \
  groupingThread/medoidStoringAlgorithm/.gitignore
