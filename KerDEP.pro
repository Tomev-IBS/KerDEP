#-------------------------------------------------
#
# Project created by QtCreator 2017-02-21T16:19:57
#
#-------------------------------------------------

QT       += core gui
QT       += datavisualization

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET      =   KerDEP
TEMPLATE    =   app
CONFIG      +=  static

QMAKE_CXXFLAGS += -std=c++11

SOURCES     +=  main.cpp\
    DESDAReservoir.cpp \
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
    DESDAReservoir.h \
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
                Reservoir_sampling/distributionDataSample.h \
                KDE/smoothingParameterCounter.h \
                KDE/weightedSilvermanSmoothingParameterCounter.h \
  UI/i_plotLabelDataPreparator.h \
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

DISTFILES   +=  .gitignore
