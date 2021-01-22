cmake_minimum_required(VERSION 3.5)

project(KerDEP LANGUAGES CXX)
set(CMAKE_PREFIX_PATH "D:\\Qt\\5.15.0\\mingw81_64\\")
set(CMAKE_PREFIX_PATH "D:\\Qwt-6.1.5\\")
set(CMAKE_PREFIX_PATH "D:\\Qwt-6.1.5\\lib\\")
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++ -lwsock32 -lws2_32 ${CMAKE_CSS_STANDARD_LIBRARIES}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive")

set(BOOST_ROOT d:/boost_1_75_0/)
FIND_PACKAGE(Boost)

IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

include_directories(include/Compressed_Cumulative_WDE_Over_Stream/)

add_executable(Compressed_Cumulative_WDE_Over_Stream include/Compressed_Cumulative_WDE_Over_Stream/CompressedCumulativeWaveletDensityEstimator.h src/CompressedCumulativeWaveletDensityEstimator.cpp include/Compressed_Cumulative_WDE_Over_Stream/TranslatedDilatedScalingFunction.h src/TranslatedDilatedScalingFunction.cpp include/Compressed_Cumulative_WDE_Over_Stream/WaveletDensityEstimator.h ../Compressed_Cumulative_WDE_Wrappers/LinearWDE.h ../Compressed_Cumulative_WDE_Wrappers/LinearWDE.cpp)

#if(NOT Boost_FOUND)
#    message(FATAL_ERROR "Could not find boost!")
#endif()
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
message(STATUS "Boost_VERSION: ${Boost_VERSION}")

find_package(QT NAMES Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets Gui Core printsupport REQUIRED)
find_package(Qwt REQUIRED)
find_package(Qt5PrintSupport REQUIRED)

add_executable(KerDEP
                main.cpp
                DESDAReservoir.cpp
                UI/QwtContourPlotUI.cpp
                UI/plot.cpp
                UI/plotLabelDoubleDataPreparator.cpp
                UI/plotLabelIntDataPreparator.cpp
                mainwindow.cpp
                QCustomPlot/qcustomplot.cpp
                KDE/kerneldensityestimator.cpp
                Distributions/normaldistribution.cpp
                Functions/Kernels/dullkernel.cpp
                Functions/Kernels/normalkernel.cpp
                Functions/Kernels/trianglekernel.cpp
                Functions/Kernels/epanecznikowkernel.cpp
                KDE/pluginsmoothingparametercounter.cpp
                Functions/multivariatenormalprobabilitydensityfunction.cpp
                Libraries/matrixoperationslibrary.cpp
                Functions/complexfunction.cpp
                Distributions/complexdistribution.cpp
                Reservoir_sampling/biasedReservoirSamplingAlgorithm.cpp
                Reservoir_sampling/distributiondataparser.cpp
                Reservoir_sampling/basicReservoirSamplingAlgorithm.cpp
                Reservoir_sampling/progressivedistributiondatareader.cpp
                Reservoir_sampling/distributionDataSample.cpp
                KDE/weightedSilvermanSmoothingParameterCounter.cpp
                groupingThread/groupingThread.cpp
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/categorical/smdCategoricalAttributesDistanceMeasure.cpp
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/numerical/gowersNumericalAttributesDistanceMeasure.cpp
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/numerical/smdNumericalAttributesDistanceMeasure.cpp
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/averageLinkClusterDistanceMeasure.cpp
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/centroidLinkClusterDistanceMeasure.cpp
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/completeLinkClusterDistanceMeasure.cpp
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/singleLinkClusterDistanceMeasure.cpp
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.cpp
                groupingThread/kMedoidsAlgorithm/categoricalAttributeData.cpp
                groupingThread/kMedoidsAlgorithm/customObjectsDistanceMeasure.cpp
                groupingThread/kMedoidsAlgorithm/kMedoidsAlgorithm.cpp
                groupingThread/kMedoidsAlgorithm/numericalAttributeData.cpp
                groupingThread/medoidStoringAlgorithm/medoidStoringAlgorithm.cpp
                groupingThread/kMeansAlgorithm.cpp
                DESDA.cpp
                StationarityTests/kpssstationaritytest.cpp
                UI/plotLabel.cpp
                mainwindow.h
                DESDAReservoir.h
                QCustomPlot/qcustomplot.h
                Functions/function.h
                KDE/kerneldensityestimator.h
                Distributions/distributions.h
                Distributions/distribution.h
                Distributions/normaldistribution.h
                Functions/Kernels/kernel.h
                Functions/Kernels/dullkernel.h
                Functions/Kernels/normalkernel.h
                Functions/Kernels/trianglekernel.h
                Functions/Kernels/epanecznikowkernel.h
                Functions/Kernels/kernels.h
                KDE/pluginsmoothingparametercounter.h
                Functions/multivariatenormalprobabilitydensityfunction.h
                Libraries/matrixoperationslibrary.h
                Functions/complexfunction.h
                Distributions/complexdistribution.h
                Reservoir_sampling/biasedReservoirSamplingAlgorithm.h
                Reservoir_sampling/dataParser.h
                Reservoir_sampling/dataReader.h
                Reservoir_sampling/reservoirSamplingAlgorithm.h
                Reservoir_sampling/distributiondataparser.h
                Reservoir_sampling/basicReservoirSamplingAlgorithm.h
                Reservoir_sampling/progressivedistributiondatareader.h
                Reservoir_sampling/distributionDataSample.h
                KDE/smoothingParameterCounter.h
                KDE/weightedSilvermanSmoothingParameterCounter.h
                UI/QwtContourPlotUI.h
                UI/i_plotLabelDataPreparator.h
                UI/plot.h
                UI/plotLabelDoubleDataPreparator.h
                UI/plotLabelIntDataPreparator.h
                groupingThread/groupingThread.h
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/categorical/smdCategoricalAttributesDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/numerical/gowersNumericalAttributesDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasures/numerical/smdNumericalAttributesDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/averageLinkClusterDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/centroidLinkClusterDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/completeLinkClusterDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/clusterDistanceMeasures/singleLinkClusterDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/distanceBasedGroupingAlgorithm.h
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/groupingAlgorithm.h
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/sample.h
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/similarityBasedGroupingAlgorithm.h
                groupingThread/kMedoidsAlgorithm/groupingAlgorithm/summarizedCluster.h
                groupingThread/kMedoidsAlgorithm/attributeData.h
                groupingThread/kMedoidsAlgorithm/attributesDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/categoricalAttributeData.h
                groupingThread/kMedoidsAlgorithm/clustersDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/customObjectsDistanceMeasure.h
                groupingThread/kMedoidsAlgorithm/dataParser.h
                groupingThread/kMedoidsAlgorithm/dataReader.h
                groupingThread/kMedoidsAlgorithm/kMedoidsAlgorithm.h
                groupingThread/kMedoidsAlgorithm/numericalAttributeData.h
                groupingThread/kMedoidsAlgorithm/objectsDistanceMeasure.h
                groupingThread/medoidStoringAlgorithm/medoidStoringAlgorithm.h
                groupingThread/kMeansAlgorithm.h
                DESDA.h
                StationarityTests/kpssstationaritytest.h
                StationarityTests/i_stationaritytest.h
                UI/plotLabel.h
                mainwindow.ui
                Benchmarking/errorsCalculator.cpp
                Benchmarking/errorsCalculator.h
        ClusterKernelWrappers/epanecznikowKernelRealValuedFunction.cpp ClusterKernelWrappers/epanecznikowKernelRealValuedFunction.h ClusterKernelWrappers/varianceBasedClusterKernel.cpp ClusterKernelWrappers/varianceBasedClusterKernel.h ClusterKernelWrappers/enhancedClusterKernelAlgorithm.cpp ClusterKernelWrappers/enhancedClusterKernelAlgorithm.h ClusterKernelWrappers/univariateStreamElement.cpp ClusterKernelWrappers/univariateStreamElement.h
        ClusterKernelsKDE/src/ClusterKernelsAlgorithm.cpp
        ClusterKernelsKDE/src/UnivariateListBasedClusterKernelAlgorithm.cpp ClusterKernelsKDE/src/WeightedUnivariateListBasedClusterKernelAlgorithm.cpp ClusterKernelsKDE/include/ClusterKernelsKDE/WeightedUnivariateListBasedClusterKernelAlgorithm.h
        Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/CompressedCumulativeWaveletDensityEstimator.h
        Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/TranslatedDilatedScalingFunction.h
        Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/WaveletDensityEstimator.h
        Compressed_Cumulative_WDE_Over_Stream/src/CompressedCumulativeWaveletDensityEstimator.cpp
        Compressed_Cumulative_WDE_Over_Stream/src/TranslatedDilatedScalingFunction.cpp
        Compressed_Cumulative_WDE_Wrappers/LinearWDE.h
        Compressed_Cumulative_WDE_Wrappers/LinearWDE.cpp
        )

include_directories(${QWT_INCLUDE_DIR})

include_directories(ClusterKernelsKDE/include/ClusterKernelsKDE/)
include_directories(ClusterKernelsKDE/src/)
include_directories(ClusterKernelsWrapper/)

include_directories(Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/)
include_directories(Compressed_Cumulative_WDE_Over_Stream/src/)
include_directories(Compressed_Cumulative_WDE_Wrappers/)

target_link_libraries(KerDEP PRIVATE
        Qt${QT_VERSION_MAJOR}::Widgets
        ${QWT_LIBRARY}
        Qt${QT_VERSION_MAJOR}::PrintSupport
        Qt${QT_VERSION_MAJOR}::Core
        Qt${QT_VERSION_MAJOR}::Gui
        )