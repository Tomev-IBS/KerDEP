cmake_minimum_required(VERSION 3.6)

project(DEDSTA LANGUAGES CXX)
set(CMAKE_PREFIX_PATH "D:/Qwt-6.1.5/lib/")
list(PREPEND CMAKE_PREFIX_PATH "k:/Libs/Qt/5.15.2/mingw81_64/lib/\;")
message("Printing CMAKE_PREFIX_PATH")
message(${CMAKE_PREFIX_PATH})
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++ -lwsock32 -lws2_32 ${CMAKE_CSS_STANDARD_LIBRARIES}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive")

IF(EXISTS d:/boost_1_75_0/)
    set(BOOST_ROOT d:/boost_1_75_0/)
ENDIF()
IF(EXISTS k:/libs/boost_1_75_0/)
    set(BOOST_ROOT k:/libs/boost_1_75_0/)
ENDIF()

FIND_PACKAGE(Boost)

IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

include_directories(include/Compressed_Cumulative_WDE_Over_Stream/)

#if(NOT Boost_FOUND)
#    message(FATAL_ERROR "Could not find boost!")
#endif()
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
message(STATUS "Boost_VERSION: ${Boost_VERSION}")

IF(EXISTS k:/Libs/Qt/5.15.2/mingw81_64/lib/cmake)
    # set(Qt5_ROOT k:/Libs/Qt/5.15.2/mingw81_64/lib/cmake)
    # set(CMAKE_MODULE_PATH k:/Libs/Qt/5.15.2/mingw81_64/lib/cmake/)
    list(PREPEND CMAKE_MODULE_PATH "k:\\Libs\\Qt\\5.15.2\\mingw81_64\\lib\\cmake\\\;")
    list(PREPEND CMAKE_MODULE_PATH "k:\\Libs\\Qt\\5.15.2\\mingw81_64\\lib\\cmake\\Qt5\\\;")
    list(PREPEND CMAKE_MODULE_PATH "k:\\Libs\\Qt\\5.15.2\\mingw81_64\\lib\\\;")
    set(Qt5_DIR k:/Libs/Qt/5.15.2/mingw81_64/lib/cmake)
ENDIF()

message("Printing CMAKE_MODULE_PATH")
message(${CMAKE_MODULE_PATH})
message("Done")

find_package(Qt5 COMPONENTS Widgets Gui Core printsupport REQUIRED)
# find_package(Qt${QT_VERSION_MAJOR} COMPONENTS Widgets Gui Core printsupport REQUIRED PATH_SUFFIXES "k:/Libs/Qt/mingw81_64/lib/cmake/")

# Nothing works, CMAKE sucks.
# set(Qt5_CONFIG "k:/Libs/Qt/5.15.2/mingw81_64/lib/cmake/Qt5/Qt5Config.cmake")
set(Qt5Core_CONFIG "k:/Libs/Qt/5.15.2/mingw81_64/lib/cmake/Qt5Core/Qt5CoreConfig.cmake")
set(Qt5Widgets_CONFIG "k:/Libs/Qt/5.15.2/mingw81_64/lib/cmake/Qt5Widgets/Qt5WidgetsConfig.cmake")
set(Qt5Gui_CONFIG "k:/Libs/Qt/5.15.2/mingw81_64/lib/cmake/Qt5Gui/Qt5GuiConfig.cmake")
set(Qt5printsupport_CONFIG "k:/Libs/Qt/5.15.2/mingw81_64/lib/cmake/Qt5PrintSupport/Qt5PrintSupportConfig.cmake")

message(" ------------- QT STATUS -----------------")
message("   -- Qt${QT_VERSION_MAJOR} PATH             = " ${Qt5_CONFIG})
message("   -- Qt${QT_VERSION_MAJOR}Core PATH         = " ${Qt5Core_CONFIG})
message("   -- Qt${QT_VERSION_MAJOR}Widgets PATH      = " ${Qt5Widgets_CONFIG})
message("   -- Qt${QT_VERSION_MAJOR}Gui PATH          = " ${Qt5Gui_CONFIG})
message("   -- Qt${QT_VERSION_MAJOR}PrintSupport PATH = " ${Qt5printsupport_CONFIG})
message(" ------------- QT STATUS END --------------")

find_package(Qwt REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(DEDSTA
        main.cpp
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
        DEDSTA.cpp
        StationarityTests/kpssstationaritytest.cpp
        UI/plotLabel.cpp
        mainwindow.h
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
        DEDSTA.h
        StationarityTests/kpssstationaritytest.h
        StationarityTests/i_stationaritytest.h
        UI/plotLabel.h
        mainwindow.ui
        Benchmarking/errorsCalculator.cpp
        Benchmarking/errorsCalculator.h
        ClusterKernelWrappers/epanecznikowKernelRealValuedFunction.cpp
        ClusterKernelWrappers/epanecznikowKernelRealValuedFunction.h
        ClusterKernelWrappers/varianceBasedClusterKernel.cpp
        ClusterKernelWrappers/varianceBasedClusterKernel.h
        ClusterKernelWrappers/enhancedClusterKernelAlgorithm.cpp
        ClusterKernelWrappers/enhancedClusterKernelAlgorithm.h
        ClusterKernelWrappers/univariateStreamElement.cpp
        ClusterKernelWrappers/univariateStreamElement.h
        ClusterKernelsKDE/src/ClusterKernelsAlgorithm.cpp
        ClusterKernelsKDE/src/UnivariateListBasedClusterKernelAlgorithm.cpp
        ClusterKernelsKDE/src/WeightedUnivariateListBasedClusterKernelAlgorithm.cpp
        ClusterKernelsKDE/include/ClusterKernelsKDE/WeightedUnivariateListBasedClusterKernelAlgorithm.h
        Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/CompressedCumulativeWaveletDensityEstimator.h
        Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/TranslatedDilatedScalingFunction.h
        Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/TranslatedDilatedWaveletFunction.h
        Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/WaveletDensityEstimator.h
        Compressed_Cumulative_WDE_Over_Stream/src/CompressedCumulativeWaveletDensityEstimator.cpp
        Compressed_Cumulative_WDE_Over_Stream/src/TranslatedDilatedScalingFunction.cpp
        Compressed_Cumulative_WDE_Over_Stream/src/TranslatedDilatedWaveletFunction.cpp
        Compressed_Cumulative_WDE_Wrappers/LinearWDE.h
        Compressed_Cumulative_WDE_Wrappers/LinearWDE.cpp
        Compressed_Cumulative_WDE_Wrappers/DEDSTACcWde.cpp
        Compressed_Cumulative_WDE_Wrappers/DEDSTACcWde.h
        Compressed_Cumulative_WDE_Wrappers/math_helpers.h
        Compressed_Cumulative_WDE_Wrappers/math_helpers.cpp
        Compressed_Cumulative_WDE_Wrappers/DEDSTAWindowedWde.cpp
        Compressed_Cumulative_WDE_Wrappers/DEDSTAWindowedWde.h
        Compressed_Cumulative_WDE_Wrappers/weightedLinearWde.cpp
        Compressed_Cumulative_WDE_Wrappers/weightedLinearWde.h
        Compressed_Cumulative_WDE_Wrappers/WeightedThresholdedWDE.cpp
        Compressed_Cumulative_WDE_Wrappers/WeightedThresholdedWDE.h
        Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/ThresholdingStrategyInterface.h
        Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/hardThresholdingStrategy.cpp
        Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/hardThresholdingStrategy.h
        Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/softThresholdingStrategy.cpp
        Compressed_Cumulative_WDE_Wrappers/ThresholdingStrategies/softThresholdingStrategy.h
        SOMKE/include/SOMKE/SOMKEAlgorithm.h
        SOMKE/include/SOMKE/SOMSequenceEntry.h
        SOMKE/include/SOMKE/Kernel.h
        SOMKE/src/SOMKEAlgorithm.cpp
        SOMKEWrappers/somkeNormalKernel.cpp SOMKEWrappers/somkeNormalKernel.h SOMKE/include/SOMKE/SOMKEMergingStrategy.h SOMKEWrappers/MergingStrategies/somkeFixedMemoryMergingStrategy.cpp SOMKEWrappers/MergingStrategies/somkeFixedMemoryMergingStrategy.h SOMKEWrappers/MergingStrategies/somkeFixedThresholdMergingStrategy.cpp SOMKEWrappers/MergingStrategies/somkeFixedThresholdMergingStrategy.h Reservoir_sampling/textDataReader.cpp Reservoir_sampling/textDataReader.h Reservoir_sampling/sinusoidalDistributionDataReader.cpp Reservoir_sampling/sinusoidalDistributionDataReader.h Distributions/alternatingSplittingDistribution.cpp Distributions/alternatingSplittingDistribution.h DataBinning/linearBinning.cpp DataBinning/linearBinning.h KDE/WeightedCVBandwidthSelector.cpp KDE/WeightedCVBandwidthSelector.h)

include_directories(${QWT_INCLUDE_DIR})

IF (EXISTS "k:\\Libs\\knnl\\include\\")
    set(knnl_include "k:\\Libs\\knnl\\include\\")
ENDIF()

IF (EXISTS "k:\\Libs\\armadillo-11.4.3\\include\\")
    set(armadillo_include "k:\\Libs\\armadillo-11.4.3\\include\\")
ENDIF()

find_library(openblas_lib libopenblas HINTS "k:\\Libs\\OpenBLAS_0.3.21\\lib\\")
message(STATUS ${openblas_lib})

include_directories(ClusterKernelsKDE/include/ClusterKernelsKDE/)
include_directories(ClusterKernelsKDE/src/)
include_directories(ClusterKernelsWrapper/)

include_directories(Compressed_Cumulative_WDE_Over_Stream/include/Compressed_Cumulative_WDE_Over_Stream/)
include_directories(Compressed_Cumulative_WDE_Over_Stream/src/)
include_directories(Compressed_Cumulative_WDE_Wrappers/)

include_directories(SOMKE/include/SOMKE/)
include_directories(SOMKE/src/)

target_include_directories(DEDSTA PUBLIC ${knnl_include})
target_include_directories(DEDSTA PUBLIC ${armadillo_include})

target_link_libraries(DEDSTA PRIVATE
        Qt5::Widgets
        Qt5::PrintSupport
        Qt5::Core
        Qt5::Gui
        ${QWT_LIBRARY}
        ${openblas_lib}
)



