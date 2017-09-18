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

QMAKE_CXXFLAGS += -std=c++11

SOURCES     +=  main.cpp\
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
    Detectors/rareElementsDetector.cpp

HEADERS     +=  mainwindow.h \
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
    Detectors/rareElementsDetector.h

FORMS       +=  mainwindow.ui

DISTFILES   +=  .gitignore
