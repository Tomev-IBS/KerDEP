#-------------------------------------------------
#
# Project created by QtCreator 2017-02-21T16:19:57
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = KerDEP
TEMPLATE = app
CONFIG += static

SOURCES += main.cpp\
        mainwindow.cpp \
    QCustomPlot/qcustomplot.cpp \
    Functions/gaussianprobabilitydensityfunction.cpp \
    KDE/kerneldensityestimator.cpp \
    Functions/trianglekernelfunction.cpp \
    Functions/epanecznikowkernelfunction.cpp \
    Functions/dullkernelfunction.cpp \
    Distributions/normaldistribution.cpp

HEADERS  += mainwindow.h \
    QCustomPlot/qcustomplot.h \
    Functions/gaussianprobabilitydensityfunction.h \
    Functions/function.h \
    KDE/kerneldensityestimator.h \
    Functions/trianglekernelfunction.h \
    Functions/functions.h \
    Functions/epanecznikowkernelfunction.h \
    Functions/dullkernelfunction.h \
    Distributions/distributions.h \
    Distributions/distribution.h \
    Distributions/normaldistribution.h

FORMS    += mainwindow.ui

DISTFILES += \
    .gitignore \
