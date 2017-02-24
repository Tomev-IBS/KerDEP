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
    Functions/trianglekernelfunction.cpp

HEADERS  += mainwindow.h \
    QCustomPlot/qcustomplot.h \
    Functions/gaussianprobabilitydensityfunction.h \
    Functions/function.h \
    KDE/kerneldensityestimator.h \
    Functions/trianglekernelfunction.h \
    Functions/functions.h

FORMS    += mainwindow.ui

DISTFILES += \
    .gitignore \
