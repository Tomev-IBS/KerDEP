#-------------------------------------------------
#
# Project created by QtCreator 2017-02-21T16:19:57
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = KerDEP
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    QCustomPlot/qcustomplot.cpp \
    Functions/gaussianprobabilitydensityfunction.cpp

HEADERS  += mainwindow.h \
    QCustomPlot/qcustomplot.h \
    customfunction.h \
    Functions/gaussianprobabilitydensityfunction.h \
    Functions/function.h

FORMS    += mainwindow.ui

DISTFILES += \
    .gitignore
