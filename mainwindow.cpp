#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "random"
#include "QDebug"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Set validators for lineedits
    QIntValidator* lineEditsValidator = new QIntValidator(1, 1000000, this);

    ui->lineEdit_sampleSize->setValidator(lineEditsValidator);
    ui->lineEdit_seed->setValidator(lineEditsValidator);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_generate_clicked()
{
    // Log that application started generating plot
    qDebug() << "Plot generation started.";
    qDebug() << "Seed: " + ui->lineEdit_seed->text() +
                ", Sample size: " + ui->lineEdit_sampleSize->text();

    // Generate a vector of values from normal distribution
}
