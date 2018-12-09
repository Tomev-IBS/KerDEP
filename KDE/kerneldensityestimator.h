#ifndef KERNELDENSITYESTIMATOR_H
#define KERNELDENSITYESTIMATOR_H

#include "../Functions/Kernels/kernels.h"
#include "../Distributions/distribution.h"

#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"

#include <QObject>

#include <memory>

enum estimatorsKernelsType
{
    PRODUCT =   0,
    RADIAL  =   1
};

class kernelDensityEstimator
{
    public:
        kernelDensityEstimator(QVector<std::shared_ptr<QVector<qreal>>>* samples, QVector<qreal>* smoothingParameter, QVector<QString>* carriersRestrictions, int kernelType, QVector<int>* kernelsIDs);

        void setSamples(QVector<std::shared_ptr<QVector<qreal>>>* samples);
        int setClusters(std::vector<std::shared_ptr<cluster> > clusters);

        int setSmoothingParameters(std::vector<double> smoothingParams);

        void setAdditionalMultipliers(std::vector<double> multipliers);

        qreal getValue(QVector<qreal>* x);


    private:

        int kernelType;
        double weight;

        std::vector<std::vector<double>> _spModifyingParameters;
        double _modificationIntensivity = 0.5;

        void updateSPModifyingParameters();
        std::vector<double> getSParameter();


        QVector<std::shared_ptr<QVector<qreal>>>   samples;
        QVector<kernelPtr>                         kernels;
        QVector<qreal>                             smoothingParameters;
        QVector<QString>                           carriersRestrictions;
        std::vector<double>                        additionalMultipliers;

        std::vector<std::shared_ptr<cluster>> clusters;

        qreal getProductKernelValue(QVector<qreal>* x);
          double getProductValuesFromClusters(QVector<qreal> *x);
            int extractSampleFromCluster(std::shared_ptr<cluster> c, QVector<qreal> *smpl);
            double getProductKernelAddendFromSample(QVector<qreal> *sample, QVector<qreal> *x);
            double getProductKernelAddendFromClusterIndex(int i, QVector<qreal> *x);
          double getProductValuesFromSamples(QVector<qreal> *x);
        qreal getRadialKernelValue(QVector<qreal>* x);

        void fillKernelsList(QVector<int>* kernelsIDs);
            void addProductKernelsToTheList(QVector<int>* kernelsIDs);
            void addRadialKernelsToTheList(QVector<int>* kernelsIDs);

        int partitionCharacteristicFunction(qreal carrier, qreal restriction);
};

#endif // KERNELDENSITYESTIMATOR_H
