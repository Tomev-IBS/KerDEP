#ifndef KERNELDENSITYESTIMATOR_H
#define KERNELDENSITYESTIMATOR_H

#include "../Functions/Kernels/kernels.h"
#include "../Distributions/distribution.h"

#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"

#include <memory>

using std::string;

enum estimatorsKernelsType
{
    PRODUCT =   0
};

class kernelDensityEstimator: public function
{
    public:
        kernelDensityEstimator(vector<std::shared_ptr<vector<double>>>* samples, vector<double>* smoothingParameter, vector<string>* carriersRestrictions, int kernelType, vector<int>* kernelsIDs);

        void setSamples(vector<std::shared_ptr<vector<double>>>* samples);
        unsigned long long setClusters(std::vector<std::shared_ptr<cluster> > clusters);

        int setSmoothingParameters(const std::vector<double> &smoothingParams);

        void setAdditionalMultipliers(std::vector<double> multipliers);

        double getValue(vector<double>* x);

        bool _shouldConsiderWeights = true;

        int getDimension();

        void updateSPModifyingParameters();

        // DEBUG ONLY
        bool _printClusters = false;
        //

    private:

        int kernelType;
        double weight;

        std::vector<std::vector<double>> _spModifyingParameters;
        double _modificationIntensivity = 0.5;


        std::vector<double> getSParameter();


        vector<std::shared_ptr<vector<double>>>   samples;
        vector<kernelPtr>                         kernels;
        vector<double>                            smoothingParameters;
        vector<string>                            carriersRestrictions;
        std::vector<double>                       additionalMultipliers;

        std::vector<std::shared_ptr<cluster>> clusters;

        double getProductKernelValue(vector<double>* x);
          double getProductValuesFromClusters(vector<double> *x);
            int extractSampleFromCluster(std::shared_ptr<cluster> c, vector<double> *smpl);
            double getProductKernelAddendFromSample(vector<double> *sample, vector<double> *x);
            double getProductKernelAddendFromClusterIndex(int i, vector<double> *x);
          double getProductValuesFromSamples(vector<double> *x);

        void fillKernelsList(vector<int>* kernelsIDs);
            void addProductKernelsToTheList(vector<int>* kernelsIDs);

        int partitionCharacteristicFunction(double carrier, double restriction);
};

#endif // KERNELDENSITYESTIMATOR_H
