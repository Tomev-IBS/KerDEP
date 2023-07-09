#ifndef KERNELDENSITYESTIMATOR_H
#define KERNELDENSITYESTIMATOR_H

#include "../Functions/Kernels/kernels.h"
#include "../Distributions/distribution.h"

#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"

#include <memory>
#include <armadillo>

using std::string;
using std::vector;
using namespace arma;

enum estimatorsKernelsType {
  PRODUCT = 0
};

class kernelDensityEstimator : public function {
  public:
    kernelDensityEstimator(vector<std::shared_ptr<vector<double>>> *samples, vector<double> *smoothingParameter,
                           vector<string> *carriersRestrictions, int kernelType, vector<int> *kernelsIDs,
                           const bool &radial = false);
    void setSamples(vector<std::shared_ptr<vector<double>>> *samples);
    unsigned long long setClusters(vector<std::shared_ptr<cluster> > clusters);
    int setSmoothingParameters(const vector<double> &smoothingParams);
    void setAdditionalMultipliers(vector<double> multipliers);
    double getValue(vector<double> *x);
    bool _shouldConsiderWeights = true;
    int getDimension();
    void updateSPModifyingParameters();

    // DEBUG ONLY
    bool _printClusters = false;

    // Radial kernel
    bool _radial;
    bool _full_transform = false;

    void updateCovarianceMatrix();
    double getRadialKernelValue(vector<double>* x) const;

  protected:

    int kernelType;
    double weight;
    vector<vector<double>> _spModifyingParameters;
    double _modificationIntensivity = 0.5;

    vector<double> getSParameter();
    vector<std::shared_ptr<vector<double>>> samples;
    vector<kernelPtr> kernels;
    vector<double> smoothingParameters;
    vector<string> carriersRestrictions;
    vector<double> additionalMultipliers;
    vector<std::shared_ptr<cluster>> clusters;
    double getProductKernelValue(vector<double> *x);
    double getProductValuesFromClusters(vector<double> *x);
    //int extractSampleFromCluster(std::shared_ptr<cluster> c, vector<double> *smpl);
    [[nodiscard]] vector<double> extractSampleFromCluster(const std::shared_ptr<cluster> &c) const;
    double getProductKernelAddendFromSample(vector<double> *sample, vector<double> *x);
    double getProductKernelAddendFromClusterIndex(int i, vector<double> *x);
    double getProductValuesFromSamples(vector<double> *x);
    void fillKernelsList(vector<int> *kernelsIDs);
    void addProductKernelsToTheList(vector<int> *kernelsIDs);
    int partitionCharacteristicFunction(double carrier, double restriction);

    // Radial kernel
    mat _covarianceMatrix;
    mat cov_inv;

    void compute_covariance_matrix();
    [[nodiscard]] double compute_weighted_covariance(const int &i, const int &j, const vector<vector<double>> &data) const;
    [[nodiscard]] double compute_weighted_mean(const int &i, const vector<vector<double>> &data) const;
};

#endif // KERNELDENSITYESTIMATOR_H
