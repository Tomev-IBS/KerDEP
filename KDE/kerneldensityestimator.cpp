#include "kerneldensityestimator.h"

#include <cmath>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


kernelDensityEstimator::kernelDensityEstimator(
    vector<std::shared_ptr<vector<double>>>* samples,
    vector<double>* smoothingParameters,
    vector<string> *carriersRestrictions,
    int kernelType,
    vector<int>* kernelsIDs)
    : kernelType(kernelType){

    this->samples               = vector<std::shared_ptr<vector<double>>>(*samples);
    this->smoothingParameters   = vector<double>(*smoothingParameters);
    this->carriersRestrictions  = vector<string>(*carriersRestrictions);

    if(kernelsIDs->size() != smoothingParameters->size())
    {
        return;
    }

    fillKernelsList(kernelsIDs);
}

void kernelDensityEstimator::setSamples(vector<std::shared_ptr<vector<double> >> *samples)
{
  this->samples = vector<std::shared_ptr<vector<double>>>(*samples);
}

unsigned long long kernelDensityEstimator::setClusters(std::vector<std::shared_ptr<cluster>> clusters)
{
  this->clusters = clusters;
  _spModifyingParameters = {{}};
  return this->clusters.size();
}

int kernelDensityEstimator::setSmoothingParameters(const std::vector<double> &smoothingParams)
{
  smoothingParameters.clear();

  for(double param : smoothingParams) smoothingParameters.push_back(param);

  return smoothingParameters.size();
}

void kernelDensityEstimator::setAdditionalMultipliers(std::vector<double> multipliers)
{
  this->additionalMultipliers = multipliers;
}

double kernelDensityEstimator::getValue(vector<double>* x)
{
    if(x == nullptr)
    {
        //qDebug() << "Argument is null pointer.";
        return -1.0;
    }

    if(x->size() == 0)
    {
        //qDebug() << "Argument is empty.";
        return -2.0;
    }

    return getProductKernelValue(x);
}

int kernelDensityEstimator::getDimension()
{
  return kernels.size();
}

void kernelDensityEstimator::updateSPModifyingParameters()
{
  int c = 1;
  double geometricMean = 0;
  double summaricWeight = 0;

  // Counting weighted geometric mean of lastKDEValue on clusters
  for(auto cl : clusters){
      summaricWeight += cl->getWeight();
      if(cl->_currentKDEValue > 0)
        geometricMean += cl->getWeight() * log(cl->_currentKDEValue);
  }

  if(summaricWeight > 0)
    geometricMean /= summaricWeight;

  geometricMean = exp(geometricMean);

  // Counting s_i, page 90 of THE BOOK
  for(auto cl : clusters){
      _spModifyingParameters[0].push_back(pow(cl->_currentKDEValue / geometricMean ,c));
  }
}

double kernelDensityEstimator::getProductKernelValue(vector<double> *x)
{
    weight = 0;

    // Check if values vector dimension is same size as kernels dimension
    if(x->size() != kernels.size())
    {
        return -4.0;
    }

    double result = getProductValuesFromClusters(x);

    for(double smoothingParameter : smoothingParameters)
        result /= smoothingParameter;

    result /= weight;

    return result;
}

double kernelDensityEstimator::getProductValuesFromClusters(vector<double>* x)
{
  double result = 0.0;

  // #pragma omp parallel for default(none) shared(x, kernels) reduction (+:result)
  for(size_t i = 0; i < clusters.size(); ++i)
  {
    auto c = clusters[i];

    double addend = getProductKernelAddendFromClusterIndex(i, x);

    if(_shouldConsiderWeights) {
      addend *= c->getCWeight();
      weight += c->getCWeight();
    } else {
      weight += 1;
    }

    if(clusters.size() == additionalMultipliers.size()) {
      addend *= additionalMultipliers[i];
    }
    result += addend;
  }

  return result;
}

/*
int kernelDensityEstimator::extractSampleFromCluster(std::shared_ptr<cluster> c, vector<double> *smpl)
{
  // This method assumes, that clustered sample has numerical values only
  smpl->clear();

  if(c.get()->getObject().get() == nullptr) return -3;

  std::shared_ptr<sample> obj = c->getObject();

  //qDebug() << (obj->attributesValues)["Val0"];

  std::unordered_map<std::string, std::string> attrVals
      = obj->attributesValues;

  for(auto attribute: *(obj->attirbutesOrder)){
    smpl->push_back(std::stod(attrVals[attribute]));
  }

  return smpl->size();
}
 */

vector<double> kernelDensityEstimator::extractSampleFromCluster(const std::shared_ptr<cluster> &c) const{
// This method assumes, that clustered sample has numerical values only
  vector<double> smpl = {};

  if(c.get()->getObject().get() == nullptr) return smpl;

  std::shared_ptr<sample> obj = c->getObject();
  std::unordered_map<std::string, std::string> attrVals
      = obj->attributesValues;

  for(auto attribute: *(obj->attirbutesOrder)) {
    smpl.push_back(std::stod(attrVals[attribute]));
  }

  return smpl;
}

double kernelDensityEstimator::getProductKernelAddendFromSample(vector<double> *sample, vector<double> *x)
{
  double result = 1.0;

  double restriction, component;
  bool hasRestriction;

  std::unique_ptr<vector<double>> tempValueHolder(new vector<double>());

  for(size_t i = 0; i < kernels.size(); ++i)
  {
      tempValueHolder->clear();
      tempValueHolder->push_back((x->at(i) - sample->at(i))/smoothingParameters.at(i));

      component = kernels.at(i)->getValue(tempValueHolder.get());

      try {
        hasRestriction = true;
        restriction = std::stod(carriersRestrictions.at(i));
      } catch (std::exception& e) {
        hasRestriction = false;
      }

      if(hasRestriction)
      {
          tempValueHolder->clear();
          tempValueHolder->push_back((x->at(i)+sample->at(i)-2*restriction)/smoothingParameters.at(i));
          component += kernels.at(i)->getValue(tempValueHolder.get());
          component *= partitionCharacteristicFunction(x->at(i), restriction);
      }

      result *= component;
  }

  return result;
}

double kernelDensityEstimator::getProductKernelAddendFromClusterIndex(int index, vector<double> *x)
{
  double result = 1.0;

  double restriction, component;
  bool hasRestriction;

  std::unique_ptr<vector<double>> tempValueHolder(new vector<double>());

  vector<double> s = extractSampleFromCluster(clusters[index]);;
  vector<double> *sample = &s;

  for(size_t i = 0; i < kernels.size(); ++i)
  {
    double argument = (x->at(i) - sample->at(i)) / smoothingParameters.at(i);
    if(_spModifyingParameters[i].size() == clusters.size())
        argument /= _spModifyingParameters[i][index];

    tempValueHolder->clear();
    tempValueHolder->push_back(argument);

    component = kernels.at(i)->getValue(tempValueHolder.get());
    if(_spModifyingParameters[i].size() == clusters.size())
        component /= _spModifyingParameters[i][index];

    try {
      hasRestriction = true;
      restriction = std::stod(carriersRestrictions.at(i));
    } catch (std::exception& e) {
      hasRestriction = false;
    }

    if(hasRestriction)
    {
        tempValueHolder->clear();
        tempValueHolder->push_back((x->at(i)+sample->at(i)-2*restriction)/smoothingParameters.at(i));
        component += kernels.at(i)->getValue(tempValueHolder.get());
        component *= partitionCharacteristicFunction(x->at(i), restriction);
    }

    result *= component;
  }

  return result;
}

double kernelDensityEstimator::getProductValuesFromSamples(vector<double> *x)
{
  double result = 0.0;

  for(std::shared_ptr<vector<double>> sample : samples)
    result += getProductKernelAddendFromSample(sample.get(), x);

  return result;
}

void kernelDensityEstimator::fillKernelsList(vector<int> *kernelsIDs)
{
  addProductKernelsToTheList(kernelsIDs);
}

void kernelDensityEstimator::addProductKernelsToTheList(vector<int> *kernelsIDs)
{
    for(int kernelID : *kernelsIDs)
    {
        switch (kernelID)
        {
            case NORMAL:
                kernels.push_back(kernelPtr(new normalKernel()));
            break;
            case TRIANGLE:
                kernels.push_back(kernelPtr(new triangleKernel()));
            break;
            case EPANECZNIKOW:
                kernels.push_back(kernelPtr(new epanecznikowKernel()));
            break;
            case DULL:
            default:
                kernels.push_back(kernelPtr(new dullKernel()));
            break;
        }
    }
}

int kernelDensityEstimator::partitionCharacteristicFunction(double carrier, double restriction)
{
    return carrier >= restriction ? 1 : 0;
}

void kernelDensityEstimator::compute_covariance_matrix(){

  std::vector<std::vector<double>> extracted_samples = {};

  for(auto c : clusters){
    extracted_samples.push_back(extractSampleFromCluster(c));
  }

  int dimension = extracted_samples[0].size();

  _covarianceMatrix = mat(dimension, dimension);

  for(int i = 0; i < dimension; ++i){
    for(int j = 0; j <= i; ++j){
      _covarianceMatrix(i,j)=_covarianceMatrix(j, i)=compute_weighted_covariance(i, j, extracted_samples);
    }
  }
}

double kernelDensityEstimator::compute_weighted_covariance(const int &i, const int &j, const vector<vector<double>> &data) const{
  double covariance = 0;
  double weights_sum = 0;

  double i_mean = compute_weighted_mean(i, data);
  double j_mean = compute_weighted_mean(j, data);

  for(int k = 0; k < data.size(); ++k){
    double weight = clusters[i]->getCWeight();
    vector<double> datum = data[k];

    covariance += weight * (datum[i] - i_mean) * (datum[j] - j_mean);

    weights_sum += weight;
  }

  return covariance / weights_sum;
}

double kernelDensityEstimator::compute_weighted_mean(const int &i, const vector<vector<double>> &data) const{
  double mean = 0;
  double weights_sum = 0;

  for(int j = 0; j < clusters.size(); ++j){
    double weight = clusters[j]->getCWeight();
    mean += data[j][i] * weight;
    weights_sum += weight;
  }

  return mean / weights_sum;
}

void kernelDensityEstimator::updateCovarianceMatrix(){
  if(clusters.size() < 2){
    return;
  }

  compute_covariance_matrix();
}

double kernelDensityEstimator::getRadialKernelValue(vector<double>* x) const{
  double value = 0;
  double weights_sum = 0;


  vec x_vec = vec(x->size());
  for(int i = 0; i < x->size(); ++i){
    x_vec(i) = (*x)[i];
  }


  // I'll only implement radial 2D kernel. This project has to be rewritten anyway.
  mat cov_inv = _covarianceMatrix.i();

  for(auto c : clusters){

    double weight = c->getCWeight();
    weights_sum += weight;

    vector<double> s = extractSampleFromCluster(c);
    vec c_vec = vec(x->size());

    for(int i = 0; i < x->size(); ++i){
      c_vec(i) = s[i];
    }

    double kernel_value = 0;

    vec v = (x_vec - c_vec);

    kernel_value = exp(as_scalar(v.t() * cov_inv * v) * smoothingParameters[0]);

    value += weight * kernel_value;
  }


  value /= weights_sum;

  value /= pow(smoothingParameters[0], x->size());

  value /= pow(2 * M_PI, - x->size() / 2);

  value /= sqrt(det(_covarianceMatrix));

  return value;

}
