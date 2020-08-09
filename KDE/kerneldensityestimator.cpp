#include "kerneldensityestimator.h"

#include <cmath>

#include <QDebug>

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
        return -1.0;
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
        return -1.0;
    }

    double result = getProductValuesFromClusters(x);

    for(double smoothingParameter : smoothingParameters)
        result /= smoothingParameter;

    result /= weight;

    return result;
}

double kernelDensityEstimator::getProductValuesFromClusters(vector<double>* x)
{
  double result = 0.0, addend;
  vector<double> sample;
  int i = 0, index = 0;

  for(std::shared_ptr<cluster> c : clusters)
  {
    addend = getProductKernelAddendFromClusterIndex(index++, x);

    if(_shouldConsiderWeights)
      addend *= c->getCWeight();

    if(clusters.size() == additionalMultipliers.size())
    {
      addend *= additionalMultipliers[i++];
    }

    if(_shouldConsiderWeights)
    {
      weight += c->getCWeight();
    }
    else
    {
      weight += 1;
    }

    result += addend;
  }

  return result;
}

int kernelDensityEstimator::extractSampleFromCluster(std::shared_ptr<cluster> c, vector<double> *smpl)
{
  // This method assumes, that clustered sample has numerical values only
  smpl->clear();

  if(c.get()->getObject().get() == nullptr) return -1;

  std::shared_ptr<sample> obj = c->getObject();

  //qDebug() << (obj->attributesValues)["Val0"];

  std::unordered_map<std::string, std::string> attrVals
      = obj->attributesValues;

  for(auto attrVal : attrVals)
    smpl->push_back(std::stod(attrVal.second));

  return smpl->size();
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

  vector<double> s;
  vector<double> *sample = &s;

  extractSampleFromCluster(clusters[index], sample);

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
