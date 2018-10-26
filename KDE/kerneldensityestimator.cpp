#include "kerneldensityestimator.h"

#include <cmath>

#include "QDebug"

kernelDensityEstimator::kernelDensityEstimator(QVector<std::shared_ptr<QVector<qreal>>>* samples, QVector<qreal>* smoothingParameters, QVector<QString> *carriersRestrictions, int kernelType, QVector<int>* kernelsIDs)
    : kernelType(kernelType){

    this->samples               = QVector<std::shared_ptr<QVector<qreal>>>(*samples);
    this->smoothingParameters   = QVector<qreal>(*smoothingParameters);
    this->carriersRestrictions  = QVector<QString>(*carriersRestrictions);

    if(kernelsIDs->size() != smoothingParameters->size())
    {
        qDebug() << "Smoothing parameters and kernels number aren't equal.";
        return;
    }

    fillKernelsList(kernelsIDs);
}

void kernelDensityEstimator::setSamples(QVector<std::shared_ptr<QVector<qreal> >> *samples)
{
  this->samples = QVector<std::shared_ptr<QVector<qreal>>>(*samples);
}

int kernelDensityEstimator::setClusters(std::vector<std::shared_ptr<cluster>> clusters)
{
  this->clusters = clusters;

  updateSPModifyingParameters();

  return this->clusters.size();
}

int kernelDensityEstimator::setSmoothingParameters(std::vector<double> smoothingParams)
{
  smoothingParameters.clear();

  for(double param : smoothingParams) smoothingParameters.append(param);

  return smoothingParameters.size();
}

void kernelDensityEstimator::setAdditionalMultipliers(std::vector<double> multipliers)
{
  this->additionalMultipliers = multipliers;
}

qreal kernelDensityEstimator::getValue(QVector<qreal>* x)
{
    if(x == NULL)
    {
        qDebug() << "Argument is null pointer.";
        return -1.0;
    }

    if(x->size() == 0)
    {
        qDebug() << "Argument is empty.";
        return -1.0;
    }

    switch(kernelType)
    {
        case PRODUCT:
            return getProductKernelValue(x);
        break;
        case RADIAL:
            return getRadialKernelValue(x);
        break;
        default:
            qDebug() << "Kernel type not precised.";
            return -1.0;
        break;
    }
}

void kernelDensityEstimator::updateSPModifyingParameters()
{
  _spModifyingParameters.clear();

  std::vector<double> s = getSParameter();

  std::unordered_map<std::string, double> currentVar;

  for(std::shared_ptr<cluster> c : clusters)
  {
    std::vector<double> modParam;
    currentVar = c->getVariation();
    int i = 0;

    if(currentVar.size() == 0)
    {
      for(int dim = 0; dim < c->dimension(); ++dim)
        modParam.push_back(1);
    }
    else
    {
      for(auto kv : currentVar)
        modParam.push_back(pow(kv.second / s[i++] , -_modificationIntensivity));
    }

    _spModifyingParameters.push_back(modParam);
  }
}

std::vector<double> kernelDensityEstimator::getSParameter()
{
  double summaricWeight = 0;

  std::unordered_map<std::string, double> s;

  // initialize variance
  for(std::shared_ptr<cluster> c : clusters)
    summaricWeight += c->getWeight();

  std::unordered_map<std::string, attributeData*> *attrsData
      = clusters[0]->getRepresentative()->attributesData;

  std::string attrType;
  bool isAttributeNumerical;

  // Initialize
  for(auto kv : *(attrsData))
  {
    attrType = kv.second->getType();
    isAttributeNumerical = attrType == "numerical";
    if(isAttributeNumerical) s[kv.first] = 0;
  }

  for(std::shared_ptr<cluster> c : clusters)
  {
    std::unordered_map<std::string, double> currentVar = c->getVariation();

    for(auto kv : currentVar)
    {
      if(currentVar[kv.first] < 1e-5) continue;
      s[kv.first] = s[kv.first] + log(currentVar[kv.first]);
    }
  }

  std::vector<double> result;

  for(auto kv : s)
    result.push_back(exp(kv.second / summaricWeight));

  return result;
}

qreal kernelDensityEstimator::getProductKernelValue(QVector<qreal> *x)
{  
    weight = 0;

    // Check if values vector dimension is same size as kernels dimension
    if(x->size() != kernels.size())
    {
        qDebug() << "Kernel size and argument dimension doesn't match.";
        qDebug() << "Argument dimension: " << x->size()
                 << " Kernels dimension: " << kernels.size();
        return -1.0;
    }

    double result = 0.0;

    result += getProductValuesFromClusters(x);

    for(double smoothingParameter : smoothingParameters)
        result /= smoothingParameter;

    result /= weight;

    return result;
}

double kernelDensityEstimator::getProductValuesFromClusters(QVector<qreal>* x)
{
  double result = 0.f, addend;
  QVector<qreal> sample;
  int i = 0, index = 0;

  for(std::shared_ptr<cluster> c : clusters)
  {
    //extractSampleFromCluster(c, &sample);
    //addend = getProductKernelAddendFromSample(&sample, x);
    addend = getProductKernelAddendFromClusterIndex(index++, x);
    addend *= c.get()->getWeight();

    if(clusters.size() == additionalMultipliers.size())
      addend *= additionalMultipliers[i++];

    weight += (c.get()->getWeight());

    result += addend;
  }

  return result;
}

int kernelDensityEstimator::extractSampleFromCluster(std::shared_ptr<cluster> c, QVector<qreal> *smpl)
{
  // This method assumes, that clustered sample has numerical values only
  smpl->clear();

  if(c.get()->getObject().get() == nullptr) return -1;

  std::shared_ptr<sample> obj = c->getObject();

  //qDebug() << (obj->attributesValues)["Val0"];

  std::unordered_map<std::string, std::string> attrVals
      = obj->attributesValues;



  for(auto attrVal : attrVals)
    smpl->append(std::stod(attrVal.second));

  return smpl->size();
}

double kernelDensityEstimator::getProductKernelAddendFromSample(QVector<qreal> *sample, QVector<qreal> *x)
{
  double result = 1.0;

  qreal restriction, component;
  bool hasRestriction;

  std::unique_ptr<QVector<qreal>> tempValueHolder(new QVector<qreal>());

  for(int i = 0; i < kernels.size(); ++i)
  {
      tempValueHolder->clear();
      tempValueHolder->append((x->at(i) - sample->at(i))/smoothingParameters.at(i));

      component = kernels.at(i)->getValue(tempValueHolder.get());

      restriction = carriersRestrictions.at(i).toDouble(&hasRestriction);

      if(hasRestriction)
      {
          tempValueHolder->clear();
          tempValueHolder->append((x->at(i)+sample->at(i)-2*restriction)/smoothingParameters.at(i));
          component += kernels.at(i)->getValue(tempValueHolder.get());

          component *= partitionCharacteristicFunction(x->at(i), carriersRestrictions.at(i).toDouble());
      }

      result *= component;
  }

  return result;
}

double kernelDensityEstimator::getProductKernelAddendFromClusterIndex(int index, QVector<qreal> *x)
{
  double result = 1.0;

  qreal restriction, component;
  bool hasRestriction;

  std::unique_ptr<QVector<qreal>> tempValueHolder(new QVector<qreal>());

  QVector<qreal> s;
  QVector<qreal> *sample = &s;

  extractSampleFromCluster(clusters[index], sample);

  for(int i = 0; i < kernels.size(); ++i)
  {
      tempValueHolder->clear();
      tempValueHolder->append((x->at(i) - sample->at(i))/
                              (smoothingParameters.at(i) * _spModifyingParameters[index][i]));

      component = kernels.at(i)->getValue(tempValueHolder.get());
      component /= _spModifyingParameters[index][i];

      restriction = carriersRestrictions.at(i).toDouble(&hasRestriction);

      if(hasRestriction)
      {
          tempValueHolder->clear();
          tempValueHolder->append((x->at(i)+sample->at(i)-2*restriction)/smoothingParameters.at(i));
          component += kernels.at(i)->getValue(tempValueHolder.get());

          component *= partitionCharacteristicFunction(x->at(i), carriersRestrictions.at(i).toDouble());
      }

      result *= component;
  }

  return result;
}

double kernelDensityEstimator::getProductValuesFromSamples(QVector<qreal> *x)
{
  double result = 0.f;

  for(std::shared_ptr<QVector<qreal>> sample : samples)
    result += getProductKernelAddendFromSample(sample.get(), x);

  return result;
}

qreal kernelDensityEstimator::getRadialKernelValue(QVector<qreal> *x)
{
    // TODO TR: Code it when radial kernel is implemented
    x->clear();

    return -1.0;
}

void kernelDensityEstimator::fillKernelsList(QVector<int> *kernelsIDs)
{
    switch(kernelType)
    {
        case PRODUCT:
            addProductKernelsToTheList(kernelsIDs);
        break;
        case RADIAL:
            addRadialKernelsToTheList(kernelsIDs);
        break;
        default:
        break;
    }
}

void kernelDensityEstimator::addProductKernelsToTheList(QVector<int> *kernelsIDs)
{
    foreach(int kernelID, *kernelsIDs)
    {
        switch (kernelID)
        {
            case NORMAL:
                kernels.append(kernelPtr(new normalKernel()));
            break;
            case TRIANGLE:
                kernels.append(kernelPtr(new triangleKernel()));
            break;
            case EPANECZNIKOW:
                kernels.append(kernelPtr(new epanecznikowKernel()));
            break;
            case DULL:
            default:
                kernels.append(kernelPtr(new dullKernel()));
            break;
        }
    }
}

void kernelDensityEstimator::addRadialKernelsToTheList(QVector<int> *kernelsIDs)
{
  // TR TODO: Fill when radial kernel will be used
  kernelsIDs->clear();
}

int kernelDensityEstimator::partitionCharacteristicFunction(qreal carrier, qreal restriction)
{
    return carrier >= restriction ? 1 : 0;
}
