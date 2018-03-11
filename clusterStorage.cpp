#include "clusterStorage.h"

clusterStorage::clusterStorage()
{

}

unsigned int clusterStorage::removeUnpromisingClusters(double weightThreshold)
{
  unsigned int removedClustersNumber = 0;

  for(unsigned int layerNum = 0; layerNum < storage.size(); ++layerNum)
  {
    for(int index = storage[layerNum].size() - 1; index >= 0; --index)
    {
      if(shouldClusterBeRemoved(storage[layerNum][index], weightThreshold))
      {
        ++removedClustersNumber;
        storage[layerNum].erase(storage[layerNum].begin() + index);
      }
    }
  }

  return removedClustersNumber;
}

unsigned int clusterStorage::updateWeights(double coefficient)
{
  unsigned int updatedClustersNumber = 0;

  for(layer l : storage)
  {
    for(std::shared_ptr<cluster> c : l)
    {
      ++updatedClustersNumber;
      c->setWeight(c->getWeight() * coefficient);
    }
  }

  return updatedClustersNumber;
}

std::vector<std::shared_ptr<cluster>> clusterStorage::getAllClusters()
{
  std::vector<std::shared_ptr<cluster> > allClusters;

  for(layer l : storage)
  {
    for(std::shared_ptr<cluster> c : l) allClusters.push_back(c);
  }

  return allClusters;
}

std::vector<std::shared_ptr<cluster> > clusterStorage::getWeightyClusters(
    double weightThreshold)
{
  std::vector<std::shared_ptr<cluster> > weightyClusters;

  for(layer l : storage)
  {
    for(std::shared_ptr<cluster> c : l)
    {
      if(c->getWeight() >= weightThreshold) weightyClusters.push_back(c);
    }
  }

  return weightyClusters;
}

unsigned int clusterStorage::size()
{
  return getAllClusters().size();
}

bool clusterStorage::shouldClusterBeRemoved(std::shared_ptr<cluster> c,
                                            double weightThreshold)
{
  if(c->getWeight() < weightThreshold &&
     c->positiveTemporalDerivativeTimesInARow < unpromisingClustersThreshold)
    return true;

  return false;
}


