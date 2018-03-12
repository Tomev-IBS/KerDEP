#ifndef CLUSTERSTORAGE_H
#define CLUSTERSTORAGE_H

#include "groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"

#include <vector>
#include <memory>

typedef std::vector<std::shared_ptr<cluster>> layer;

class clusterStorage
{
  public:
    clusterStorage();

    unsigned int removeUnpromisingClusters(double weightThreshold);
    unsigned int updateWeights(double coefficient);
    std::vector<std::shared_ptr<cluster>> getAllClusters();
    std::vector<std::shared_ptr<cluster>> getWeightyClusters(
        double weightThreshold);

    unsigned int size();
    void clear();
    void clearLevel(unsigned int level);

    int addCluster(std::shared_ptr<cluster> c, unsigned int level);

    layer operator [](int layerIndex);

  protected:

    unsigned int unpromisingClustersThreshold= 5;

    std::vector<layer> storage;

    bool shouldClusterBeRemoved(std::shared_ptr<cluster> c,
                                double weightThreshold);


};

#endif // CLUSTERSTORAGE_H
