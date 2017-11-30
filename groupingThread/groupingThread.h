#ifndef GROUPINGTHREAD_H
#define GROUPINGTHREAD_H

#include "../Reservoir_sampling/sample.h"
#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"
#include "../groupingThread/kMedoidsAlgorithm/attributeData.h"

#include "medoidStoringAlgorithm/medoidStoringAlgorithm.h"

#include <QThread>

#include <unordered_map>
#include <vector>
#include <memory>

class groupingThread  : public QThread
{
  public:

    groupingThread(std::shared_ptr<std::vector<std::vector<std::shared_ptr<cluster> >> > medoidsStorage);
    void run();
    int initialize();

    int getObjectsForGrouping(std::vector<std::shared_ptr<sample> > samples);
    int setAttributesData(std::unordered_map<std::string, attributeData*>* attributesData);

  protected:

    std::vector<std::shared_ptr<sample>> objects;
    std::unordered_map<std::string, attributeData*>* attributesData;

    std::unique_ptr<medoidStoringAlgorithm> storingAlgorithm;

    std::shared_ptr<std::vector<std::vector<std::shared_ptr<cluster>>>> medoidsStorage;


};

#endif // GROUPINGTHREAD_H
