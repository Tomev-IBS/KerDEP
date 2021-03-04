#ifndef GROUPINGTHREAD_H
#define GROUPINGTHREAD_H

#include "../Reservoir_sampling/sample.h"
#include "../groupingThread/kMedoidsAlgorithm/groupingAlgorithm/cluster.h"
#include "../groupingThread/kMedoidsAlgorithm/attributeData.h"

#include "medoidStoringAlgorithm/medoidStoringAlgorithm.h"
#include "kMeansAlgorithm.h"

#include <QThread>

#include <unordered_map>
#include <vector>
#include <memory>

class groupingThread : public QThread
{
  public:

    groupingThread(std::vector<std::shared_ptr<cluster> > *medoidsStorage,
                   std::shared_ptr<dataParser> parser);
    void run();
    int initialize(int medoidsNumber, int bufferSize);

    int getObjectsForGrouping(std::vector<std::shared_ptr<sample> > samples);
    int getClustersForGrouping(std::vector<std::shared_ptr<cluster> > clusters);
    int setAttributesData(std::unordered_map<std::string, attributeData*>* attributesData);

  protected:

    std::vector<std::shared_ptr<sample>> objects;
    std::vector<std::shared_ptr<cluster>> clusters;

    std::shared_ptr<dataParser> parser;

    std::unordered_map<std::string, attributeData*>* attributesData;

    std::unique_ptr<medoidStoringAlgorithm> storingAlgorithm;

    std::vector<std::shared_ptr<cluster>> *medoidsStorage;


};

#endif // GROUPINGTHREAD_H
