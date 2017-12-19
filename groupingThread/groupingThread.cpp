#include "groupingThread.h"

#include <QDebug>
#include <memory>

#include "kMedoidsAlgorithm/kMedoidsAlgorithm.h"
#include "kMedoidsAlgorithm/customObjectsDistanceMeasure.h"

#include "kMedoidsAlgorithm/attributesDistanceMeasures/numerical/gowersNumericalAttributesDistanceMeasure.h"
#include "kMedoidsAlgorithm/attributesDistanceMeasures/categorical/smdCategoricalAttributesDistanceMeasure.h"
#include "kMedoidsAlgorithm/clusterDistanceMeasures/completeLinkClusterDistanceMeasure.h"



groupingThread::groupingThread(std::shared_ptr<std::vector<std::vector<std::shared_ptr<cluster>>>> medoidsStorage)
{
  this->medoidsStorage = medoidsStorage;
}

int groupingThread::initialize()
{
  int NUMBER_OF_MEDOIDS = 10;
  int MEDOIDS_FINDING_STRATEGY = RANDOM_ACCORDING_TO_DISTANCE;

  attributesDistanceMeasure* CADM = new smdCategoricalAttributesDistanceMeasure();
  attributesDistanceMeasure* NADM = new gowersNumericalAttributesDistanceMeasure(attributesData);
  objectsDistanceMeasure* ODM = new customObjectsDistanceMeasure(CADM, NADM, attributesData);
  std::shared_ptr<clustersDistanceMeasure> CDM(new completeLinkClusterDistanceMeasure(ODM));

  std::shared_ptr<groupingAlgorithm> algorithm(new kMedoidsAlgorithm
                          (
                            NUMBER_OF_MEDOIDS,
                            CDM,
                            MEDOIDS_FINDING_STRATEGY
                          )
                        );

  storingAlgorithm.reset(new medoidStoringAlgorithm(algorithm));

  return 0;
}

void groupingThread::run()
{
  qDebug() << "Wątek totalnie działa.";

  //storingAlgorithm->findAndStoreMedoidsFromObjects(&objects, medoidsStorage);
  storingAlgorithm->findAndStoreMedoidsFromClusters(&clusters, medoidsStorage);

  for(unsigned int i = 0; i < medoidsStorage->size(); ++i)
  {
    qDebug() << "Level: " << i << ". Custers number: "
             << medoidsStorage->at(i).size();
  }

  qDebug() << "Grouping finished and medoids stored.";
}

int groupingThread::getObjectsForGrouping(std::vector<std::shared_ptr<sample>> samples)
{
  objects.clear();

  for(auto object : samples)
    objects.push_back(std::shared_ptr<sample>(object));

  return objects.size();
}

int groupingThread::getClustersForGrouping(std::vector<std::shared_ptr<cluster> > clusters)
{
  this->clusters.clear();

  this->clusters.insert(this->clusters.end(), clusters.begin(), clusters.end());

  return this->clusters.size();
}

int groupingThread::setAttributesData(std::unordered_map<std::string, attributeData *> *attributesData)
{
  this->attributesData = attributesData;

  return attributesData->size();
}


