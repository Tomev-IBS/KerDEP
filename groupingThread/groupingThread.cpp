#include "groupingThread.h"

#include <QDebug>
#include <memory>

#include "kMedoidsAlgorithm/kMedoidsAlgorithm.h"
#include "kMedoidsAlgorithm/customObjectsDistanceMeasure.h"

#include "kMedoidsAlgorithm/attributesDistanceMeasures/numerical/gowersNumericalAttributesDistanceMeasure.h"
#include "kMedoidsAlgorithm/attributesDistanceMeasures/categorical/smdCategoricalAttributesDistanceMeasure.h"
#include "kMedoidsAlgorithm/clusterDistanceMeasures/completeLinkClusterDistanceMeasure.h"

#include "medoidStoringAlgorithm/medoidStoringAlgorithm.h"

groupingThread::groupingThread(std::vector<std::vector<std::shared_ptr<cluster>> >  *medoidsStorage)
{
  this->medoidsStorage = medoidsStorage;
}

void groupingThread::run()
{
  qDebug() << "Wątek totalnie działa.";

  int NUMBER_OF_MEDOIDS = 10;
  int MEDOIDS_FINDING_STRATEGY = RANDOM_ACCORDING_TO_DISTANCE;

  attributesDistanceMeasure* CADM = new smdCategoricalAttributesDistanceMeasure();
  attributesDistanceMeasure* NADM = new gowersNumericalAttributesDistanceMeasure(attributesData);
  objectsDistanceMeasure* ODM = new customObjectsDistanceMeasure(CADM, NADM, attributesData);
  clustersDistanceMeasure* CDM = new completeLinkClusterDistanceMeasure(ODM);

  std::shared_ptr<groupingAlgorithm> algorithm(new kMedoidsAlgorithm
                          (
                            NUMBER_OF_MEDOIDS,
                            CDM,
                            MEDOIDS_FINDING_STRATEGY
                          )
                        );

  medoidStoringAlgorithm* storingAlgorithm = new medoidStoringAlgorithm(algorithm);

  // TODO TR: It generates some errors. I'll work on it later on.
  //std::unique_ptr<medoidStoringAlgorithm> storingAlgorithm(new medoidStoringAlgorithm(algorithm));

  storingAlgorithm->findAndStoreMedoids(&objects, medoidsStorage);

  for(unsigned int i = 0; i < medoidsStorage->size(); ++i)
  {
    qDebug() << "Level: " << i << ". Custers number: "
             << medoidsStorage->at(i).size();
  }

  qDebug() << "Grouping finished and stored.";
}

int groupingThread::getObjectsForGrouping(std::vector<std::shared_ptr<sample>> samples)
{
  objects.clear();

  for(auto object : samples)
  {
    objects.push_back(std::shared_ptr<sample>(object));
  }

  return objects.size();
}

int groupingThread::setAttributesData(std::unordered_map<std::string, attributeData *> *attributesData)
{
  this->attributesData = attributesData;

  return attributesData->size();
}
