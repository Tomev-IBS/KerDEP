#include "groupingThread.h"

#include <QDebug>
#include <memory>

#include "kMedoidsAlgorithm/kMedoidsAlgorithm.h"
#include "kMedoidsAlgorithm/customObjectsDistanceMeasure.h"

#include "kMedoidsAlgorithm/attributesDistanceMeasures/numerical/gowersNumericalAttributesDistanceMeasure.h"
#include "kMedoidsAlgorithm/attributesDistanceMeasures/categorical/smdCategoricalAttributesDistanceMeasure.h"
#include "kMedoidsAlgorithm/clusterDistanceMeasures/completeLinkClusterDistanceMeasure.h"

groupingThread::groupingThread()
{

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

  groupingAlgorithm* algorithm = new kMedoidsAlgorithm( NUMBER_OF_MEDOIDS,
                                                        CDM,
                                                        MEDOIDS_FINDING_STRATEGY);


  std::vector<cluster> target;
  std::vector<sample*> samples;

  for(auto object : objects)
  {
    samples.push_back(object.get());
  }

  algorithm->groupObjects(&samples, &target);

  qDebug() << "Znaleziono " << target.size() << " klastrów.";

  for(cluster c : target)
  {
    c.getRepresentative()->print();
    qDebug() << c.size();
  }
}

int groupingThread::getObjectsForGrouping(std::vector<sample *> samples)
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
