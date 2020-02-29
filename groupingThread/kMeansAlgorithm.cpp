#include "kMeansAlgorithm.h"

#include <random>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <limits>

#include <QDebug>

kMeansAlgorithm::kMeansAlgorithm(int numberOfClusters,
    std::shared_ptr<clustersDistanceMeasure> clusDistanceMeasure,
    int initialMeansFindingStrategy,
    std::shared_ptr<dataParser> parser)
{
  this->numberOfClusters = numberOfClusters;
  this->initialMeansFindingStrategy = initialMeansFindingStrategy;
  this->clusDistanceMeasure = clusDistanceMeasure;
  this->parser = parser;
}

void kMeansAlgorithm::groupObjects(
    std::vector<std::shared_ptr<sample> > *objects,
    std::vector<std::shared_ptr<cluster> > *target)
{
  if(canGroupingBePerformed(objects->size()))
  {
    clusterObjects(objects);
    performGrouping(target);
  }
}

void kMeansAlgorithm::groupClusters(
    std::vector<std::shared_ptr<cluster> > *clusters,
    std::vector<std::shared_ptr<cluster> > *target)
{
  if(canGroupingBePerformed(clusters->size()))
  {
    this->clusters = *clusters;
    performGrouping(target);
  }
}

bool kMeansAlgorithm::canGroupingBePerformed(unsigned int samplesSize)
{
  if(samplesSize < numberOfClusters)
  {
    std::cout << "Number of medoids is greater than objects number.\n";
    return false;
  }

  if(numberOfClusters <= 0)
  {
    std::cout << "Number of medoids is lower or equal 0.\n";
    return false;
  }

  return true;
}

void kMeansAlgorithm::clusterObjects(std::vector<std::shared_ptr<sample>> *objects)
{
  clusters.clear();

  for(unsigned long clusterNumber = 0; clusterNumber < objects->size(); ++clusterNumber)
  {
    clusters.push_back(
      std::make_shared<cluster>(cluster(clusterNumber + 1, objects->at(clusterNumber))));
  }
}

int kMeansAlgorithm::performGrouping(
    std::vector<std::shared_ptr<cluster> > *target)
{
  findInitialMeans();

  double oldError = 0.0;
  double newError = std::numeric_limits<double>::max();
  double errorThreshold = 1.0e-5;

  do
  {
    //std::cout << "Appling new means...\n";
    //qDebug() << "Appling new means...\n";

    /*
    qDebug() << "============Clusters==============";
    for(auto c : clusters)
      qDebug() << QString::fromStdString(c->getRepresentative()->attributesValues["Val0"]);

    qDebug() << "";

    qDebug() << "========== MEANS ===========";
    for(auto c : means){
        qDebug() << QString::fromStdString(c->getRepresentative()->attributesValues["Val0"]);
    }

    qDebug() << "";
    */

    applyNewMeans(target);

    //std::cout << "Assigning clusters to means...\n";
    //qDebug() << "Assigning clusters to means...\n";

    assignClustersToMeans(target);

    //std::cout << "Counting errors...\n";
    //qDebug() << "Counting errors...\n";

    oldError = newError;
    newError = countAssigmentError(target);

    //std::cout << "Finding new means...\n";
    //qDebug() << "Finding new means...\n";

    findNewMeans(target);

    /*
    for(auto c : *target){

      qDebug() << "Value" << QString::fromStdString(c->getRepresentative()->attributesValues["Val0"]);
      qDebug() << "Weight: " << c->getWeight();
      qDebug() << "Size: " << c->size();

      if(c->size() > 1){

        qDebug() << "============== SUBCLUSTERS =============";

        std::vector<clusterPtr> subclusters;

        c->getSubclusters(&subclusters);

        for(auto sc: subclusters){
            qDebug() << "Value: " << QString::fromStdString(sc->getRepresentative()->attributesValues["Val0"]);
            qDebug() << "Weight: " << sc->getWeight();

        }

        qDebug() << "=========== END SUBCLUSTERS ==============";

      }


        //std::cout
        //  << "Value: "  << c->getMean()->attributesValues["Val0"]
        //  << " Size: "  << c->size()
        //  << "Weight: " << c->getWeight() << std::endl;
    }
    */

    /*
    std::cout << "While condition: "
              << (oldError - newError > errorThreshold) << std::endl;

    std::cout << "Old error: " << oldError << std::endl
              << "New error: " << newError << std::endl;
    */

  } while( oldError - newError > errorThreshold);

  //std::cout << "Grouping finished.\n";

  std::cout << "Grouping finished.\nClusters:\n";

  return target->size();
}

int kMeansAlgorithm::findInitialMeans()
{
  switch(this->initialMeansFindingStrategy)
  {
    case LATEST_M_CLUSTERS:
      return getMeansFromNewestClusters();
    case FIRST_M_CLUSTERS:
      return getMeansFromFirstMClusters();
    case RANDOM_ACCORDING_TO_DISTANCE:
      return findMeansAccordingToDistance();
    case RANDOM:
    default:
      return findRandomMeans();
  }
}

int kMeansAlgorithm::getMeansFromFirstMClusters()
{
    means.clear();

    // Create distinct, new clusters for means
    for(int i = 0; i < numberOfClusters; ++i)
    {
      means.push_back(std::shared_ptr<cluster>(new cluster(i, clusters.at(i)->getObject())));
      means.at(means.size() - 1)->setRepresentative(clusters.at(i)->getRepresentative());
    }

    return means.size();
}

int kMeansAlgorithm::getMeansFromNewestClusters()
{
  means.clear();

  int oldestClusterIndex = 0;
  double smallestTimestamp = clusters[0]->timestamp;

  // Look for oldest cluster, and fill means in one loop
  for(int i = 0; i < clusters.size(); ++i){
    means.push_back(std::shared_ptr<cluster>(new cluster(i, clusters.at(i)->getObject())));
    means.at(means.size() - 1)->setRepresentative(clusters.at(i)->getRepresentative());

    if(clusters[i]->timestamp < smallestTimestamp){
      oldestClusterIndex = i;
      smallestTimestamp = clusters[i]->timestamp;
    }
  }

  // Then remove oldest cluster from means
  means.erase(means.begin() + oldestClusterIndex, means.begin() + oldestClusterIndex + 1);

  return means.size();
}

int kMeansAlgorithm::findRandomMeans()
{
  std::vector<std::shared_ptr<cluster>> possibleMeans = clusters;
  means.clear();

  std::random_device rd;     // only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

  int meanPosition;

  while(means.size() < numberOfClusters)
  {
    std::uniform_int_distribution<int>
      uniform_distribution(0, possibleMeans.size());

    meanPosition = uniform_distribution(rng);

    means.push_back(possibleMeans[meanPosition]);
    possibleMeans.erase(possibleMeans.begin() + meanPosition);
  }

  // Create distinct, new clusters for means
  for(int i = 0; i < numberOfClusters; ++i)
  {
    means.push_back(std::shared_ptr<cluster>(new cluster(i, means.at(0)->getObject())));
    means.at(means.size() - 1)->setRepresentative(means.at(0)->getRepresentative());
    means.erase(means.begin());
  }

  return means.size();
}

int kMeansAlgorithm::findMeansAccordingToDistance()
{
  means.clear();

  // Create and fill vector of clusters indexes
  std::vector<int> nonMeanIndexes;

  for(unsigned int i = 0; i < clusters.size(); ++i)
    nonMeanIndexes.push_back(i);

  // Select first medoid at random (uniformly) and add it to vector
  int nonMeanPosition = rand() % nonMeanIndexes.size();
  int newMeanIndex = nonMeanIndexes.at(nonMeanPosition);

  //qDebug() << "Adding mean.";
  addNewMeanToMeansVector(newMeanIndex);

  //qDebug() << "Erasing index.";
  nonMeanIndexes.erase(nonMeanIndexes.begin() + nonMeanPosition);

  //qDebug() << "Filling means indices.";
  while(means.size() < numberOfClusters)
  {
    addNewMedoidAccordingToDistanceToMeansVector(&nonMeanIndexes);
  }

  //qDebug() << "Creating means clusters.";
  // Create distinct, new clusters for means
  for(int i = 0; i < means.size(); ++i)
  {
    means.push_back(std::shared_ptr<cluster>(new cluster(i, means.at(0)->getObject())));
    means.at(means.size() - 1)->setRepresentative(means.at(0)->getRepresentative());
    means.erase(means.begin());
  }

  return means.size();
}

int kMeansAlgorithm::addNewMeanToMeansVector(int meanIndex)
{
  means.push_back(clusters.at(meanIndex));

  return means.size();
}

int kMeansAlgorithm::addNewMedoidAccordingToDistanceToMeansVector(
    std::vector<int> *nonMeanIndexes)
{
  std::vector<double> weights, probabilities;

  fillMeansWeights(&weights, nonMeanIndexes);
  fillMeansProbabilities(&weights, &probabilities);

  double r = ((double) rand() / (RAND_MAX));
  int newMeanIndexPosition = 0;

  //qDebug() << "Getting means position.";
  while(r > probabilities.at(newMeanIndexPosition))
  {
    ++newMeanIndexPosition;
    if(newMeanIndexPosition == probabilities.size() - 1)
      break;
  }

  //qDebug() << "Adding to means vector.";
  addNewMeanToMeansVector(nonMeanIndexes->at(newMeanIndexPosition));
  //qDebug() << "Erasing from nonmean indexes.";
  nonMeanIndexes->erase(nonMeanIndexes->begin() + newMeanIndexPosition);

  return 0;
}

int kMeansAlgorithm::fillMeansWeights(std::vector<double> *weights,
                                      std::vector<int> *nonMeansIndexes)
{
  double weight;

  weights->clear();

  for(int nonMedoidIndex : *nonMeansIndexes)
  {
    weight = 0;

    for(std::shared_ptr<cluster> mean : means)
      weight += clusDistanceMeasure->countClustersDistance(
        mean.get(),
        clusters.at(nonMedoidIndex).get());

    weight /= means.size();
    weights->push_back(weight);
  }

  return weights->size();
}

int kMeansAlgorithm::fillMeansProbabilities(std::vector<double> *weights,
                                            std::vector<double> *probabilities)
{
  probabilities->clear();

  double weightsSum = 0;

  for(double weight : *weights) weightsSum += weight;

  probabilities->push_back(weights->at(0) / weightsSum);

  double probability;

  for(unsigned int i = 1; i < weights->size(); ++i)
  {
    probability = probabilities->at(i - 1);
    probability += weights->at(i) / weightsSum;
    probabilities->push_back(probability);
  }

  return probabilities->size();
}

int kMeansAlgorithm::applyNewMeans(
    std::vector<std::shared_ptr<cluster>> *target)
{
  target->clear();

  for(std::shared_ptr<cluster> mean : means)
    target->push_back(mean);

  return target->size();
}

int kMeansAlgorithm::assignClustersToMeans(
  std::vector<std::shared_ptr<cluster> > *target)
{
  double currentMeanDistance = 0.0;
  double minMeanDistance = 1.0;
  int closestMeanIndex = 0;

  for(std::shared_ptr<cluster> c : clusters)
  {
    minMeanDistance =
      clusDistanceMeasure->countClustersDistance(target->at(0).get(),
                                                 c.get());
    closestMeanIndex = 0;

    for(int meanIndex = 1; meanIndex < target->size(); ++meanIndex)
    {
      currentMeanDistance =
        clusDistanceMeasure->countClustersDistance(target->at(meanIndex).get(),
                                                   c.get());

      if(currentMeanDistance < minMeanDistance)
      {
        minMeanDistance = currentMeanDistance;
        closestMeanIndex = meanIndex;
      }
    }

    target->at(closestMeanIndex)->addSubcluster(c);
  }
}

double kMeansAlgorithm::countAssigmentError(
  std::vector<std::shared_ptr<cluster> > *target)
{
  double error = 0.0;

  std::vector<std::shared_ptr<cluster>> subclusters;

  for(std::shared_ptr<cluster> c : *target)
  {
    c->getSubclusters(&subclusters);

    for(std::shared_ptr<cluster> subcluster : subclusters)
      error +=
        clusDistanceMeasure->countClustersDistance(subcluster.get(), c.get());
  }

  return error;
}

int kMeansAlgorithm::findNewMeans(std::vector<std::shared_ptr<cluster> > *target)
{
  std::vector<std::string> keys;
  std::vector<std::shared_ptr<sample>> newMeans;
  std::shared_ptr<sample> currentMean;
  std::vector<std::shared_ptr<sample>> currentClusterObjects;
  std::vector<std::shared_ptr<cluster>> currentClusterSubclusters;
  double attributesMean = 0.0;

  means.clear();

  for(int i = 0; i < target->size(); ++i)
  {
    currentClusterObjects.clear();
    target->at(i)->getObjects(&currentClusterObjects);
    target->at(i)->getSubclusters(&currentClusterSubclusters);

    keys.clear();

    getAttributesKeysFromObjects(&keys, &currentClusterObjects);

    parser->addDatumToContainer(&newMeans);
    currentMean = newMeans[newMeans.size() - 1];
    currentMean->attributesData = currentClusterObjects[0]->attributesData;

    /* As kmeans only works for continous (numerical) values, it can be
     * assumed that only numerical values are considered and I can count
     * means directly, for each attribute. */

    for(std::string key : keys)
    {
      if(currentClusterSubclusters.size() == 0)
        attributesMean =  std::stod(currentClusterObjects[0]->attributesValues[key]);
      else
        attributesMean = getAttributesMeanFromSubclusters(key, &currentClusterSubclusters);

      currentMean->attributesValues[key] = std::to_string(attributesMean);
    }

  }

  for(std::shared_ptr<sample> mean : newMeans)
    means.push_back(std::shared_ptr<cluster>(new cluster(std::shared_ptr<sample>(mean), true)));

  return means.size();
}

int kMeansAlgorithm::getAttributesKeysFromObjects(std::vector<std::string> *keys,
    std::vector<std::shared_ptr<sample> > *currentClusterObjects)
{
  for(std::shared_ptr<sample> currentObject : *currentClusterObjects)
  {
    for(auto kv : currentObject->attributesValues)
    {
      if(std::find(keys->begin(), keys->end(), kv.first) == keys->end())
        keys->push_back(kv.first); // add key if it doesn't
    }
  }

  return keys->size();
}

double kMeansAlgorithm::getAttributesMeanFromObjects(std::string attributesName,
       std::vector<std::shared_ptr<sample> > *currentClusterObjects)
{
  double mean = 0.0;

  for(std::shared_ptr<sample> currentObject : *currentClusterObjects)
    mean += stod(currentObject->attributesValues[attributesName]);

  return mean / currentClusterObjects->size();
}

double kMeansAlgorithm::getAttributesMeanFromSubclusters(std::string attributesName,
       std::vector<std::shared_ptr<cluster>> *currentClusterSubclusters)
{
  double mean = 0.0;
  double weightsSum = 0.0;
  double currentClusterWeight = 0.0;

  for(std::shared_ptr<cluster> c : *currentClusterSubclusters)
  {
    if(!c->representsObject()) continue;

    currentClusterWeight = c->getWeight();
    weightsSum += currentClusterWeight;
    mean += currentClusterWeight * stod(c->getObject()->attributesValues[attributesName]);
  }

  if(weightsSum != 0)
    mean /= weightsSum;

  return mean;
}

