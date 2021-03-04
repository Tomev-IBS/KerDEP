#ifndef KMEANSALGORITHM_H
#define KMEANSALGORITHM_H


#include "./kMedoidsAlgorithm/groupingAlgorithm/distanceBasedGroupingAlgorithm.h"
#include "./kMedoidsAlgorithm/objectsDistanceMeasure.h"
#include "./kMedoidsAlgorithm/clustersDistanceMeasure.h"

#include <unordered_map>
#include "./kMedoidsAlgorithm/dataParser.h"


class kMeansAlgorithm : public distanceBasedGroupingAlgorithm
{
  public:

    enum initialMeansFindingStrategy
    {
      RANDOM = 0,
      RANDOM_ACCORDING_TO_DISTANCE = 1,
      FIRST_M_CLUSTERS = 2,
      LATEST_M_CLUSTERS = 3
    };

    kMeansAlgorithm(    int numberOfClusters,
                        std::shared_ptr<clustersDistanceMeasure> clusDistanceMeasure,
                        int initialMeansFindingStrategy,
                        std::shared_ptr<dataParser> parser);

    void groupObjects(std::vector<std::shared_ptr<sample>>* objects,
                      std::vector<std::shared_ptr<cluster>>* target);

    void groupClusters(std::vector<std::shared_ptr<cluster>>* clusters,
                       std::vector<std::shared_ptr<cluster>>* target);

    void setMedoids(std::vector<std::shared_ptr<cluster>>* newMedoids);

    std::vector<std::shared_ptr<cluster> > getMedoids(std::vector<std::shared_ptr<sample> > *objects);

    void generateClusteringFromMedoids(std::vector<std::shared_ptr<sample> > *objects,
                                       std::vector<std::shared_ptr<cluster>>* target);

  protected:

    unsigned int  numberOfClusters = 1;
    int           initialMeansFindingStrategy = RANDOM;

    time_t countStart;

    std::shared_ptr<dataParser> parser;

    std::vector<std::shared_ptr<cluster>> clusters;
    std::vector<std::shared_ptr<cluster>> means;

    std::unordered_map<std::string, double> similarityData;

    bool canGroupingBePerformed(unsigned int samplesSize);
    void clusterObjects(std::vector<std::shared_ptr<sample>> *objects);
    int performGrouping(std::vector<std::shared_ptr<cluster>>* target);
      int findInitialMeans();
        int getMeansFromFirstMClusters();
        int getMeansFromNewestClusters();
        int findRandomMeans();
        int findMeansAccordingToDistance();
          int addNewMeanToMeansVector(int meanIndex);
          int addNewMedoidAccordingToDistanceToMeansVector(
            std::vector<int>* nonMeanIndexes);
            int fillMeansWeights(std::vector<double> *weights,
                                 std::vector<int> *nonMeansIndexes);
            int fillMeansProbabilities(std::vector<double> *weights,
                                       std::vector<double> *probabilities);
      int applyNewMeans(std::vector<std::shared_ptr<cluster> > *target);
      int assignClustersToMeans(std::vector<std::shared_ptr<cluster> > *target);
      double countAssigmentError(std::vector<std::shared_ptr<cluster> > *target);
      int findNewMeans(std::vector<std::shared_ptr<cluster> > *target);
        int getAttributesKeysFromObjects(std::vector<std::string> *keys,
            std::vector<std::shared_ptr<sample> > *currentClusterObjects);
        double getAttributesMeanFromObjects(std::string attributesName,
               std::vector<std::shared_ptr<sample> > *currentClusterObjects);
        double getAttributesMeanFromSubclusters(std::string attributesName,
               std::vector<std::shared_ptr<cluster> > *currentClusterSubclusters);
};

#endif // KMEANSALGORITHM_H
