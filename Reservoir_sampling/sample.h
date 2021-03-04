#ifndef MEDOIDSTORINGALGORITHM_SAMPLE_H
#define MEDOIDSTORINGALGORITHM_SAMPLE_H

#include <string>
#include <unordered_map>
#include <vector>

#include "../groupingThread/kMedoidsAlgorithm/attributeData.h"

class sample
{
  public:
    std::unordered_map<std::string, std::string> attributesValues;
    std::unordered_map<std::string, attributeData*> *attributesData;
    std::vector<std::string> *attributesOrder;
};


#endif //MEDOIDSTORINGALGORITHM_SAMPLE_H
