#ifndef MEDOIDSTORINGALGORITHM_SAMPLE_H
#define MEDOIDSTORINGALGORITHM_SAMPLE_H

#include <string>
#include <unordered_map>

class sample
{
  public:
    std::unordered_map<std::string, std::string> attributesValues;
    virtual void print() = 0;
};


#endif //MEDOIDSTORINGALGORITHM_SAMPLE_H
