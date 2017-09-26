#ifndef RESERVOIRALGORITHM_DATAREADER_H
#define RESERVOIRALGORITHM_DATAREADER_H

#include <string>
#include <vector>
#include <unordered_map>

#include "../groupingThread/kMedoidsAlgorithm/attributeData.h"

class dataReader
{
  public:

    virtual void getNextRawDatum(void *target) = 0;
    virtual void gatherAttributesData(void *attributes) = 0;
    virtual bool hasMoreData() = 0;

    virtual std::vector<std::string>* getAttributesOrder() = 0;
};


#endif //RESERVOIRALGORITHM_DATAREADER_H
