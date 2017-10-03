#ifndef RESERVOIRALGORITHM_DATAPARSER_H
#define RESERVOIRALGORITHM_DATAPARSER_H

#include <vector>
#include <memory>

#include "sample.h"

class dataParser
{
  public:

    virtual void parseData(void *target) = 0;
    virtual int addDatumToContainer(std::vector<std::shared_ptr<sample>> *container) = 0;
    virtual void writeDatumOnPosition(std::vector<std::shared_ptr<sample>> *container, int position) = 0;
    virtual void setAttributesOrder(std::vector<std::string> *attributesOrder) = 0;

    void *buffer;

  protected:

    std::vector<std::string> *attributesOrder;
};


#endif //RESERVOIRALGORITHM_DATAPARSER_H
