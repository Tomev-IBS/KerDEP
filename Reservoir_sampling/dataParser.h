//
// Created by Tomev on 29.05.2017.
//

#ifndef RESERVOIRALGORITHM_DATAPARSER_H
#define RESERVOIRALGORITHM_DATAPARSER_H

class dataParser
{
  public:

    virtual void parseData(void* source, void *target) = 0;
    virtual int addDatumToContainer(void *container) = 0;
    virtual void writeDatumOnPosition(void *datum, void *container, int position) = 0;
};


#endif //RESERVOIRALGORITHM_DATAPARSER_H
