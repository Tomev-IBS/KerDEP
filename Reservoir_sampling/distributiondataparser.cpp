#include <QDebug>

#include "distributiondataparser.h"

distributionDataParser::distributionDataParser()
{
    buffer = new QVector<qreal>();
}

void distributionDataParser::parseData(void *target)
{
    QVector<qreal> *data    = static_cast<QVector<qreal>*>(buffer);
    distributionDataSample *sample  = static_cast<distributionDataSample*>(target);

    sample->values.clear();

    foreach(qreal datum, *data) sample->values.append(datum);
}

int distributionDataParser::addDatumToContainer(std::vector<sample*> *container)
{
    container->push_back(new distributionDataSample());

    return container->size();
}

void distributionDataParser::writeDatumOnPosition(std::vector<sample *> *container, int position)
{
    parseData(container->at(position));
}
