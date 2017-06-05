#include "distributiondataparser.h"

distributionDataParser::distributionDataParser(){}

void distributionDataParser::parseData(void *source, void *target)
{
    QVector<qreal> *data    = static_cast<QVector<qreal>*>(source);
    QVector<qreal> *sample  = static_cast<QVector<qreal>*>(target);

    sample->clear();

    foreach(qreal datum, *data) sample->append(datum);
}

int distributionDataParser::addDatumToContainer(void *container)
{
    QVector<QVector<qreal>*> *samples = static_cast<QVector<QVector<qreal>*>*>(container);
    samples->push_back(new QVector<qreal>());

    return samples->size();
}

void distributionDataParser::writeDatumOnPosition(void *datum, void *container, int position)
{
    parseData(datum, static_cast<QVector<QVector<qreal>*>*>(container)->at(position));
}
