#include "distributiondataparser.h"

#include <QDebug>

distributionDataParser::distributionDataParser()
{
    buffor = new QVector<qreal>();
}

void distributionDataParser::parseData(void *target)
{
    QVector<qreal> *data    = static_cast<QVector<qreal>*>(buffor);
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

void distributionDataParser::writeDatumOnPosition(void *container, int position)
{
    parseData(static_cast<QVector<QVector<qreal>*>*>(container)->at(position));
}
