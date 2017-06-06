#include "distributiondatareader.h"


distributionDataReader::distributionDataReader(distribution *sourceDistribution) : sourceDistribution(sourceDistribution){}

void distributionDataReader::getNextRawDatum(void *target)
{
    QVector<qreal>* targetPtr = static_cast<QVector<qreal>*>(target);
    targetPtr->clear();

    sourceDistribution->getValue(targetPtr);
}
