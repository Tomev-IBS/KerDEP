#include "distributiondatareader.h"

distributionDataReader::distributionDataReader(distribution *sourceDistribution) : sourceDistribution(sourceDistribution){}

void distributionDataReader::getNextRawDatum(void *target)
{
    QVector<qreal>* targetPtr = static_cast<QVector<qreal>*>(target);
    targetPtr->clear();

    sourceDistribution->getValue(targetPtr);
}

void distributionDataReader::gatherAttributesData(void *attributes)
{
    // There are no attributes in distribution, just numbers;
    return;
}

bool distributionDataReader::hasMoreData()
{
    // One can always generate more data from distribution.
    return true;
}
