#include "distributiondatareader.h"


distributionDataReader::distributionDataReader(distribution *sourceDistribution) : sourceDistribution(sourceDistribution){}

distributionDataReader::getNextRawDatum(void *target)
{
    sourceDistribution->getValue(static_cast<QVector<qreal>*>(target));
}
