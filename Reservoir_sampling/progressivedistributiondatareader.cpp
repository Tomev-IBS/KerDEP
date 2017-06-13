#include "progressivedistributiondatareader.h"

progressiveDistributionDataReader::progressiveDistributionDataReader(distribution *source, qreal progressionSize) :
    sourceDistribution(source), progressionSize(progressionSize)
{}

void progressiveDistributionDataReader::getNextRawDatum(void *target)
{
    QVector<qreal>* targetPtr = static_cast<QVector<qreal>*>(target);
    targetPtr->clear();

    sourceDistribution->getValue(targetPtr);

    sourceDistribution->increaseMeans(progressionSize);
}

void progressiveDistributionDataReader::gatherAttributesData(void *attributes)
{
    // There are no attributes in distribution, just numbers;
    return;
}

bool progressiveDistributionDataReader::hasMoreData()
{
    // One can always generate more data from distribution.
    return true;
}
