#include "plotLabelIntDataPreparator.h"

plotLabelIntDataPreparator::plotLabelIntDataPreparator(){ }

/** plotLabelIntDataPreparator::prepareValue
 * @brief Prepares int data. Ints are expected to be simply converted into
 * QString.
 * @param value -- int pointer holding value to be converted.
 * @return Converted int.
 */
QString plotLabelIntDataPreparator::prepareValue(void *value)
{
  // Note that according to DbC I'll not check if the input is correct.
  // I demand it is.
  return QString::number(*(static_cast<int*>(value)));
}


