#include "plotLabelDoubleDataPreparator.h"

plotLabelDoubleDataPreparator::plotLabelDoubleDataPreparator(){ }

/** plotLabelDoubleDataPreparator::prepareValue
 * @brief Prepares double data. Real number are to be formatted with 6 digit
 * precission, thus formatNumberForDisplay method is used.
 * @param value -- int pointer holding value to be converted.
 * @return Converted int.
 */
QString plotLabelDoubleDataPreparator::prepareValue(void *value)
{
  // Note that according to DbC I'll not check if the input is correct.
  // I demand it is.
  auto number = *(static_cast<double*>(value));
  return formatNumberForDisplay(number);
}

/** plotLabelDoubleDataPreparator::formatNumberForDisplay
 * @brief Implementation of desired method for real values formatting.
 * According to PK the number should be displayed as #.######
 * @param number -- Double to be formatted.
 * @return Number formated as #.######.
 */
QString plotLabelDoubleDataPreparator::formatNumberForDisplay(double number)
{
  QString result = " ";

  if(number < 0) result = "";

  QStringList splitNumber = QString::number(number, 'f', 7).split(".");
  result += splitNumber[0];

  if(splitNumber.size() == 1) return result;

  result += ".";

  for(int i = 0; i < 6 && i < splitNumber[1].size(); ++i)
    result += splitNumber[1][i];

  return result;
}
