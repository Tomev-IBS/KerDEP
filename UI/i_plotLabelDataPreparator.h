#ifndef I_PLOTLABELDATAPREPARATOR_H
#define I_PLOTLABELDATAPREPARATOR_H

#include <QObject>

/**
 * @brief The i_plotLabelDataPreparator class will be used for updating plot
 * labels, that can have various data to display (especially ints and reals
 * formatted in certain way). It's necessary for
 */

class i_plotLabelDataPreparator
{
public:
  // From given value type, create a QString to be printed.
  virtual QString prepareValue(void *value) = 0;
};

#endif // I_PLOTLABELDATAPREPARATOR_H
