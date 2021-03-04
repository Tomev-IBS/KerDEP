#ifndef PLOTLABELDOUBLEDATAPREPARATOR_H
#define PLOTLABELDOUBLEDATAPREPARATOR_H

#include "i_plotLabelDataPreparator.h"

class plotLabelDoubleDataPreparator : public i_plotLabelDataPreparator
{
public:
  plotLabelDoubleDataPreparator();
  QString prepareValue(void *value) override;
protected:
  QString formatNumberForDisplay(double number);
};

#endif // PLOTLABELDOUBLEDATAPREPARATOR_H
