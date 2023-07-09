#ifndef PLOTLABELDOUBLEDATAPREPARATOR_H
#define PLOTLABELDOUBLEDATAPREPARATOR_H

#include "i_plotLabelDataPreparator.h"

class plotLabelDoubleDataPreparator : public i_plotLabelDataPreparator
{
public:
  plotLabelDoubleDataPreparator(const unsigned int &decimal_numbers= 3, const bool &possibly_negative = true);
  QString prepareValue(void *value) override;
protected:
  QString formatNumberForDisplay(double number);
  bool possibly_negative_;
  unsigned int decimal_numbers_;
};

#endif // PLOTLABELDOUBLEDATAPREPARATOR_H
