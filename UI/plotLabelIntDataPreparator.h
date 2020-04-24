#ifndef PLOTLABELINTDATAPREPARATOR_H
#define PLOTLABELINTDATAPREPARATOR_H


#include "i_plotLabelDataPreparator.h"

class plotLabelIntDataPreparator : public i_plotLabelDataPreparator
{
public:
  plotLabelIntDataPreparator();
  QString prepareValue(void *value) override;
};

#endif // PLOTLABELINTDATAPREPARATOR_H
