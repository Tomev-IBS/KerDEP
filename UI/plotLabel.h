#ifndef PLOTLABEL_H
#define PLOTLABEL_H

#include "QCustomPlot/qcustomplot.h"

class plotLabel
{
  public:
    plotLabel(QCustomPlot* plot, const double &hOffset, const double &vOffset, QString text);
    plotLabel(const plotLabel &pl);
    void setText(QString text);
  private:
    QFont _label_font = QFont("Courier New", 18);
    QCPItemText _label;
};

#endif // PLOTLABEL_H
