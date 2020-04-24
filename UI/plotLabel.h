#ifndef PLOTLABEL_H
#define PLOTLABEL_H

#include "QCustomPlot/qcustomplot.h"

class plotLabel
{
  public:
    plotLabel(QCustomPlot* plot, const double &hOffset,
            const double &vOffset, QString text, void* value = nullptr);
    plotLabel(const plotLabel &pl);
    void setText(QString text);
    void setFont(const QFont &newFont);
    void updateText();
  private:
    QFont _label_font = QFont("Courier New", 18);
    QCPItemText _label;
    void* _value;
};

#endif // PLOTLABEL_H
