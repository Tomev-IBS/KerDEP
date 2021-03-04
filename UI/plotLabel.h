#ifndef PLOTLABEL_H
#define PLOTLABEL_H

#include "QCustomPlot/qcustomplot.h"
#include "i_plotLabelDataPreparator.h"

class plotLabel
{
  public:
    plotLabel(QCustomPlot* plot, const double &hOffset,
            const double &vOffset, QString text, void* value = nullptr,
            std::shared_ptr<i_plotLabelDataPreparator> dataPreparator = std::shared_ptr<i_plotLabelDataPreparator>());
    plotLabel(const plotLabel &pl);
    void setText(QString text);
    void setFont(const QFont &newFont);
    void updateText();
  private:
    QString _text;
    QFont _label_font = QFont("Courier New", 18);
    QCPItemText _label;
    void* _value;
    std::shared_ptr<i_plotLabelDataPreparator> _dataPreparator;
};

#endif // PLOTLABEL_H
