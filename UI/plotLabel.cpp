#include "plotLabel.h"

plotLabel::plotLabel(QCustomPlot *plot, const double &hOffset,
                     const double &vOffset, QString text, void* value)
  : _label(plot), _value(value)
{
  _label_font.setStyleHint(QFont::TypeWriter);
  _label.setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  _label.position->setType(QCPItemPosition::ptAxisRectRatio);
  _label.position->setCoords(hOffset, vOffset);
  _label.setFont(_label_font);
  _label.setText(text);
}

plotLabel::plotLabel(const plotLabel &pl) :
  _label(pl._label.parentPlot())
{
  _label_font.setStyleHint(QFont::TypeWriter);
  _label.setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  _label.position->setType(QCPItemPosition::ptAxisRectRatio);
  _label.position->setCoords(pl._label.position->coords().x(),
                             pl._label.position->coords().y());
  _label.setFont(_label_font);
  _label.setText(pl._label.text());
}

void plotLabel::setText(QString text)
{
  _label.setText(text);
}

void plotLabel::setFont(const QFont &newFont)
{
  _label.setFont(newFont);
}

void plotLabel::updateText()
{
  if(_value == nullptr) return;
}
