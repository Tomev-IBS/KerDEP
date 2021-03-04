#include "plotLabel.h"

plotLabel::plotLabel(QCustomPlot *plot, const double &hOffset,
                     const double &vOffset, QString text, void* value,
                     std::shared_ptr<i_plotLabelDataPreparator> dataPreparator)
  : _label(plot), _text(text), _value(value), _dataPreparator(dataPreparator)
{
  _label_font.setStyleHint(QFont::TypeWriter);
  _label.setPositionAlignment(Qt::AlignTop|Qt::AlignLeft);
  _label.position->setType(QCPItemPosition::ptAxisRectRatio);
  _label.position->setCoords(hOffset, vOffset);
  _label.setFont(_label_font);

  if(_value == nullptr){
    _label.setText(text);
  } else {
    _label.setText(_text + _dataPreparator->prepareValue(_value));
  }
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
  _label.setText(_text + _dataPreparator->prepareValue(_value));
}
