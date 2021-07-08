#include <qprinter.h>
#include <qprintdialog.h>
#include <qnumeric.h>
#include <qwt_color_map.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_scale_widget.h>
#include <qwt_scale_draw.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_layout.h>
#include <qwt_plot_renderer.h>
#include <qwt_compat.h>
#include "plot.h"


class ColorMap: public QwtLinearColorMap{
  public:
    ColorMap():
        QwtLinearColorMap( Qt::white, Qt::black )
    {
      addColorStop(0.001, Qt::darkBlue);
      addColorStop( 0.005, Qt::blue);
      addColorStop( 0.01, Qt::cyan);
      addColorStop( 0.02, Qt::green);
      addColorStop( 0.035, Qt::yellow);
      addColorStop( 0.05, Qt::red);
      addColorStop(0.5, Qt::darkRed);
    }
};

Plot::Plot( QWidget *parent ): QwtPlot( parent ), d_alpha(0)
{
  plotLayout()->setAlignCanvasToScales( true );

  if(show_color_map_) {
    QwtDoubleInterval range(0, 0.2);
    QwtScaleWidget *scale = this->axisWidget(QwtPlot::yRight);
    scale->setColorBarEnabled(true);

    scale->setColorMap(range, new ColorMap());

    this->setAxisScale(QwtPlot::yRight, 0, 0.2);
    this->enableAxis(QwtPlot::yRight);
  }
}

void Plot::clearSpectrograms()
{
  for(int j = spectrograms.size() - 1; j >= 0; --j){
    spectrograms[j]->detach();
    delete spectrograms[j];
    spectrograms.pop_back();
  }
}

void Plot::setContours(const QList<double> &contourLevels)
{
  for(auto spectrogram: spectrograms){
    spectrogram->setContourLevels(contourLevels);
  }
}

void Plot::addQwtPlotSpectrogram(SpectrogramData *data, const QPen &pen)
{
  spectrograms.push_back(new QwtPlotSpectrogram());
  spectrograms.back()->setRenderThreadCount( 1 ); // use system specific thread count
  spectrograms.back()->setCachePolicy( QwtPlotRasterItem::PaintCache );
  spectrograms.back()->setDefaultContourPen(pen);
  spectrograms.back()->setData(data);
  if(show_color_map_) {
    spectrograms.back()->setColorMap(new ColorMap());
  }
  spectrograms.back()->attach(this);
}

void Plot::setAxesLimit(const double &limit)
{
  // Only one argument, because it has to be square.

  // First column - HORIZONTAL
  //setAxisScale(Qt::ZAxis, -40, 40); // 25, 31, 51, for classical, slower and lazy paths. Pm 6 for v = 0
  //setAxisTitle(Qt::ZAxis, "Temperature [C]");
  setAxisScale(Qt::ZAxis, -5, 55); // 25, 31, 51, for classical, slower and lazy paths. Pm 6 for v = 0
  setAxisTitle(Qt::ZAxis, "");

  // Second column - VERTICAL
  setAxisScale(Qt::XAxis, -5, 55); // 11, 14, 24, for classical, slower and lazy paths. Pm 3 for v = 0
  setAxisTitle(Qt::XAxis, "");
  //setAxisScale(Qt::XAxis, 0, 100);
  //setAxisTitle(Qt::XAxis, "Humidity [%]");
}

void Plot::showContour( bool on )
{
  for(auto spectrogram : spectrograms){
    spectrogram->setDisplayMode( QwtPlotSpectrogram::ContourMode, on );
  }
}

void Plot::showSpectrogram( bool on )
{
  for(auto spectrogram: spectrograms){
    spectrogram->setDefaultContourPen( on ? QPen( Qt::black, 0 ) : QPen( Qt::NoPen ) );
  }
}

void Plot::setAlpha( int alpha )
{
  // setting an alpha value doesn't make sense in combination with a color map interpolating the alpha value
  d_alpha = alpha;
  for(auto spectrogram: spectrograms){
    spectrogram->setAlpha(alpha);
  }
}

void Plot::replot() {
  if(spectrograms.size() > 0 && show_color_map_) {
    spectrograms.back()->setColorMap(new ColorMap());
  }
  QwtPlot::replot();
}

void Plot::ShowColorMap(const bool &show_color_map) {
  show_color_map_ = show_color_map;
}



