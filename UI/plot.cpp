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
#include "plot.h"


Plot::Plot( QWidget *parent ): QwtPlot( parent ), d_alpha(0)
{
  plotLayout()->setAlignCanvasToScales( true );
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
  spectrograms.back()->setRenderThreadCount( 0 ); // use system specific thread count
  spectrograms.back()->setCachePolicy( QwtPlotRasterItem::PaintCache );
  spectrograms.back()->setDefaultContourPen(pen);
  spectrograms.back()->setData(data);
  spectrograms.back()->attach(this);
}

void Plot::setAxesLimit(const double &limit)
{
  // Only one argument, because it has to be square.

  setAxisScale(Qt::XAxis, -3, 11);
  setAxisScale(Qt::ZAxis, -3, 25);
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
    spectrogram->setDisplayMode( QwtPlotSpectrogram::ImageMode, on );
    spectrogram->setDefaultContourPen(
          on ? QPen( Qt::black, 0 ) : QPen( Qt::NoPen ) );
  }
}

void Plot::setAlpha( int alpha )
{
  // setting an alpha value doesn't make sense in combination
  // with a color map interpolating the alpha value
  d_alpha = alpha;
  for(auto spectrogram: spectrograms){

    spectrogram->setAlpha(alpha);
  }
}



