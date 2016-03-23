

#ifndef TICMARKS_H
#define TICMARKS_H

#include "Whiteboard.h"
#include <cmath>
#include <sstream>

namespace ns_whiteboard 
{

enum tic_style { linear, semilogy, semilogx, loglog };

class TickingStrategyBase
{
public:
TickingStrategyBase( line * x_axis,
		     line * y_axis,
		    const double x_scaling = 1.0,
		    const double x_offset = 0.0,
		    const double y_scaling = 1.0,
		    const double y_offset = 0.0,
		    const double size = 5.0)
  : m_xAxis(x_axis),
    m_yAxis(y_axis),
    m_size(size),
    m_xScale(x_scaling),
    m_xOffset(x_offset),
    m_yScale(y_scaling),
    m_yOffset(y_offset) {}

virtual void CreateTicking() = 0;

virtual ~TickingStrategyBase() {}

double GetSize() { return m_size; }

vector<line>  GetXLines() { return m_xAxisTics; }
vector<line>  GetYLines() { return m_yAxisTics; }
vector<text>  GetXText() { return m_xTxt; }
vector<text>  GetYText() { return m_yTxt; }

protected:
line * m_xAxis;
line * m_yAxis;
double m_size;
double m_xScale, m_xOffset, m_yScale, m_yOffset;

vector<line> m_xAxisTics, m_yAxisTics;
vector<text> m_xTxt, m_yTxt;


};

class LinearTicking : public TickingStrategyBase
{
public:
LinearTicking( line * x_axis,
	       line * y_axis,
	      const double x_scaling = 1.0,
	      const double x_offset = 0.0,
	      const double y_scaling = 1.0,
	      const double y_offset = 0.0,
	      const double size = 5.0)
  : TickingStrategyBase(x_axis,
			y_axis,
			x_scaling,
			x_offset,
			y_scaling,
			y_offset,
			size) {}

virtual ~LinearTicking() {}

void CreateTicking();

};


void LinearTicking::CreateTicking()
{

  // in pixels
  double x_tic_width = m_xAxis->GetWidth();  
  double y_tic_width = m_yAxis->GetWidth();  

  // in pixels
  xy_coords x_start = m_xAxis->StartCoords();
  xy_coords x_stop = m_xAxis->StopCoords();

  xy_coords y_start = m_yAxis->StartCoords();
  xy_coords y_stop = m_yAxis->StopCoords();
  
  // in pixels
  double x_len = x_stop.first-x_start.first;
  double y_len = y_stop.second-y_start.second;
  
  // in data units
  double max_pos_x = m_xScale*x_stop.first + m_xOffset;
  double max_pos_y = m_yScale*y_stop.second + m_yOffset;

  // in data units
  double x_fac(0), y_fac(0), tic_x_spacing(0), tic_y_spacing(0);
  x_fac = std::floor(std::log10(max_pos_x));
  y_fac = std::floor(std::log10(max_pos_y));
  
  if ( max_pos_x > 1 )
    tic_x_spacing = x_fac*std::pow(10,x_fac-1);
  else
    tic_x_spacing = -1*x_fac*std::pow(10,x_fac);

  if ( max_pos_y > 1 )
    tic_y_spacing = y_fac*std::pow(10,y_fac-1);
  else
    tic_y_spacing = -1*y_fac*std::pow(10,y_fac);


  if (tic_x_spacing <= 0)
    tic_x_spacing = 0.1;
  if (tic_y_spacing <= 0)
    tic_y_spacing = 0.1;

  // in data units 
  double x_start_x = m_xScale*x_start.first + m_xOffset;
  double x_start_y = m_yScale*x_start.second + m_yOffset;

  double x_stop_x = m_xScale*x_stop.first + m_xOffset;
  double x_stop_y = m_yScale*x_stop.second + m_yOffset;

  double y_start_x = m_xScale*y_start.first + m_xOffset;
  double y_start_y = m_yScale*y_start.second + m_yOffset;

  double y_stop_x = m_xScale*y_stop.first + m_xOffset;
  double y_stop_y = m_yScale*y_stop.second + m_yOffset;



  // get the tic marks for the x-axis
  double pos(x_start_x+tic_x_spacing);
  while ( pos <= m_xScale*x_len+x_start_x )
  {
    // in pixels
    double pos_pix = (pos-m_xOffset)/m_xScale;
    
    xy_coords start = make_pair(pos_pix-x_tic_width/2.0, y_start.second-m_size/2.0);
    xy_coords stop = make_pair(pos_pix-x_tic_width/2.0, y_start.second+m_size/2.0);

    line l( start,
	    stop,
	    x_tic_width );

    m_xAxisTics.push_back(l);

    ostringstream ost;
    ost << (pos-x_start_x);

    string tic_label(ost.str());
    double font_size = 10.0;
    string font("Times-Roman");
    color c(black);

    xy_coords txt_start = make_pair(start.first,start.second-2.0*m_size);

    text t(txt_start,
	   tic_label,
	   c,
	   font_size,
	   font);

    m_xTxt.push_back(t);

    pos += tic_x_spacing;

  }


  // get the tic marks for the y-axis
  pos = y_start_y+tic_y_spacing;
  while ( pos <= m_yScale*y_len+y_start_y )
  {
    double pos_pix = (pos-m_yOffset)/m_yScale;
    xy_coords start = make_pair( x_start.first-m_size/2.0, pos_pix-y_tic_width/2.0);
    xy_coords stop = make_pair( x_start.first+m_size/2.0, pos_pix-y_tic_width/2.0);

    line l( start,
	    stop,
	    y_tic_width );

    m_yAxisTics.push_back(l);


    ostringstream ost;
    ost << (pos-y_start_y);

    string tic_label(ost.str());
    double font_size = 10.0;
    string font("Times-Roman");
    color c(black);    

    xy_coords txt_start = make_pair(start.first-1.5*m_size,start.second);

    text t(txt_start,
	   tic_label,
	   c,
	   font_size,
	   font,
	   90);

    m_yTxt.push_back(t);


    pos += tic_y_spacing;




  }

}

}


#endif
