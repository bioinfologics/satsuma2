

#ifndef AXES_H
#define AXES_H

#include "Whiteboard.h"
#include "TicMarks.h"


namespace ns_whiteboard 
{



class axes : public graphic
{

public:

axes() 
  : m_xLine(NULL),
    m_yLine(NULL),
    m_xLabel(NULL),
    m_yLabel(NULL),
    m_ticking(NULL),
    self_owned_lines(true) {}

axes( line * x_line, 
      line * y_line )
  : m_xLine(x_line),
    m_yLine(y_line),
    m_xLabel(NULL),
    m_yLabel(NULL),
    m_ticking(NULL),
    self_owned_lines(false) {}


void Draw( display_type *d );

void CreateTics(const double x_scaling=1.0,
		const double x_offset=0.0,
		const double y_scaling=1.0,
		const double y_offset=0.0,
		const tic_style t = linear);


void SetXAxis( const xy_coords &startCoords, 
	       const xy_coords &stopCoords,
	       const double width = 1.0,
	       const color &c = black );

void SetYAxis( const xy_coords &startCoords, 

	       const xy_coords &stopCoords,
	       const double width = 1.0,
	       const color &c = black );


void RescaleAxes() {}

void SetXAxisLabel( const string &label,
		    const color &c = black,
		    const double fontSize = 12.0,
		    const string &font = "Times-Roman");

void SetYAxisLabel( const string &label,
		    const color &c = black,
		    const double fontSize = 12.0,
		    const string &font = "Times-Roman");
 

protected:

line *m_xLine, *m_yLine;

text *m_xLabel, *m_yLabel;

TickingStrategyBase * m_ticking;

bool self_owned_lines;

public:
virtual ~axes() 
{
  if ( self_owned_lines )
  {
    delete m_xLine;
    delete m_yLine;
  }
  delete m_xLabel;
  delete m_yLabel;
  delete m_ticking;
}


};


void axes::CreateTics(const double x_scaling,
		      const double x_offset,
		      const double y_scaling,
		      const double y_offset,
		      const tic_style t)
{

  if ( m_xLine == NULL || m_yLine == NULL )
  {
    cout << "Please set the axes first" << endl;
    exit(-1);
  }
  
  double tic_size = 8*m_xLine->GetWidth();

  switch ( t )
  {
    case linear:
      m_ticking = new LinearTicking(m_xLine,
				    m_yLine,
				    x_scaling,
				    x_offset,
				    y_scaling,
				    y_offset,
				    tic_size);
      break;
      
    // to be implemented
    case semilogx:
      break;
    case semilogy:
      break;
    case loglog:
      break;

  }

  m_ticking->CreateTicking();
}


void axes::SetXAxisLabel( const string &label,
			  const color &c,
			  const double fontSize,
			  const string &font )
{

  if ( m_xLine == NULL )
  {
    cout << "SetXAxis first, please." << endl;
    exit(-1);
  }

  if ( m_ticking == NULL )
  {
    cout << "due to stupidity, please CreateTics(...) first" << endl;
    exit(-1);
    //    CreateTics();
  }

  xy_coords axis_start = m_xLine->StartCoords();
  xy_coords axis_stop = m_xLine->StopCoords();


  double label_x_pos = axis_start.first + (double) (axis_stop.first-axis_start.first)/2 - (double) label.size()/2;
  double label_y_pos = axis_start.second - 2.0*m_ticking->GetSize() - 2.0*fontSize;

  xy_coords label_coords = make_pair(label_x_pos,label_y_pos);  

  m_xLabel = new text( label_coords,
		       label,
		       c,
		       fontSize,
		       font);
		       

}


void axes::SetYAxisLabel( const string &label,
			  const color &c,
			  const double fontSize,
			  const string &font )
{

  if ( m_yLine == NULL )
  {
    cout << "SetYAxis first, please." << endl;
    exit(-1);
  }

  if ( m_ticking == NULL )
  {
    cout << "due to stupidity, please CreateTics(...) first" << endl;
    exit(-1);
    //    CreateTics();
  }


  xy_coords axis_start = m_yLine->StartCoords();
  xy_coords axis_stop = m_yLine->StopCoords();

  double  label_x_pos = axis_start.first - 1.0*m_ticking->GetSize() - 2.0*fontSize;
  double label_y_pos = axis_start.second 
    + (double) (axis_stop.second - axis_start.second)/2 - (double) label.size()/2;

  xy_coords label_coords = make_pair(label_x_pos,label_y_pos);  

  m_yLabel = new text( label_coords,
		       label,
		       c,
		       fontSize,
		       font,
		       90);
		       

}

void axes::SetXAxis( const xy_coords &startCoords, 
		     const xy_coords &stopCoords,
		     const double width,
		     const color &c )
{
  m_xLine = new line(startCoords,
		     stopCoords,
		     width,
		     c);
  
}

void axes::SetYAxis( const xy_coords &startCoords, 
		     const xy_coords &stopCoords,
		     const double width,
		     const color &c)
{
  m_yLine = new line(startCoords,
		     stopCoords,
		     width,
		     c);
  
}



void axes::Draw( display_type *d )
{

  if ( m_xLine == NULL || m_yLine == NULL || m_xLabel == NULL || m_yLabel == NULL )
  {
    cout << "bad form" << endl;
    exit(-1);
  }

  vector<line>  xtics = m_ticking->GetXLines();
  vector<line>  ytics = m_ticking->GetYLines();
  
  vector<text>  xticlabs = m_ticking->GetXText();
  vector<text>  yticlabs = m_ticking->GetYText();

  if ( xtics.size() != xticlabs.size() )
  {
    cout << "xtic irregularity" << endl;
    exit(-1);
  }
  if ( ytics.size() != yticlabs.size() )
  {
    cout << "ytic irregularity" << endl;
    exit(-1);
  }


  m_xLine->Draw(d);
  m_yLine->Draw(d);
  m_xLabel->Draw(d);
  m_yLabel->Draw(d);


  for ( int i=0; i<(int) xtics.size(); ++i )
  {
    xtics[i].Draw(d);
    xticlabs[i].Draw(d);
  }

  for ( int i=0; i<(int) ytics.size(); ++i )
  {
    ytics[i].Draw(d);
    yticlabs[i].Draw(d);
  }


}






}



#endif
