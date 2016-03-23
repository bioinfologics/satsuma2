#ifndef WHITEBOARD_H
#define WHITEBOARD_H

#include "Color.h"
#include "Eps.h"

#include "../base/SVector.h"
#include <string>

namespace ns_whiteboard {

  typedef pair<double,double> xy_coords;

  enum vert_align { align_top, align_middle, align_bottom };
  enum horiz_align { align_left, align_center, align_right };

class display_type
{
 public:
  virtual void DrawPoint( const xy_coords &coord ) = 0;
  virtual void DrawPoint( double x, double y ) = 0;
  virtual void DrawLine( const xy_coords &startCoord,
			 const xy_coords &stopCoord ) = 0;
  virtual void DrawRect( const xy_coords &topLeftCoord,
			 const xy_coords &bottomRightCoord ) = 0;
  virtual void DrawText( const xy_coords& coords,
                         const string& chars,
                         double rot_angle,
                         double point_size,
                         const vert_align vAlign,
                         const horiz_align hAlign ) = 0;
  virtual void DrawArc( const xy_coords &coords,
			const double radius,
			const double startAngle,
			const double stopAngle ) = 0;

  virtual const color& GetColor() const = 0;
  virtual void SetColor( const color& C ) = 0;

  virtual double GetLineWidth() const = 0;
  virtual void SetLineWidth( const double width ) = 0;

  virtual const string& GetFontName() const = 0;
  virtual double         GetFontSize() const = 0;
  
  virtual void SetFont( const string& font, const double size ) = 0;

  virtual void Finish() = 0;
public:
  virtual ~display_type(){}
};


class graphic 
{
 public:
  graphic() {}

  virtual void Draw( display_type *d ) = 0;
public:
  virtual ~graphic(){}
};


class whiteboard
{
 public:
  whiteboard( ) {}
 
  void Add( graphic *g )
    {
      m_graphics.push_back(g);
    }

  void DisplayOn( display_type *d )
    {
      for ( int i=0; i<(int) m_graphics.size(); ++i )
	m_graphics[i]->Draw( d );
      d->Finish();
    }

  ///Clear all data, deleting all added pointers.
  void DeletePointers() {
    for (int i=0; i != m_graphics.isize(); ++i) {
      delete m_graphics[i];
    }
    m_graphics.clear();
  }

  /// Number of graphics pointers contained.
  int size() const { return m_graphics.size(); }

  ///Clear all data but do not delete pointers: creator is responsible 
  /// for them.
  void clear() { m_graphics.clear(); }

 private:

  svec<graphic *> m_graphics;
};



class point : public graphic
{
 public:
  point( const xy_coords& coords,
	 const double size = 1.0,
	 const color& c = black ) 
    : m_coords(coords),
      m_size(size),
      m_color(c) {}

  void Draw( display_type *d )
    {
      d->SetLineWidth(m_size);
      d->SetColor(m_color);
      d->DrawPoint(m_coords);
    }

 protected:
  xy_coords m_coords;
  double m_size;
  color m_color;
public:
  virtual ~point(){}
};

///Point class with double coordinates to avoid rounding errors.
class dpoint : public graphic
{
 public:
  dpoint( double x, double y,
	 const double size = 1.0,
	 const color& c = black ) 
    : x_(x), y_(y),
      m_size(size),
      m_color(c) {}

  void Draw( display_type *d )
    {
      d->SetLineWidth(m_size);
      d->SetColor(m_color);
      d->DrawPoint(x_,y_);
    }

 protected:
  double x_;
  double y_;
  double m_size;
  color m_color;
public:
  virtual ~dpoint(){}
};


class line : public graphic
{
 public:
  line( const xy_coords& startCoords, 
	const xy_coords& stopCoords,
	const double width = 1.0,
	const color& c = black ) 
    : m_startCoords(startCoords),
      m_stopCoords(stopCoords),
      m_width(width),
      m_color(c)
    {}

  void Draw( display_type *d )
    {
      d->SetLineWidth(m_width);
      d->SetColor(m_color);
      d->DrawLine(m_startCoords,m_stopCoords);
    }


xy_coords StartCoords() const { return m_startCoords; }
xy_coords StopCoords() const { return m_stopCoords; }
double GetWidth()  const { return m_width; }

 protected:
  xy_coords m_startCoords, m_stopCoords;
  double m_width;
  color m_color;
public:
  virtual ~line(){}
};



class rect : public graphic
{
 public:
  rect( const xy_coords& topLeftCoords, 
	const xy_coords& bottomRightCoords,
	const color& c = black ) 
    : m_topLeftCoords(topLeftCoords),
      m_bottomRightCoords(bottomRightCoords),
      m_color(c)
    {}

  void Draw( display_type *d )
    {
      d->SetColor(m_color);
      d->DrawRect(m_topLeftCoords,m_bottomRightCoords);
    }

 protected:
  xy_coords m_topLeftCoords, m_bottomRightCoords;
  color m_color;
public:
  virtual ~rect(){}
};



// coords : x-y coordinates of the arc's center of curvature
// radius : radius of curvature
// startAngle, stopAngle : measured counterclockwise from positive x-axis
class arc : public graphic
{
 public:
  arc( const xy_coords & coords,
       const double radius,
       const double startAngle,
       const double stopAngle,
       const double width = 1.0,
       const color &c = black )
    : m_coords(coords),
      m_radius(radius),
      m_startAngle(startAngle),
      m_stopAngle(stopAngle),
      m_width(width),
      m_color(c)
    {}

  void Draw( display_type *d )
    {
      d->SetLineWidth(m_width);
      d->SetColor(m_color);
      d->DrawArc(m_coords,m_radius,m_startAngle,m_stopAngle);
    }

 protected:
  xy_coords m_coords;
  double m_radius, m_startAngle, m_stopAngle, m_width;
  color m_color;
public:
  virtual ~arc(){}
};


class text : public graphic
{
 public:
  text( const xy_coords& coords,
	const string& chars,
	const color& c = black,
	const double fontSize = 12.0,
	const string& fontName = "Times-Roman",
	const double rotAngle = 0.0,
	bool alignLeft = false)
    : m_coords(coords),
      m_chars(chars),
      m_color(c),
      m_fontSize(fontSize),
      m_fontName(fontName),
      m_rotAngle(rotAngle),
      m_vertAlign(align_bottom),
      m_horizAlign(align_center)
  {
    if (alignLeft)
      m_horizAlign = align_left;
  }

  void SetVertAlign( const vert_align a ) { m_vertAlign = a; }
  void SetHorizAlign( const horiz_align a ) { m_horizAlign = a; }

  virtual ~text() {}

  void Draw( display_type *d )
    {
      d->SetFont(m_fontName,m_fontSize);
      d->SetColor(m_color);
      d->DrawText(m_coords,m_chars,m_rotAngle,m_fontSize,m_vertAlign,m_horizAlign);
    }

 protected:
  xy_coords m_coords;
  string m_chars;
  color m_color;
  double m_fontSize;
  string m_fontName;
  double m_rotAngle;
  vert_align m_vertAlign;
  horiz_align m_horizAlign;
};




class decorator : public graphic
{
 public:
  decorator(graphic *g) : mp_graphic(g) {}
  
  virtual void Draw( display_type *d )
    { 
      mp_graphic->Draw(d);
    }

 protected:
  graphic *mp_graphic;

public:
  virtual ~decorator(){}
};



class ps_display : public display_type
{
 public:
  ps_display( ostream &ostrm,
	      double horizSize, double vertSize,
	      double border = 0.0 )
    : m_outstrm(ostrm),
      m_color(white),
      m_lineWidth(0.0)
    {
      PrintEpsHeader( m_outstrm, horizSize, vertSize, border ); 
    }

  virtual ~ps_display() {}
  
  void DrawPoint( const xy_coords& coords)
    {
      m_outstrm <<coords.first-m_lineWidth/2<<" "<<coords.second<<" moveto" << "\n";
      m_outstrm <<m_lineWidth<<" 0 rlineto"<< "\n";
      m_outstrm <<"stroke" << "\n";
    }

  void DrawPoint( double x, double y)
    {
      m_outstrm << x - m_lineWidth/2<<" "<< y <<" moveto" << "\n";
      m_outstrm << m_lineWidth << " 0 rlineto"<< "\n";
      m_outstrm << "stroke" << "\n";
    }

  void DrawLine( const xy_coords& startCoords,
		 const xy_coords& stopCoords )
    {
      m_outstrm <<startCoords.first<<" "<<startCoords.second<<" moveto" << "\n";
      m_outstrm <<stopCoords.first<<" "<<stopCoords.second<<" lineto" << "\n";
      m_outstrm <<"stroke" << "\n";
    }

  void DrawText( const xy_coords& coords,
		 const string& chars,
		 double rot_angle,
                 double point_size,
                 const vert_align vAlign,
                 const horiz_align hAlign )
    {
      m_outstrm <<"("<<chars<<") newpath" << "\n";
      m_outstrm <<coords.first<<" "<<coords.second<<" moveto" << "\n";

      if( rot_angle != 0 )
	m_outstrm << "gsave " << rot_angle << " rotate ";

      switch( hAlign ) {
        case align_left:
          // do nothing
          break;
        case align_center:
          m_outstrm <<"("<<chars<<") stringwidth pop 2 div neg 0 rmoveto ";
          break;
        case align_right:
          m_outstrm <<"("<<chars<<") stringwidth pop neg 0 rmoveto ";
          break;
      }

      switch ( vAlign ) {
        case align_top:
          // move down 80% of the font size
          m_outstrm << "0 " << point_size << " 0.8 mul neg rmoveto ";
          break;
        case align_middle:
          // move down 40% of the font size
          m_outstrm << "0 " << point_size << " 0.4 mul neg rmoveto ";
          break;
        case align_bottom:
          // do nothing
          break;
      }

      m_outstrm << "show ";

      if( rot_angle != 0 )
	m_outstrm << " grestore";

      m_outstrm << "\n";
    }

  void DrawArc( const xy_coords &coords,
		const double radius,
		const double startAngle,
		const double stopAngle )
    {
      m_outstrm <<"newpath" << "\n";
      m_outstrm << coords.first <<" "<<coords.second<<" "<<radius<<" "<<startAngle<<" "
		<< stopAngle<<" arc"<< "\n";
      m_outstrm <<"stroke" << "\n";
    }

  void DrawRect( const xy_coords &topLeftCoords,
                 const xy_coords &bottomRightCoords )
    {
      m_outstrm <<"newpath" << "\n";
      m_outstrm << topLeftCoords.first <<" "<< topLeftCoords.second << " moveto" << "\n";
      m_outstrm << bottomRightCoords.first <<" "<< topLeftCoords.second << " lineto" << "\n";
      m_outstrm << bottomRightCoords.first <<" "<< bottomRightCoords.second << " lineto" << "\n";
      m_outstrm << topLeftCoords.first <<" "<< bottomRightCoords.second << " lineto" << "\n";
      m_outstrm << topLeftCoords.first <<" "<< topLeftCoords.second << " lineto" << "\n";
      m_outstrm <<"fill" << "\n";
    }
		

  const color& GetColor() const { return m_color; }

  void SetColor( const color& C )
    {
      if ( C != m_color )
	{
	  m_outstrm <<C.R()<<" "<<C.G()<<" "<<C.B()<<" setrgbcolor" << "\n";
	  m_color = C;
	}
    }

  double GetLineWidth() const { return m_lineWidth; }

  void SetLineWidth( const double width )
    {
      if ( width != m_lineWidth )
	{
	  m_outstrm << width << " setlinewidth" << "\n";
	  m_lineWidth = width;
	}
    }

  const string& GetFontName() const { return m_fontName; }
  double         GetFontSize() const { return m_fontSize; }
  
  void SetFont( const string& font, const double size )
    {
      if ( font != m_fontName || size != m_fontSize )
	{
	  m_outstrm << "/"<<font<<" findfont "<<size<<" scalefont setfont" << "\n";
	  m_fontName = font;
	  m_fontSize = size;
	}
    }

  void Finish() { m_outstrm << flush; }

 protected:
  ostream& m_outstrm;
  color m_color;
  double m_lineWidth;
  string m_fontName;
  double m_fontSize;
};

}

#endif
