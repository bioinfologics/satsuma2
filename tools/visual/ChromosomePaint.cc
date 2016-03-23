#include "Whiteboard.h"

#include "../base/CommandLineParser.h"
#include "../base/FileParser.h"
#include "../base/SVector.h"
#include "Color.h"

#include "../analysis/SequenceMatch.h"
#include "Axes.h"


class Chr
{
public:
  Chr() {
    m_y = 0.;  
    m_len = 0;
    m_width = 10.;
  }

  Chr(double y, double l, const string &c, double w, color col) {
    m_y = y;
    m_name = c;
    m_len = l;
    m_width = w;
    m_color = col;
  }

  void Set(double y, double l, const string &c, double w = 10.) {
    m_y = y;
    m_name = c;
    m_len = l;
    m_width = w;
  }

  

  const string & Name() const {return m_name;}
  color Color() const {return m_color;}
  
  void DrawFill(ns_whiteboard::whiteboard & board, double scale, double x_offset, double from, double to, color col) {
    //cout << "Filling." << endl;
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset + from/scale, m_y + m_width ), 
                                        ns_whiteboard::xy_coords(x_offset + to/scale, m_y),
                                        col) );

  }

  double DrawBorder(ns_whiteboard::whiteboard & board, double scale, double x_offset) {
    double w = 1.;
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset-w, m_y + m_width + w), 
                                        ns_whiteboard::xy_coords(x_offset + m_len/scale+w, m_y - w),
                                        color(0, 0, 0)) );
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset, m_y + m_width), 
                                        ns_whiteboard::xy_coords(x_offset + m_len/scale, m_y),
                                        color(0.99, 0.99, 0.99)) );
    
    board.Add( new ns_whiteboard::text( ns_whiteboard::xy_coords(x_offset + m_len/scale + 8, m_y + m_width/2),
					m_name, color(0,0,0), 8., "Times-Roman", 0, true));
 
    
    double shift = 13;
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset-w-1.5*m_width-shift, m_y + m_width + w), 
					ns_whiteboard::xy_coords(x_offset+w-1.5*m_width-shift + m_width, m_y - w),
					color(0, 0, 0)) );
    board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset-1.5*m_width-shift, m_y + m_width), 
					ns_whiteboard::xy_coords(x_offset-1.5*m_width-shift + m_width, m_y),
					m_color) );

   return x_offset + m_len/scale;

  }

private:
  double m_y;
  double m_len;
  string m_name;
  double m_width;

  color m_color;
  
};

color GetColor(int i) 
{
  //cout << "GetColor" << endl;
  return MakeUpColor(i);

  double r, g, b;
  r = g = b = 0.;
  
  int v = i % 13;

  switch (v) {
  case 1:
    r = 1.;
    g = 0.;
    b = 0.;
    break;
  case 2:
    r = 0.;
    g = .8;
    b = 0.;
    break;
  case 3:
    r = 0.;
    g = 0.;
    b = .8;
    break;
  case 4:
    r = .7;
    g = 0.7;
    b = 0.;
    break;
  case 5:
    r = 0.7;
    g = 0.;
    b = 0.7;
    break;
    
  case 6:
    r = 0.5;
    g = 0.5;
    b = 0.5;
    break;
  case 9:
    r = 0.95;
    g = 0.5;
    b = 0.1;
    break;
  case 8:
    r = 0.99;
    g = 0.7;
    b = 0.7;
    break;
  case 7:
    r = 0.4;
    g = 0.4;
    b = 0.99;
    break;
    //=============================
  case 10:
    r = 0.99;
    g = 0.4;
    b = 0.0;
    break;
  case 11:
    r = 0.85;
    g = 0.85;
    b = 0.85;
    break;
 case 12:
    r = 0.99;
    g = 0.8;
    b = 0.8;
    break;

    
  default:
    r = g = b = 0.;
    break;
  }

  return color(r, g, b);
}


double ReadHeader(string & name, svec<Chr> &chr, FlatFileParser & parser, double off)
{
  //cout << "Off: " << off << endl;
  while(parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.GetItemCount() == 2 && parser.AsString(0) == "genome") {
      name == parser.AsString(1);
      continue;
    }
    if (parser.GetItemCount() == 2 && parser.AsString(0) == "chromosomes") {
      int n = parser.AsInt(1);
      for (int i=0; i<n; i++) {
	parser.ParseLine();
	//cout << off << endl;
	chr.push_back(Chr(off, parser.AsFloat(1), parser.AsString(0), 10., GetColor(i)));
	off += 14.;
      }
      return off;
    }
  }

  return off;
}
 

int FindChr(const svec<Chr> & c, const string & name) 
{
  for (int i=0; i<c.isize(); i++) {
    if (c[i].Name() == name)
      return i;
  }
  
  return -1;
}
 
int main( int argc, char** argv )
{
 
  commandArg<string> aStringI1("-i","MizBee file");
  commandArg<string> aStringO("-o","outfile (post-script)");
  commandArg<double> dotSize("-d","dot size", 1.);
  commandArg<double> aScale("-s","scale", 60000.);
  commandArg<int> cTarget("-t","target id", -1);
  commandArg<bool> Det("-d","print indivisual matchs", false);
  commandArg<bool> bF("-f","forward only", 0);

  
  commandLineParser P(argc,argv);
  P.SetDescription("Comparative chromosome painter");

  P.registerArg(aStringI1);
  P.registerArg(aStringO);
  P.registerArg(dotSize);
  P.registerArg(aScale);
  P.registerArg(cTarget);
  P.registerArg(Det);
  P.registerArg(bF);

  P.parse();

  string i1 = P.GetStringValueFor(aStringI1);
  string o = P.GetStringValueFor(aStringO);
  double dd = P.GetDoubleValueFor(dotSize);
  int targetID = P.GetIntValueFor(cTarget);
  bool fwOnly = P.GetBoolValueFor(bF);
  bool bDet = P.GetBoolValueFor(Det);


  int i, j;

  
  ns_whiteboard::whiteboard board;

  int x_offset = 100;
  int y_offset = 100;

  double scale = P.GetDoubleValueFor(aScale);
  double x_scale = scale;
  double x_max = 0;
  double y_max = 0;

  FlatFileParser parser;

  parser.Open(i1);
  
  svec<Chr> target;
  svec<Chr> query;

  string name1, name2;
  double d = ReadHeader(name1, target, parser, x_offset);

  d = ReadHeader(name2, query, parser, d + 20.);
  
  y_max = d;

  for (i=0; i<target.isize(); i++) {
    d = target[i].DrawBorder(board, scale, x_offset);
    if (d > x_max)
      x_max = d;
  }
  for (i=0; i<query.isize(); i++) {
    d = query[i].DrawBorder(board, scale, x_offset);
    if (d > x_max)
      x_max = d;
  }
  

  cout << "Target chr: " << target.isize() << endl;
  cout << "Query chr: " << query.isize() << endl;

  string nT, nQ; 

  //parser.ParseLine();
  while(parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;

    if (!bDet && parser.GetItemCount() < 10)
      continue;

    if (parser.GetItemCount() >= 10) {
      nT = parser.AsString(0);
      nQ = parser.AsString(4);
	
      if (!bDet) {
	int fromT = parser.AsInt(1); 
	int toT = parser.AsInt(2); 
	int fromQ = parser.AsInt(5); 
	int toQ = parser.AsInt(6); 
	
	int iT = FindChr(target, nT);
	int iQ = FindChr(query, nQ);
	
	//cout << nT << " " << iT << "\t" << nQ << " " << iQ << endl;
	
	//cout << "Doing something/" << endl;
	target[iT].DrawFill(board, scale, x_offset, fromT, toT,  query[iQ].Color());
	query[iQ].DrawFill(board, scale, x_offset, fromQ, toQ,  target[iT].Color());
      }
    } else {
      if (bDet) {
	int fromT = parser.AsInt(0); 
	int toT = parser.AsInt(1); 
	int fromQ = parser.AsInt(3); 
	int toQ = parser.AsInt(4); 
	
	int iT = FindChr(target, nT);
	int iQ = FindChr(query, nQ);
	
	//cout << nT << " " << iT << "\t" << nQ << " " << iQ << endl;
	
	//cout << "Doing something/" << endl;
	target[iT].DrawFill(board, scale, x_offset, fromT, toT,  query[iQ].Color());
	query[iQ].DrawFill(board, scale, x_offset, fromQ, toQ,  target[iT].Color());

      }
    }
  }

  ofstream out (o.c_str());
  ns_whiteboard::ps_display display(out, x_max + x_offset, y_max + x_offset);
  board.DisplayOn(&display);
 
  //board.DeletePointers();
 
  return 0;
}
