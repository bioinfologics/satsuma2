#define FORCE_DEBUG

#include "visual/Whiteboard.h"

#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/SVector.h"
#include "visual/Color.h"

#include "analysis/SequenceMatch.h"
#include "visual/Axes.h"
#include "analysis/DNAVector.h"

class Chr
{
public:
  Chr() {
    m_y = 0.;  
    m_len = 0;
    m_width = 10.;
    m_res = 4000;
  }

  Chr(double y, double l, const string &c, double w, color col) {
    m_res = 10000;
    m_y = y;
    m_name = c;
    m_len = l;
    m_width = w;
    m_color = col;
    int n = (int)(l/m_res+0.5);
    m_dens.resize(n, 0.);
  }

  int isize() const {return m_dens.isize();}
  //const double operator[] (int i) const {return m_dens[i]/m_res;}
  

  void Set(double y, double l, const string &c, double w = 10.) {
    m_y = y;
    m_name = c;
    m_len = l;
    m_width = w;
    int n = (int)(l/m_res+0.5);
    m_dens.resize(n, 0.);
 }

  

  const string & Name() const {return m_name;}
  color Color() const {return m_color;}
  
  void Add(int from, int to) {
    //cout << "in Add" << endl;
    int index = (int)((double)from/m_res + 0.5);
    //cout << "Set " << from << " - " << to << " -> " << index << " of " <<  m_dens.isize() << endl;
    //cout << ": " << m_dens[index] << endl;
    if (index >= m_dens.isize())
      return;
    m_dens[index] += to - from;
  }

  
  void Draw(ns_whiteboard::whiteboard & board, double scale, double x_offset, const vecDNAVector & dna) {
    cout << "Draw!" << endl;
    int i;
    
    for (i=0; i<m_dens.isize(); i++) {
      double c = 0.999 - m_dens[i]/m_res;
      if (c < 0)
	c = 0;
      color col(c, c, c);
      double from = i*m_res;
      double to = (i+1)*m_res;
      
      cout << i << " " << m_dens[i] << " " << c << endl;

      board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset + from/scale, m_y + m_width ), 
					  ns_whiteboard::xy_coords(x_offset + to/scale, m_y),
					  col) );
    }
    int index = dna.NameIndex(m_name);
    if (index == -1) {
      char tmp[512];
      strcpy(tmp, m_name.c_str());
      tmp[strlen(tmp)-1] = 0;
      index = dna.NameIndex(tmp);
    }
    cout << "Plotting " << m_name << " " << index << endl;
    const DNAVector & d = dna[index];
    double start = 0;
    double stop = 0;
    for (i=1; i<d.isize(); i++) {
      if (d[i-1] != 'N' && d[i] == 'N')
	start = i;
      if (d[i-1] == 'N' && d[i] != 'N') {
	stop = i;
	start /= 10.;
	stop /= 10.;
	board.Add( new ns_whiteboard::rect( ns_whiteboard::xy_coords(x_offset + start/scale, m_y + m_width ), 
					    ns_whiteboard::xy_coords(x_offset + stop/scale, m_y),
					    color(0.99, 0.54, 0.0)) );
	
      }
     
    }
    
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
  svec<double> m_dens;
  double m_res;
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
	if (parser.AsFloat(1) > 1000000) {
	  chr.push_back(Chr(off, parser.AsFloat(1), parser.AsString(0), 10., GetColor(i)));
	  off += 14.;
	}
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
  commandArg<string> fString("-f","fasta file");
  commandArg<string> aStringO("-o","outfile (post-script)");
  commandArg<double> dotSize("-d","dot size", 1.);
  commandArg<double> aScale("-s","scale", 60000.);
  commandArg<int> cTarget("-t","target id", -1);
  commandArg<bool> bF("-f","forward only", 0);
  
  
  commandLineParser P(argc,argv);
  P.SetDescription("Comparative cromosome painter (density)");

  P.registerArg(aStringI1);
  P.registerArg(fString);
  P.registerArg(aStringO);
  P.registerArg(dotSize);
  P.registerArg(aScale);
  P.registerArg(cTarget);
  P.registerArg(bF);

  P.parse();

  string i1 = P.GetStringValueFor(aStringI1);
  string fasta = P.GetStringValueFor(fString);
  string o = P.GetStringValueFor(aStringO);
  double dd = P.GetDoubleValueFor(dotSize);
  int targetID = P.GetIntValueFor(cTarget);
  bool fwOnly = P.GetBoolValueFor(bF);


  int i, j;


  vecDNAVector dna;
  dna.Read(fasta);
  
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

  //d = ReadHeader(name2, query, parser, d + 20.);
  d = ReadHeader(name2, query, parser, x_offset);
  
  y_max = d;

  /*for (i=0; i<target.isize(); i++) {
    d = target[i].DrawBorder(board, scale, x_offset);
    if (d > x_max)
      x_max = d;
      }*/
  for (i=0; i<query.isize(); i++) {
    d = query[i].DrawBorder(board, scale, x_offset);
    if (d > x_max)
      x_max = d;
  }
  

  cout << "Target chr: " << target.isize() << endl;
  cout << "Query chr: " << query.isize() << endl;

  //parser.ParseLine();
  string chrT, chrQ;
  while(parser.ParseLine()) {
    if (parser.GetItemCount() == 10) {
      chrT = parser.AsString(0);
      chrQ = parser.AsString(4);
      continue;
    }
    if (parser.GetItemCount() == 0) {
      continue;
    }
    //const string & nT = parser.AsString(0);
    //const string & nQ = parser.AsString(4);
	
    int fromT = parser.AsInt(0); 
    int toT = parser.AsInt(1); 
    int fromQ = parser.AsInt(3); 
    int toQ = parser.AsInt(4); 

    int iT = FindChr(target, chrT);
    int iQ = FindChr(query, chrQ);
    
    //cout << "Now: " << iT << " " << iQ << "\t" << fromQ << " " << toQ << endl;

    //cout << "Doing something/" << endl;
    
    //target[iT].Add(fromT, toT);
    //cout << "Call Add " << iQ << " " << query.isize() << endl;
    if (iQ != -1)
      query[iQ].Add(fromQ, toQ);

  }

  //cout << "Start drawing." << endl;
  /*
  for (i=0; i<target.isize(); i++) {
    target[i].Draw(board, scale, x_offset);
    }*/
  cout << "Start draw." << endl;
  for (i=0; i<query.isize(); i++) {
    query[i].Draw(board, scale, x_offset, dna);
  }

  ofstream out (o.c_str());
  ns_whiteboard::ps_display display(out, x_max + x_offset, y_max + x_offset);
  board.DisplayOn(&display);
 
  //board.DeletePointers();
 
  return 0;
}
