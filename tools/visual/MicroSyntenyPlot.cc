#include "Whiteboard.h"

#include "../base/CommandLineParser.h"
#include "../base/FileParser.h"
#include "../base/SVector.h"
#include "Color.h"

#include "../analysis/SequenceMatch.h"
#include "Axes.h"


class ToPlot
{
public:
  ToPlot() {
    m_x1 = m_x2 = m_y1 = m_y2 = 0;;
    m_ident = 0.;
    m_q = 0;
    m_ta = 0;
  }


  void Set(int x1, int y1, int x2, int y2, double ident, int q, int t)
  {
    m_x1 = x1;
    m_y1 = y1;
    m_x2 = x2;
    m_y2 = y2;
    m_ident = ident;
    m_q = q;
    m_ta = t;
  }

  int X1() const {return m_x1;}
  int Y1() const {return m_y1;}
  int X2() const {return m_x2;}
  int Y2() const {return m_y2;}

  double Ident() const {return m_ident;}
  int Query() const {return m_q;}
  int Target() const {return m_ta;}

  bool operator < (const ToPlot & t) const {
    return (m_ident < t.m_ident);
  }

private:
  int m_x1, m_x2, m_y1, m_y2;
  double m_ident;
  int m_q;
  int m_ta;
};


double Intensity(double val, double strength)
{

  val += strength;
  if (val > 1.)
    val = 1.;
  if (val < 0.)
    val = 0.;
  return val;
}

double IntensityVal(double val)
{


  double ret = 1.0 - 2. * (val - 0.5);;

  if (ret > 0.9999)
    ret = 0.9999;
  if (ret < 0.)
    ret = 0.;

  //cout << "Translated " << val << " into " << ret << endl;

  return ret;
}


void GetUniqueList(svec<int> & unique, const MultiMatches & multi, bool bForce) 
{
  int i, j;
  
  int lastTarget = -1;
  svec<int> count;
  svec<int> qID;

  unique.resize(multi.GetMatchCount(), 0);

  if (bForce) {
    for (i=0; i<multi.GetMatchCount(); i++) {
      unique[i] = 1;
    }
    return;
  }


  for (i=0; i<multi.GetMatchCount(); i++) {
    const SingleMatch & m = multi.GetMatch(i);

    if (i % 1000000 == 0)
      cout << i << endl;

    if (m.GetTargetID() != lastTarget || i+1 == multi.GetMatchCount()) {
      for (j=0; j<count.isize(); j++) {
	if (count[j] == 1)
	  unique[qID[j]] = 1;
      }

      count.clear();
      qID.clear();
      lastTarget = m.GetTargetID();
      count.resize(multi.GetTargetSize(lastTarget), 0);
      qID.resize(multi.GetTargetSize(lastTarget), -1);      
    }

    int start = m.GetStartTarget();
    int end = m.GetStartTarget() + m.GetLength();
    if (m.GetProbability() > 0.9999) {
      for (j=start; j<end; j++) {
	count[j]++;
	qID[j] = i;
      }
    }
  }
  
}


int main( int argc, char** argv )
{
 
  commandArg<string> aStringI1("-i","HomologyByXCorr output file");
  commandArg<string> aStringO("-o","outfile (post-script)");
  commandArg<double> dotSize("-d","dot size", 1.);
  commandArg<double> aScale("-s","scale", 60000.);
  commandArg<int> cTarget("-t","target id", -1);
  commandArg<bool> bF("-f","forward only", 0);

  
  commandLineParser P(argc,argv);
  P.SetDescription("Micro-synteny plotter");

  P.registerArg(aStringI1);
  P.registerArg(aStringO);
  P.registerArg(dotSize);
  P.registerArg(aScale);
  P.registerArg(cTarget);
  P.registerArg(bF);

  P.parse();

  string i1 = P.GetStringValueFor(aStringI1);
  string o = P.GetStringValueFor(aStringO);
  double dd = P.GetDoubleValueFor(dotSize);
  int targetID = P.GetIntValueFor(cTarget);
  bool fwOnly = P.GetBoolValueFor(bF);


  int i, j;


  MultiMatches multi;
  multi.Read(i1);
  cout << "Done loading." << endl;


  cout << "Using dot size " << dd << endl;
  //cout << "Sorting..." << endl;
  //multi.Sort();
  //cout << "done." << endl;

  int x_offset = 100;
  int y_offset = 100;

  double scale = P.GetDoubleValueFor(aScale);
  double x_scale = scale;
  double x_max = 0;
  double y_max = 0;


  
  for (i=0; i<multi.GetMatchCount(); i++) {
    const SingleMatch & m = multi.GetMatch(i);

    //cout << m.GetTargetID() << " " << m.GetQueryID() << endl;

    if (targetID != -1 && m.GetTargetID() != targetID)
      continue;

    if (m.GetStartTarget() + m.GetLength() > x_max)
      x_max = m.GetStartTarget() + m.GetLength();
    if (m.GetStartQuery() + m.GetLength() > y_max)
      y_max = m.GetStartQuery() + m.GetLength();
  }

  //int endY = y_max;
  

  x_max /= x_scale;
  y_max /= scale;

  x_max += 2 * x_offset;
  y_max += 2 * y_offset;

  cout << "x_max=" << x_max << "  y_max=" << y_max << endl;

  ns_whiteboard::whiteboard board;
  
  svec<ToPlot> dots;
  dots.resize(multi.GetMatchCount());

  int k = 0;

  double minProb = 0.5;


  cout << "Total matches: " << multi.GetMatchCount() << endl;
  svec<int> unique;
  GetUniqueList(unique, multi, true);

  int n = 0;

  for (i=0; i<multi.GetMatchCount(); i++) {
    //if (unique[i] == 0)
    //continue;

    const SingleMatch & m = multi.GetMatch(i);

    if (targetID != -1 && m.GetTargetID() != targetID) {
      //cout << "ERROR: " << endl;
      //cout << "target=" << m.GetTargetID() << " query=" << m.GetQueryID() << " t=" << m.GetStartTarget() << " q=" << m.GetStartQuery() << " len=" << m.GetLength() << endl;
      continue;
    }

    int x1 = m.GetStartTarget();
    int y1 = m.GetStartQuery();
    int x2 = m.GetStartTarget() + m.GetLength();
    int y2 = m.GetStartQuery() + m.GetLength();
    
    //if (y1 > 16000000)
    //continue;

    //x1 -= 5 * m.GetLength();
    //y1 -= 5 * m.GetLength();
    //x2 += 5 * m.GetLength();
    //y2 += 5 * m.GetLength();


    double ident = m.GetProbability();

    //if (ident < minProb)
    //continue;
    
    if (m.GetLength() == 0) {
      cout << "ERROR: len=0" << endl;
      continue;
    }

    n++;

    if (m.IsRC()) {
      if (fwOnly)
	continue;

      int s = multi.GetQuerySize(m.GetQueryID());
      //if (s < 200000)
      //continue;
      //cout << "Rc'ing... " << endl;
      //cout << "Old: " << y1 << " " << y2 << endl;
      int tmp = y1;
      y1 = s - y1;
      y2 = s - y2;
      //cout << "New: " << y1 << " " << y2 << endl;
    } 
    
    //cout << "Adding ident: " << ident << endl;

    dots[k].Set(x1, y1, x2, y2, ident, m.GetQueryID(), m.GetTargetID());
    
    k++;
  }

  dots.resize(k);

  Sort(dots);

  for (i=0; i<dots.isize(); i++) {
  
    double x1 = dots[i].X1();
    double y1 = dots[i].Y1();
    double x2 = dots[i].X2();
    double y2 = dots[i].Y2();
    
    //if (m.IsRC()) {
     
    //}

    x1 /= x_scale;
    y1 /= scale;
    x2 /= x_scale;
    y2 /= scale;

    //x2 += dd - 1;
    //y2 += dd - 1;
    
    double ident = dots[i].Ident();
    double v = ident;
    
    
    v /= (1. - minProb);
    if (v < 0.)
      v = 0.;
    if (v > 1.)
      v = 1.;
    
    v = 0.;

    //cout << "Ident=" << ident << " v=" << v << endl;

    double r, g, b;
    
    int c = (dots[i].Query() + dots[i].Target()) % 10;
    //cout << "c=" << c << endl;
    r = g = b = v;
    
    switch (c) {
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

    default:
      r = g = b = 0.;
      break;
    }
    
      //cout << "r=" << Intensity(r, v) << " g=" <<  Intensity(g, v) << " b=" <<  Intensity(b, v) << endl;

    //color black(Intensity(r, v), Intensity(g, v), Intensity(b, v));
    color black(r, g, b);

    //color black(IntensityVal(v), IntensityVal(v), IntensityVal(v));
    //color black(r, g, b);
    //cout << "color=" << v << endl;
    board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x1 + x_offset, y1 + y_offset), 
					ns_whiteboard::xy_coords(x2 + x_offset, y2 + y_offset),
					dd, black) );
  }

  /*
  for (i=0; i<20000000; i+= 1000000) {
    int x1 = (double)i / x_scale;
    board.Add( new ns_whiteboard::line( ns_whiteboard::xy_coords(x1 + x_offset, y_offset), 
					ns_whiteboard::xy_coords(x1 + x_offset, y_offset + y_max),
					1, black) );
  }
  */
  ns_whiteboard::line x_axis(make_pair(x_offset,y_offset),
                             make_pair(x_max-x_offset,y_offset),
                             1.0);
  ns_whiteboard::line y_axis(make_pair(x_offset, y_offset),
                             make_pair(x_offset,y_max-y_offset),
                             1.0);
  ns_whiteboard::axes A(&x_axis,&y_axis);

  string xlab("Genome 1");
  string ylab("Genome 2");


  A.CreateTics(scale, 0., scale, 0.);
  
  A.SetXAxisLabel(xlab,
                  red,
                  18);
  A.SetYAxisLabel(ylab,
                  red,
                  18);


  
  board.Add(&A);




  ofstream out (o.c_str());
  ns_whiteboard::ps_display display(out, x_max, y_max);
  board.DisplayOn(&display);
 
  //board.DeletePointers();
 
  return 0;
}
