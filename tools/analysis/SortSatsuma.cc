#include <string>
#include "../base/CommandLineParser.h"
#include "../base/FileParser.h"


class Line
{
public:
  Line() {}

  void Set(const string & s,
	   const string & chr,
	   int start) {
    m_start = start;
    m_chr = chr;
    m_line = s;
  }


  bool operator < (const Line & l) const {
    if (m_chr != l.m_chr)
      return (m_chr < l.m_chr);
    return (m_start < l.m_start);
  }

  const string & GetLine() const {return m_line;}

private:
  string m_line;
  string m_chr;
  int m_start;
  

};


int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Sorts a Satsuma file by chr & start.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  svec<Line> line;

  int k = 0;

  while (parser.ParseLine()) {
    if (line.isize() <= k)
      line.resize(k+2000000);
    line[k].Set(parser.Line(), parser.AsString(0), parser.AsInt(1));
    k++;
        
  }

  line.resize(k);

  Sort(line);
  for (int i=0; i<line.isize(); i++)
    cout << line[i].GetLine() << endl;
  return 0;
}
