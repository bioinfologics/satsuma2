
#include "analysis/SatsumaAlign.h"
#include "base/FileParser.h"


void SAlignVec::Read(const string & fileName)
{
  int chunk = 65536;

  FlatFileParser parser;
  
  parser.Open(fileName);

  while (parser.ParseLine()) {
    if (m_count >= m_aligns.isize()) {
      m_aligns.resize(m_count + chunk);
    }
    m_aligns[m_count].Set(parser.AsString(0),
			  parser.AsString(3),
			  parser.AsInt(1),
			  parser.AsInt(2),
			  parser.AsInt(4),
			  parser.AsInt(5),
			  parser.AsFloat(6),
			  parser.AsString(7));

    m_count++;
    
  }
  m_aligns.resize(m_count);
}

void SAlignVec::Write(const string & fileName)
{

  FILE * pOut = fopen(fileName.c_str(), "w");

  for (int i=0; i<m_count; i++) {
    const SAlign & s = m_aligns[i];
    
    fprintf(pOut, "%s\t%d\t%d\t%s\t", s.Target().c_str(), s.StartTarget(), s.EndTarget(), s.Query().c_str());
    fprintf(pOut, "%d\t%d\t%f\t", s.StartQuery(), s.EndQuery(), s.Identity());
 
    if (s.Direction() == -1) {
      fprintf(pOut, "-\n");
    } else {
      fprintf(pOut, "+\n");
    }
  }


  fclose(pOut);
  
}
