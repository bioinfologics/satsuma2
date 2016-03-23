#define FORCE_DEBUG

#include <string>
#include "../base/CommandLineParser.h"
#include "../base/FileParser.h"
#include "DNAVector.h"
#include "../base/StringUtil.h"

class Position
{
public:
  Position() {
    m_pos = -1;
    m_posQ = -1;
    m_dir = 0;
    m_end = 0;
    m_scaff = -1;
  }

  Position(const string & chr, int pos, int dir, int posQ) {
    Set(chr, pos, dir, posQ);
  } 

  void Set(const string & chr, int pos, int dir, int posQ) {
    m_chr = chr;
    m_pos = pos;
    m_posQ = posQ;
    m_dir = dir;
  }

  void SetScaff(int i, int end) {
    m_scaff= i;
    m_end = end;
  }

  void Flip() {
    m_dir *= -1;
  }

  const string & Chr() const {return m_chr;}
  int Pos() const {return m_pos;}
  int PosQ() const {return m_posQ;}
  int Dir() const {return m_dir;}
  int Scaffold() {return m_scaff;}
  int End() {return m_end;}

  bool operator < (const Position & pos) const {
    if (m_chr != pos.m_chr)
      return m_chr < pos.m_chr;
    return m_pos < pos.m_pos;
  }

private:
  int m_pos;
  int m_posQ;
  int m_dir;
  string m_chr;
  int m_scaff;
  int m_end;
};

class Scaffold
{
public:
  Scaffold() {
    m_dir = 1;
    m_left = -1;
    m_right = -1;
    m_leftDir = 0;
    m_rightDir = 0;
    m_rightEnd = -1;
    m_leftEnd = -1;
  }

  void SetName(const string & n) {
    m_name = n;
  }

  void Add(const string & chr, int pos, int ori, int posQ) {
    m_all.push_back(Position(chr, pos, ori, posQ));
  }
  
  void LoHi() {
    int i;
    //cout << m_name << endl;
    for (i=0; i<m_all.isize(); i++) {
      AddCount(m_all[i].Chr());
    }
    int min = m_all.isize() / 5;
    for (i=0; i<m_all.isize(); i++) {
      if (Count(m_all[i].Chr()) < min)
	continue;
      //if (m_name == "scaffold90_33.9") {
      //cout << m_all[i].Chr() << " pos " << m_all[i].Pos() << endl;
      //}
      if (m_lo.Chr() == "") {
	m_lo = m_all[i];
      }
      m_hi = m_all[i];
    }
    m_all.clear();
    m_dir = m_lo.Dir();
  }

  int Dir() const {return m_dir;}
  const string & Name() const {return m_name;}
  const Position & Lo() const {return m_lo;} 
  const Position & Hi() const {return m_hi;} 
  int Left() const {return m_left;}
  int Right() const {return m_right;}
  int LeftDir() const {return m_leftDir;}
  int RightDir() const {return m_rightDir;}
  int LeftEnd() const {return m_leftEnd;}
  int RightEnd() const {return m_rightEnd;}


  void Flip() {
    m_dir *= -1;
    Position tmp = m_lo;
    m_lo = m_hi;
    m_hi = tmp;
    m_lo.Flip();
    m_hi.Flip();
    int t = m_right;
    m_right = m_left;
    m_left = t;
    t = m_rightDir;
    m_rightDir = m_leftDir;
    m_leftDir = t;  
    t = m_rightEnd;
    m_rightEnd = m_leftEnd;
    m_leftEnd = t;  
  }

  void SetRight(int i, int dir, int end) {
    m_right = i;
    m_rightDir = dir;
    m_rightEnd = end;
  }

  void SetLeft(int i, int dir, int end) {
    m_left = i;
    m_leftDir = dir;
    m_leftEnd = end;
  }

  void Print() const {
    cout << m_name << " " << m_dir << endl;
    cout << "Lo: " << m_lo.Chr() << " " << m_lo.Pos() << " " << m_lo.Dir() << endl;
    cout << "Hi: " << m_hi.Chr() << " " << m_hi.Pos() << " " << m_hi.Dir() << endl;
  }

private:
  int Count(const string & chr) const {
    int i;
    for (i=0; i<m_diff.isize(); i++) {
      if (m_diff[i] == chr) {
	return m_count[i];      
      }
    }
    return 0;
  }

  void AddCount(const string & chr) {
    int i;
    for (i=0; i<m_diff.isize(); i++) {
      if (m_diff[i] == chr) {
	m_count[i]++;
	return;
      }
    }
    m_diff.push_back(chr);
    m_count.push_back(1);
  }

  svec<string> m_diff;
  svec<int> m_count;

  string m_name;
  svec<Position> m_all;
  
  Position m_lo;
  Position m_hi;

  int m_dir;

  int m_right;
  int m_left;
  int m_rightDir;
  int m_leftDir;
  int m_rightEnd;
  int m_leftEnd;
};

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file (satsuma summary)");
  commandArg<string> fastaCmmd("-f","fasta file", "");
  commandArg<string> outCmmd("-o","output fasta file", "superscaffolds.fasta");
  commandLineParser P(argc,argv);
  P.SetDescription("Orders and orients scaffolds according to a synteny map.");
  P.registerArg(fileCmmd);
  P.registerArg(fastaCmmd);
  P.registerArg(outCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  

  vecDNAVector dna;
  vecDNAVector out;
  if (fastaName != "") 
    dna.Read(fastaName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  string last;
  svec<Scaffold> scaff;
  int i, j;


  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    int dir = 1;
    if (parser.AsString(7) == "-")
      dir = -1;
    if (parser.AsString(0) == last) {
      scaff[scaff.isize()-1].Add(parser.AsString(3), parser.AsInt(4), dir, parser.AsInt(1));
    } else {
      Scaffold s;
      s.SetName(parser.AsString(0));
      s.Add(parser.AsString(3), parser.AsInt(4), dir, parser.AsInt(1));
      scaff.push_back(s);
      last = parser.AsString(0);
    }
  }
  

  svec<Position> pos;
  cout << "Total scaffolds: " << scaff.isize() << endl;
  for (i=0; i<scaff.isize(); i++) {
    scaff[i].LoHi();
    Position a = scaff[i].Lo();
    a.SetScaff(i, 0);
    if (a.Chr() != "") {
      pos.push_back(a);
      cout << "Adding: " << a.Chr() << " " << a.Pos() << " " << scaff[i].Name() << " A" << endl;
    } else {
      //cout << "Skip " << a.Chr() << " " << a.Pos() << " " << scaff[i].Name() << " A" << endl;
    }
    Position b = scaff[i].Hi();
    b.SetScaff(i, 1);
    if (b.Chr() != "") {
      pos.push_back(b);
      cout << "Adding: " << b.Chr() << " " << b.Pos() << " " << scaff[i].Name() << " B" << endl;
     
    }
  }

  Sort(pos);

  Position lastPos;
  last = "";

  svec<int> used;
  used.resize(scaff.isize(), 0);
  int supercount = 0;
  DNAVector tmp;
  int lastMappedPos = -1;
  for (i=0; i<pos.isize(); i++) {
    if (last != pos[i].Chr()) {
      string name = ">Pseudo";
      name += last;
      if (tmp.isize() > 0) {
	out.push_back(tmp, name);
      }
      last =  pos[i].Chr();     
      tmp.resize(0);
      lastMappedPos = -1;
      //continue;
    } 
    cout << "Checking: " << pos[i].Chr() << " " << pos[i].Pos() << " " << scaff[pos[i].Scaffold()].Dir() << " " << scaff[pos[i].Scaffold()].Name() << endl;
    DNAVector & nn = dna(scaff[pos[i].Scaffold()].Name());
    if (nn.isize() == 0) {
      cout << "Skipping..." << endl;
      continue;
    }

    const Position & theLo = scaff[pos[i].Scaffold()].Lo();
    const Position & theHi = scaff[pos[i].Scaffold()].Hi();
    int gapEst = 0;
    int hiPosQ = theHi.PosQ();
    int loPosQ = theLo.PosQ();
    if (scaff[pos[i].Scaffold()].Dir() == -1) {
      loPosQ = theHi.PosQ();
      hiPosQ = theLo.PosQ();
      //gapEst = 100;
    } 
    if (pos[i].End() == 1) {
      int firstPos = pos[i].Pos() - theHi.PosQ();
      gapEst = firstPos - lastMappedPos;
      lastMappedPos = pos[i].Pos() + nn.isize() - hiPosQ;
      cout << "END gap: " << gapEst << " last: " << lastMappedPos << " len: " << nn.isize() << endl;
      //gapEst = 100;
    } else {
      int firstPos = pos[i].Pos() - theLo.PosQ();
      gapEst = firstPos - lastMappedPos;
      lastMappedPos = pos[i].Pos() + nn.isize() - loPosQ;
      cout << "START gap: " << gapEst << " last: " << lastMappedPos  << " len: " << nn.isize() << endl;
      //gapEst = 100;
    }
    
    

    if (tmp.isize() > 0) {
      DNAVector spacer;
      
      cout << "Estimate gap (raw): " << gapEst << endl;
      if (gapEst < 10)
	gapEst = 10;
      if (gapEst > 1000)
	gapEst = 1000;
      //gapEst = 100;
      cout << "Estimate gap (capped): " << gapEst << endl;
      spacer.resize(gapEst);
      for (int x=0; x<spacer.isize(); x++)
	spacer[x] = 'n';
      tmp += spacer;
    }
    if (scaff[pos[i].Scaffold()].Dir() == -1)
      nn.ReverseComplement();
    tmp += nn;

    
    nn.resize(0);

  }
  cout << "Done here." << endl;

  string name2 = ">Pseudo";
  name2 += last;
  if (tmp.isize() > 0) {
    out.push_back(tmp, name2);
  }

  cout << "Push" << endl;
  for (i=0; i<dna.isize(); i++) {
    const DNAVector & dd = dna[i];
    if (dd.isize() > 0)
      out.push_back(dd, dna.Name(i));
  }

  cout << "Write to " << outName << endl;
  
  out.Write(outName);

  return 0;
}
