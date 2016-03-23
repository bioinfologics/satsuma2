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
    m_dir = 0;
    m_end = 0;
    m_scaff = -1;
  }

  Position(const string & chr, int pos, int dir) {
    Set(chr, pos, dir);
  } 

  void Set(const string & chr, int pos, int dir) {
    m_chr = chr;
    m_pos = pos;
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

  void Add(const string & chr, int pos, int ori) {
    m_all.push_back(Position(chr, pos, ori));
  }
  
  void LoHi() {
    int i, j;
    int minSame = 5;
    //cout << m_name << endl;
    for (i=0; i<m_all.isize(); i++) {
      AddCount(m_all[i].Chr());
    }
    int min = m_all.isize() / 5;
    for (i=1; i<m_all.isize()-1; i++) {
      if (Count(m_all[i].Chr()) < min)
	continue;      


      //if (m_name == "scaffold90_33.9") {
      //cout << m_all[i].Chr() << " pos " << m_all[i].Pos() << endl;
      //}
      if (m_lo.Chr() == "") {
	bool b = false;
	if (i+minSame <= m_all.isize())
	  b = true;
	if (b) {
	  for (j=i+1; j<i+minSame; j++) {
	    if (m_all[i].Dir() != m_all[j].Dir()) {
	      b = false;
	      break;
	    }
	  }
	}
	if (b)
	  m_lo = m_all[i];
      }
      bool bL = false;
      if (i > minSame) 
	bL= true;
      
      if (bL) {
	for (j=i-minSame; j<i; j++) {
	  if (m_all[i].Dir() != m_all[j].Dir()) {
	    bL = false;
	    break;
	  }
	}	

      }

      if (bL)
	m_hi = m_all[i];
    }
    m_all.clear();
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
  P.SetDescription("Merges scaffolds according to a synteny map.");
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

  DNAVector spacer;
  spacer.resize(100);
  for (i=0; i<spacer.isize(); i++)
    spacer[i] = 'n';

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    int dir = 1;
    if (parser.AsString(7) == "-")
      dir = -1;
    if (parser.AsString(0) == last) {
      scaff[scaff.isize()-1].Add(parser.AsString(3), parser.AsInt(4), dir);
    } else {
      Scaffold s;
      s.SetName(parser.AsString(0));
      s.Add(parser.AsString(3), parser.AsInt(4), dir);
      scaff.push_back(s);
      last = parser.AsString(0);
    }
  }
  

  svec<Position> pos;
  for (i=0; i<scaff.isize(); i++) {
    scaff[i].LoHi();
    Position a = scaff[i].Lo();
    a.SetScaff(i, 0);
    if (a.Chr() != "") {
      pos.push_back(a);
      //cout << "Adding: " << a.Chr() << " " << a.Pos() << " " << scaff[i].Name() << " A" << endl;
    }
    Position b = scaff[i].Hi();
    b.SetScaff(i, 1);
    if (b.Chr() != "") {
      pos.push_back(b);
      //cout << "Adding: " << b.Chr() << " " << b.Pos() << " " << scaff[i].Name() << " B" << endl;
     
    }
  }

  Sort(pos);

  Position lastPos;
  last = "";

  svec<int> used;
  used.resize(scaff.isize(), 0);
  int supercount = 0;
  for (i=0; i<pos.isize(); i++) {
    cout << "Checking: " << pos[i].Chr() << " " << pos[i].Pos() << " " << scaff[pos[i].Scaffold()].Name() << endl;
    if (last != pos[i].Chr()) {
      last =  pos[i].Chr();
      lastPos = pos[i];
      continue;
    }
 

    if (lastPos.Scaffold() !=  pos[i].Scaffold()) {
      cout << "Connecting: " << scaff[lastPos.Scaffold()].Name() << "\t" << lastPos.Chr() << "\t" << lastPos.Pos() << "\t" << lastPos.Dir() << endl;
      cout << "        to: " << scaff[pos[i].Scaffold()].Name() << "\t" << pos[i].Chr() << "\t" << pos[i].Pos() << "\t" << pos[i].Dir() << endl;
      cout << "          : " << lastPos.End() << "\t<->\t" << pos[i].End() << endl << endl;      

      if (lastPos.Dir() != pos[i].Dir()) {
	//cout << "Orientation mismatch, not merging." << endl;
	//continue;
      }

      if (lastPos.End() == 1) {
	if (pos[i].End() == 0) {
	  scaff[lastPos.Scaffold()].SetRight(pos[i].Scaffold(), pos[i].Dir(), 0);
	  scaff[pos[i].Scaffold()].SetLeft(lastPos.Scaffold(), lastPos.Dir(), 1);
	} else {
	  scaff[lastPos.Scaffold()].SetRight(pos[i].Scaffold(), pos[i].Dir(), 1);
	  scaff[pos[i].Scaffold()].SetRight(lastPos.Scaffold(), lastPos.Dir(), 1);
	}
      } else {
	if (pos[i].End() == 0) {
	  scaff[lastPos.Scaffold()].SetLeft(pos[i].Scaffold(), pos[i].Dir(), 0);
	  scaff[pos[i].Scaffold()].SetLeft(lastPos.Scaffold(), lastPos.Dir(), 0);
	} else {
	  scaff[lastPos.Scaffold()].SetLeft(pos[i].Scaffold(), pos[i].Dir(), 1);
	  scaff[pos[i].Scaffold()].SetRight(lastPos.Scaffold(), lastPos.Dir(), 0);
	}

      }

    }
    lastPos = pos[i];
  }

  for (int i=0; i<scaff.isize(); i++) {
    if (used[i] != 0)
      continue;
    svec<int> left;
    left.push_back(i);
    used[i] = 1;
    int ext = scaff[i].Left();
    int end = scaff[i].LeftEnd();


    cout << "Final Extend: " << scaff[i].Name() << endl;

    while (ext != -1 && used[ext] == 0) {
      cout << "  <- " << scaff[ext].Name() << endl;
      if (end == 0)
	scaff[ext].Flip();
      left.push_back(ext);
      used[ext] = 1;
      end = scaff[ext].LeftEnd();
      ext = scaff[ext].Left();
    }
    svec<int> right;
   
    ext = scaff[i].Right();
    end = scaff[i].RightEnd();
    while (ext != -1 && used[ext] == 0) {
      cout << "  -> " << scaff[ext].Name() << endl;
      //cout << "      -> " << scaff[ext].Right() << endl;
      //cout << "      <- " << scaff[ext].Left() << endl;
      
     if (end == 1)
	scaff[ext].Flip();
      right.push_back(ext);
      used[ext] = 1;
      end = scaff[ext].RightEnd();
      ext = scaff[ext].Right();
      //cout << "      -> ext: " << ext << endl;
    }
    cout << "Superscaffold: " << endl;

    svec<int> toJoin;
    for (i=left.isize()-1; i>=0; i--) {
      //scaff[left[i]].Print();
      toJoin.push_back(left[i]);
    }
    for (i=0; i<right.isize(); i++) {
      //scaff[right[i]].Print();      
      toJoin.push_back(right[i]);
    }
    for (i=0; i<toJoin.isize(); i++) {
      scaff[toJoin[i]].Print();      
    }

    if (dna.isize() > 0) {
      DNAVector super;
      for (i=0; i<toJoin.isize(); i++) {
	scaff[toJoin[i]].Print();      
	DNAVector & dd = dna(scaff[toJoin[i]].Name());
	if (scaff[toJoin[i]].Dir() == -1)
	  dd.ReverseComplement();
	if (super.isize() > 0)
	  super += spacer;
	super += dd;
	dd.resize(0);
      }
      string supername = ">Superscaffold_" + Stringify(supercount);
      supercount++;
      out.push_back(super, supername);
    }
  }


  for (i=0; i<dna.isize(); i++) {
    const DNAVector & dd = dna[i];
    if (dd.isize() > 0)
      out.push_back(dd, dna.Name(i));
  }
  
  out.Write(outName);

  return 0;
}
