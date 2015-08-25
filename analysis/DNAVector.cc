//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif

#include "analysis/DNAVector.h"
#include "base/FileParser.h"
#include "util/mutil.h"
#include "analysis/Coordinate.h"


#define ONE_THIRD 1./3.






DNACodec::DNACodec() 
{
  int i;
  for (i=0; i<256; i++) {
    m_A[i] = m_C[i] = m_G[i] = m_T[i] = 0;
    m_rc[i] = 0;
  }

  Set('A', 1., 0., 0., 0., 'T');
  Set('C', 0., 1., 0., 0., 'G');
  Set('G', 0., 0., 1., 0., 'C');
  Set('T', 0., 0., 0., 1., 'A');



  Set('K', 0., 0., 0.5, 0.5, 'M');
  Set('M', 0.5, 0.5, 0., 0., 'K');

  Set('R', 0.5, 0., 0.5, 0., 'Y');
  Set('Y', 0., 0.5, 0., 0.5, 'R');

  Set('S', 0., 0.5, 0.5, 0., 'S');
  Set('W', 0.5, 0., 0., 0.5, 'W');

  //------------------------------------------
  //Set('K', 0., 0., 0.75, 0.75, 'M');
  //Set('M', 0.75, 0.75, 0., 0., 'K');

  //Set('R', 0.75, 0., 0.75, 0., 'Y');
  //Set('Y', 0., 0.75, 0., 0.75, 'R');

  //Set('S', 0., 0.75, 0.75, 0., 'S');
  //Set('W', 0.75, 0., 0., 0.75, 'W');
  //-----------------------------------------

  Set('B', 0., ONE_THIRD, ONE_THIRD, ONE_THIRD, 'V');
  Set('V', ONE_THIRD,ONE_THIRD , ONE_THIRD, 0., 'B');

  Set('H', ONE_THIRD, ONE_THIRD, 0., ONE_THIRD, 'D');
  Set('D', ONE_THIRD, 0., ONE_THIRD, ONE_THIRD, 'H');

  Set('-', 0., 0., 0., 0., '-');
  Set('N', 0.25, 0.25, 0.25, 0.25, 'N');
  Set('X', 0.25, 0.25, 0.25, 0.25, 'X');

}

char ResolveAmbiguous(const string & bases)
{
  const char * p = (const char*)bases.c_str();
  int n = strlen(p);
  if (n != 2)
    return GetAmbiguous(bases);

  bool bNuke0 = false;
  if (p[0] == 'A' || p[0] == 'C' || p[0] == 'G' || p[0] == 'T') {
    bNuke0 = true;
  }
  bool bNuke1 = false;
  if (p[1] == 'A' || p[1] == 'C' || p[1] == 'G' || p[1] == 'T') {
    bNuke1 = true;
  }

  if (bNuke0 && bNuke1)
    return GetAmbiguous(bases);

  if (bNuke0)
    return p[0];
  
  if (bNuke1)
    return p[1];
  
  
  return p[0];
}



char GetAmbiguous(const string & bases)
{
  bool a = false;
  bool c = false;
  bool g = false;
  bool t = false;

  const char * p = (const char*)bases.c_str();
  int i;
  int n = (int)strlen(p);
  for (i=0; i<n; i++) {
    if (p[i] == 'A')
      a = true;
    if (p[i] == 'C')
      c = true;
    if (p[i] == 'G')
      g = true;
    if (p[i] == 'T')
      t = true;
  }


  if (a && c && g && t)
    return 'N';

  if (a && g && t)
    return 'D';
  if (a && c && t)
    return 'H';
  if (a && c && g)
    return 'V';
  if (c && g && t)
    return 'B';


  if (a && t)
    return 'W';
  if (c && g)
    return 'S';
  if (c && t)
    return 'Y';
  if (a && g)
    return 'R';
  if (a && c)
    return 'M';
  if (g && t)
    return 'K';


  if (a)
    return 'A';
  if (c)
    return 'C';
  if (g)
    return 'G';
  if (t)
    return 'T';

  return '?';
}


void DNACodec::Set(int letter, double a, double c,double g, double t, int rc)
{
  m_rc[letter] = rc;
  m_A[letter] = a;
  m_C[letter] = c;
  m_G[letter] = g;
  m_T[letter] = t;
}


DNACodec theCodec;

char GetRC(const char c)
{
  return theCodec.GetRC(c);
}



AACodons::AACodons()
{
  m_table.resize(256, -1);
  m_bases.resize(25);

  Set("A", "GCT", 0); 
  Set("A", "GCC", 0); 
  Set("A", "GCA", 0); 
  Set("A", "GCG", 0); 
  
  Set("L", "TTA", 1); 
  Set("L", "TTG", 1); 
  Set("L", "CTT", 1); 
  Set("L", "CTC", 1); 
  Set("L", "CTA", 1); 
  Set("L", "CTG", 1);
  
  Set("R", "CGT", 2); 
  Set("R", "CGC", 2); 
  Set("R", "CGA", 2); 
  Set("R", "CGG", 2); 
  Set("R", "AGA", 2); 
  Set("R", "AGG", 2);
  
  Set("K", "AAA", 3); 
  Set("K", "AAG", 3); 
  
  Set("N", "AAT", 4); 
  Set("N", "AAC", 4); 
  
  Set("M", "ATG", 5); 
  
  Set("D", "GAT", 6); 
  Set("D", "GAC", 6); 
  
  Set("F", "TTT", 7); 
  Set("F", "TTC", 7); 
  
  Set("C", "TGT", 8); 
  Set("C", "TGC", 8); 
  
  Set("P", "CCT", 9); 
  Set("P", "CCC", 9); 
  Set("P", "CCA", 9); 
  Set("P", "CCG", 9); 
  
  Set("Q", "CAA", 10); 
  Set("Q", "CAG", 10); 
  
  Set("S", "TCT", 11); 
  Set("S", "TCC", 11); 
  Set("S", "TCA", 11); 
  Set("S", "TCG", 11); 
  Set("S", "AGT", 11); 
  Set("S", "AGC", 11); 
  
  Set("E", "GAA", 12); 
  Set("E", "GAG", 12); 
  
  Set("T", "ACT", 13); 
  Set("T", "ACC", 13); 
  Set("T", "ACA", 13); 
  Set("T", "ACG", 13); 
  
  Set("G", "GGT", 14); 
  Set("G", "GGC", 14); 
  Set("G", "GGA", 14); 
  Set("G", "GGG", 14); 
  
  Set("W", "TGG", 15); 
  
  Set("H", "CAT", 16); 
  Set("H", "CAC", 16); 
  
  Set("Y", "TAT", 17); 
  Set("Y", "TAC", 17); 
  
  
  Set("I", "ATT", 18); 
  Set("I", "ATC", 18); 
  Set("I", "ATA", 18); 
  
  Set("V", "GTT", 19); 
  Set("V", "GTC", 19); 
  Set("V", "GTA", 19); 
  Set("V", "GTG", 19); 
  
  Set("*", "TAG", 20); 
  Set("*", "TGA", 20); 
  Set("*", "TAA", 20); 


  //Special/ambiguous 

  //B (N or D)
  Set("B", "AAT", 21); 
  Set("B", "AAC", 21); 
  Set("B", "GAT", 21); 
  Set("B", "GAC", 21); 


  //Z (Q or E)
  Set("Z", "GAA", 22); 
  Set("Z", "GAG", 22); 
  Set("Z", "CAA", 22); 
  Set("Z", "CAG", 22); 


  //J (L or I)
  Set("J", "ATT", 23); 
  Set("J", "ATC", 23); 
  Set("J", "ATA", 23); 
  Set("J", "TTA", 23); 
  Set("J", "TTG", 23); 
  Set("J", "CTT", 23); 
  Set("J", "CTC", 23); 
  Set("J", "CTA", 23); 
  Set("J", "CTG", 23);
  
 

  //X anything
  Set("X", "AAA", 24); 
  Set("X", "CCC", 24); 
  Set("X", "GGG", 24); 
  Set("X", "TTT", 24); 

}

char AACodons::GetCodon(const DNAVector & d, int pos)
{
  char tmp[64];
  tmp[0] = d[pos];
  tmp[1] = d[pos+1];
  tmp[2] = d[pos+2];
  tmp[3] = 0;
  int i;
  
  //cout << "Asking for " << tmp << endl;
  
  for (i=0; i<m_codonBases.size(); i++) {
    if (m_codonBases[i] == tmp)
      return m_aminoAcids[i];
  }
  //cout << "WEIRD ERROR!" << endl;
  return 'Q';
}





AACodons trans;

//==========================================
void AminoAcidToBases(char * out, char aa)
{

  //cout << "Query " << aa << endl;
  int i;
  for (i=0; i<3; i++) {
    const string & b = trans.GetBases(aa, i);
    out[i] = GetAmbiguous(b);
    if (out[i] == 'N')
      out[i] = 'X';
  }
  out[4] = 0;
}



char BasesToAminoAcid(char * b)
{
  return trans.GetCodon(b);
}




double DNA_A(char l)
{
  return theCodec.A_Letter(l);
}

double DNA_C(char l)
{
  return theCodec.C_Letter(l);
}

double DNA_G(char l)
{
  return theCodec.G_Letter(l);
}

double DNA_T(char l)
{
  return theCodec.T_Letter(l);
}


double DNA_EqualAmb(char letter1, char letter2)
{
  if (letter1 == letter2) {
    if (letter1 != 'N')
      return 1.;
  }
  double a = theCodec.A_Letter(letter1) * theCodec.A_Letter(letter2);  
  double c = theCodec.C_Letter(letter1) * theCodec.C_Letter(letter2);  
  double g = theCodec.G_Letter(letter1) * theCodec.G_Letter(letter2);  
  double t = theCodec.T_Letter(letter1) * theCodec.T_Letter(letter2);  

  double sum = a + c + g + t;
  double plus = 0.;

  return (sum + plus);
}


double DNA_Equal(char letter1, char letter2)
{
  double a = theCodec.A_Letter(letter1) * theCodec.A_Letter(letter2);  
  double c = theCodec.C_Letter(letter1) * theCodec.C_Letter(letter2);  
  double g = theCodec.G_Letter(letter1) * theCodec.G_Letter(letter2);  
  double t = theCodec.T_Letter(letter1) * theCodec.T_Letter(letter2);  

  double sum = a + c + g + t;
  double plus = 0.;
  //if (sum < 1. && letter1 == letter2)
  //plus = 0.25;

  return (sum + plus);
}



double DNA_EqualEmph(char letter1, char letter2)
{
  double v = DNA_Equal(letter1, letter2);
  if (v > 0.26)
    v += 0.2;
  if (v > 1.0)
    v = 1.;
  return v;
}


double DNA_Diff(char letter1, char letter2)
{
  return (1 - DNA_Equal(letter1,letter2));
}


double DNA_DiffAmb(char letter1, char letter2)
{
  return (1-DNA_EqualAmb(letter1,letter2));
}


//====================================================
bool DNAVector::SetToSubOf(const DNAVector & v, int start, int len)
{
  //if(v.size() < start+len) { return false; } //TODO WARN 

  // Let's fix this differently...
  if (len < 0) {
    cout << "WARNING, negative length provided in DNAVector::SetToSubOf -> " << len << endl;
    return false;
  }
  if (start < 0) {
    cout << "WARNING, adjusting boundaries in DNAVector::SetToSubOf (1)" << endl;
    len += start;
    start = 0;
  }
  if (start+len > v.size()) {
    cout << "WARNING, adjusting boundaries in DNAVector::SetToSubOf (2)" << endl;
    len = v.size() - start;
  }

  m_data.resize(len);
  int i;
  for (i=start; i<start+len; i++)
    m_data[i-start] = v[i];

  if (v.QualSize() > 0) {
    m_qual.resize(len);
    for (i=start; i<start+len; i++)
      m_qual[i-start] = v.Qual(i);
    
  }
  return true;
}

bool DNAVector::SetToSubOf(const DNAVector & v, int start)
{
  if(start>v.size()) { return false; } // TODO WARN 
  return SetToSubOf(v, start, v.size()-start);
}

void DNAVector::RemoveGaps() 
{
  int i, k = 0;
  for (i=0; i<size(); i++) {
    if (m_data[i] != '-' && m_data[i] != '.') {
      m_data[k] = m_data[i];
      k++;
    }
  }
  resize(k);
}

void DNAVector::ReverseComplement() {
  int n = m_data.size();
  int i = 0;
  int j = n-1;

  n = (n+1)/2;

  for (i=0; i<n; i++) {
    char one = m_data[i];
    char two = m_data[j];
    one = theCodec.GetRC(one);
    two = theCodec.GetRC(two);

    m_data[i] = two;
    m_data[j] = one;
    
    j--;
    
  }

  if (m_qual.size() > 0) {
    n = m_data.size();
    i = 0;
    j = n-1;

    n = (n+1)/2;

    for (i=0; i<n; i++) {
      char one = m_qual[i];
      char two = m_qual[j];           
      m_qual[i] = two;
      m_qual[j] = one;
      
      j--;
    
    }
  }


}

void DNAVector::Reverse() {
  int n = m_data.size();
  int i = 0;
  int j = n-1;

  n = (n+1)/2;

  for (i=0; i<n; i++) {
    char one = m_data[i];
    char two = m_data[j];

    m_data[i] = two;
    m_data[j] = one;
    
    j--;
    
  }

  if (m_qual.size() > 0) {
    n = m_data.size();
    i = 0;
    j = n-1;

    n = (n+1)/2;

    for (i=0; i<n; i++) {
      char one = m_qual[i];
      char two = m_qual[j];           
      m_qual[i] = two;
      m_qual[j] = one;
      
      j--;
    
    }
  }


}

void DNAVector::SetFromBases(const string & s)
{
  int n = strlen(s.c_str());
  m_data.resize(n);

  const char * p = s.c_str();
  int i;
  for (i=0; i<n; i++)
    m_data[i] = p[i];

}

void DNAVector::Proteinize()
{
  char * p = new char[m_data.size()+1];
  int i;
  for (i=0; i<m_data.size(); i++)
    p[i] = m_data[i];
  p[m_data.size()] = 0;
  SetFromProteins(p);
  delete [] p;
}


void DNAVector::ToUpper()
{
  int i;
  for (i=0; i<m_data.size(); i++)
    m_data[i] = (char)toupper(m_data[i]);
}

void DNAVector::SetFromProteins(const string & s)
{
  int n = strlen(s.c_str());
  m_data.resize(n*3);

  char tmp[64];

  const char * p = s.c_str();

  int i;
  for (i=0; i<n; i++) {
    AminoAcidToBases(tmp, p[i]);
    m_data[3*i] = tmp[0];
    m_data[3*i+1] = tmp[1];
    m_data[3*i+2] = tmp[2];
  }

}

void DNAVector::Write(FILE * p) const
{
  int i;
  for (i=0; i<size(); i++) {
    if (i > 0 && i % 80 == 0)
      fprintf(p, "\n");
    fprintf(p, "%c", m_data[i]);
    
  }
  fprintf(p, "\n");
}


void DNAVector::Write(ostream &s) const
{
  int i;
  for (i=0; i<size(); i++) {
    if (i > 0 && i % 80 == 0)
      s << endl;
    s << m_data[i];
  }
  s << endl;
}

void DNAVector::WriteOneLine(ostream &s) const
{
  int i;
  for (i=0; i<size(); i++) {
    s << m_data[i];
  }
  s << endl;
}



void DNAVector::WriteQual(ostream &s) const
{
  int i;
  for (i=0; i<size(); i++) {
    if (i > 0 && i % 80 == 0)
      s << endl;
    s << (int) m_qual[i]<<" ";
  }
  s << endl;
}
void DNAVector::WriteQual(FILE * p) const
{
  int i;
  for (i=0; i<size(); i++) {
    if (i > 0 && i % 80 == 0)
      fprintf(p, "\n");;
    fprintf(p, "%d ", m_qual[i]);
  }
  fprintf(p, "\n");;
 
}


void DNAVector::ExtendWithString(const string& extension) {
  ExtendWithString(extension, size());
}

void DNAVector::ExtendWithString(const string& extension, int extendFrom) {
  //FILE_LOG(logDEBUG3) << "Requesting extension from: " << extendFrom << " extension size: " << extension.size();
  if(extendFrom>size() || extendFrom<0) { 
    //FILE_LOG(logWARNING) << "Requesting extension of sequence given wrong index";
    return; 
  }
  int n = extension.size();
  resize(n+extendFrom);
  for (int i=0; i<n; i++)
    (*this)[i+extendFrom] = extension[i];
}

bool DNAVector::Append(const DNAVector & d, int min, int max, double ident)
{
  int i, j;
  for (i=max; i>=min; i--) {
    int mis = (int)((1. - ident) * (double)i + 0.5);
    int no = 0;
    for (j=0; j<i; j++) {
      int x = size()-i+j;
      if (x < 0)
	return false;
      //cout << (*this)[x] << " " << d[j] << " " << x << " " << j << endl;
      if ((*this)[x] != d[j]) {
	no++;
	if (no > mis)
	  break;
      }
    }
    if (no > mis)
      continue;

    //cout << "Merge, lap=" << i << endl;
    // Merge and get out.
    for (j=0; j<i; j++) {
      int x = size()-i+j;
      if ((*this)[x] != d[j]) {
	string bases;
	bases = (*this)[x];
	bases += d[j];
	(*this)[x] = GetAmbiguous(bases);
      }
    }
    int oldSize = size();
    resize(size()-i+d.size());
    //cout << "other size=" << d.size() << " my size=" << size() << endl;
    for (j=i; j<d.size(); j++) {
      int x = oldSize-i+j;
      //cout << "Assign " << j << " to " << x << endl;
      (*this)[x] = d[j];
    }
    return true;
  }
  return false;
}

const string &DNAVector::getName() const {
	return name;
}

void DNAVector::setName(const string &newName) {
	name = newName;
}

float DNAVector::FindIdent(const DNAVector& other) const {
  int len      = min(size(), other.size());
  int matchCnt = 0;
  for (int i=0; i<len; i++){
    if((*this)[i]==other[i]) { matchCnt++; } 
  }
  return (float)matchCnt/len;
}

float DNAVector::FindIdentHP(const DNAVector& other, int max, int totaldiff) const
{
  int i = 0;
  int j = 0;
  int len      = min(size(), other.size());
  int misMatchCnt = 0;
  char a = 'N';
  char b = 'N';
  while (i<size() && j < other.size()) {    
    a = (*this)[i];
    b = other[j];
    //cout << a << " " << b << " " << i << " " << j << endl;
    if (a == b) {
      //cout << "Match." << endl;
      int ii = i;
      int jj = j;
      do {
	i++;
      } while(i<size() && (*this)[i] == a);
      do {
	j++;
      } while(j<other.size() && other[j] == b);
       
      int diff = (j-jj) - (i-ii);
      if (diff < 0)
	diff = -diff;
      if (diff >= max)
	misMatchCnt+= min(j-jj, i-ii);
    } else {
      misMatchCnt++;
      //cout << "MisMatch." << endl;
      j++;
      i++;
    }
  }
  return 1.-(float)misMatchCnt/len;
}

//////////////////////////////////////////////////////

vecDNAVector::vecDNAVector() {
	default_name_index = 0;
}

vecDNAVector::vecDNAVector(const string& file) {
	default_name_index = 0;
        Read(file);
}
//references intentionally not copied
vecDNAVector::vecDNAVector(const vecDNAVector &other) : m_data(other.m_data), default_name_index(other.default_name_index), m_name2index(other.m_name2index) {}

vecDNAVector::~vecDNAVector() {
	for(std::vector<DNAVector>::iterator currVector = m_data.begin(); currVector != m_data.end(); currVector++)
		invalidateReferences(currVector->getName());
}

//references intentionally not copied
vecDNAVector &vecDNAVector::operator = (const vecDNAVector &other) {
	if(this != &other) {
		for(std::vector<DNAVector>::iterator currVector = m_data.begin(); currVector != m_data.end(); currVector++)
			invalidateReferences(currVector->getName());
		m_data = other.m_data;
		default_name_index = other.default_name_index;
		m_name2index = other.m_name2index;
	}
	return *this;
}

const DNAVector &vecDNAVector::operator [] (int i) const {
  return m_data[i];
}

DNAVector &vecDNAVector::operator [] (int i) {
  return m_data[i];
}

const DNAVector &vecDNAVector::operator () (const string &name) const {
  return ((vecDNAVector*)this)->operator() (name);
}

DNAVector &vecDNAVector::operator () (const string &name) {
	  int index = NameIndex(name);
	  if(index == -1) {
        cout << "FATAL ERROR: could not find sequence name " << name << endl;
        exit(-1);
	  }
	  return (*this)[index];
}

bool vecDNAVector::HasChromosome(const string &name) const {
          int index = NameIndex(name);
          return (index != -1);
}

vecDNAVector::const_DNAVectorRef vecDNAVector::getDNAVectorRef(int i) const {
	return const_DNAVectorRef(this, (*this)[i].getName());
}

vecDNAVector::DNAVectorRef vecDNAVector::getDNAVectorRef(int i) {
	return DNAVectorRef(this, (*this)[i].getName());
}

vecDNAVector::const_DNAVectorRef vecDNAVector::getDNAVectorRef(const string &name) const {
	int index = NameIndex(name);
	if(index == -1)
		return const_DNAVectorRef(); //null reference
	else
		return const_DNAVectorRef(this, string(name));
}

vecDNAVector::DNAVectorRef vecDNAVector::getDNAVectorRef(const string &name) {
	int index = NameIndex(name);
	if(index == -1)
		return DNAVectorRef(); //null reference
	else
		return DNAVectorRef(this, string(name));
}

void vecDNAVector::resize(int n) {
	int prevSize = m_data.size();
  //If the size is decreased, erase the map entries for
  //indices beyond the new size and invalidate their references.
	if(n < prevSize)
		for(int i = n; i < prevSize; i++) {
			m_name2index.erase((*this)[i].getName());
			invalidateReferences((*this)[i].getName());
		}

  m_data.resize(n);
  
  	//If the size is increased, make new names for all the
	//new elements and add them to the names array and the map.
	if(n > prevSize) {
		for(int i = prevSize; i < n; i++) {
			stringstream currIndex;
			currIndex << default_name_index++;
			string currName = ">s_" + currIndex.str();
			(*this)[i].setName(currName);
			m_name2index[currName] = i;
		}
	}
}

void vecDNAVector::reserve(unsigned int newSize) {
	m_data.reserve(newSize);
}

void vecDNAVector::clear() {
	for(std::vector<DNAVector>::iterator currVector = m_data.begin(); currVector != m_data.end(); currVector++)
		invalidateReferences(currVector->getName());

  m_data.clear();
  m_name2index.clear();
}

int vecDNAVector::NameIndex(const string & name) const {
	map<string,int>::const_iterator mIter = m_name2index.find(name);
	if(mIter == m_name2index.end()) {
		string realName;
		string::size_type foundIndex = name.find('>');
		if(foundIndex == string::npos)
			realName = ">" + name;
		else if(foundIndex == 0)
			realName = name.substr(1);

		mIter = m_name2index.find(realName);
		if (mIter == m_name2index.end())
			return -1;
	}

	return mIter->second;
}

const string & vecDNAVector::Name(int i) const {
	return (*this)[i].getName();
}

const char * vecDNAVector::NameClean(int i) const {
  const char * p = (*this)[i].getName().c_str();
  return &p[1];
}

void vecDNAVector::SetName(int i, const string & s) {
	m_name2index.erase((*this)[i].getName());
	invalidateReferences((*this)[i].getName());
	m_name2index[s] = i;
	(*this)[i].setName(s);
}

void vecDNAVector::push_back(const DNAVector & v) {
  m_data.push_back(v);

  m_name2index.insert(make_pair(v.getName(), m_data.size()-1));
}

void vecDNAVector::push_back(const DNAVector & v, const string & name) {
  m_data.push_back(v);
  (*this)[m_data.size() - 1].setName(name);
  m_name2index.insert(make_pair(name, m_data.size()-1));
}

void vecDNAVector::erase(int index) {
	  m_name2index.erase((*this)[index].getName());
	  invalidateReferences(m_data[index].getName());
	  m_data.erase(m_data.begin() + index);
	  for(int i = index; i < m_data.size(); i++)
		  m_name2index[(*this)[i].getName()] = i;
}

bool vecDNAVector::erase(const string &name) {
	int index = NameIndex(name);
	if(index != -1) {
		erase(index);
		return true;
	}
	return false;
}

void vecDNAVector::fast_erase(int index) {
	  m_name2index.erase((*this)[index].getName());
	  invalidateReferences((*this)[index].getName());
	  swap((*this)[m_data.size()-1], (*this)[index]);
	  m_data.pop_back();
	  if(index < m_data.size())
		  m_name2index[(*this)[index].getName()] = index;
}

bool vecDNAVector::fast_erase(const string &name) {
	  int index = NameIndex(name);
	  if(index != -1) {
		  fast_erase(index);
		  return true;
	  }
	  return false;
}

int vecDNAVector::size() const {
	return m_data.size();
}

long long vecDNAVector::totalBases() const {
  long long sum(0);
  for (int i=0; i<(int) m_data.size(); ++i)
    sum += (*this)[i].size();

  return sum;
}


void vecDNAVector::DoProteins(bool b) 
{
  int i, j;
  for (i=0; i<m_data.size(); i++) {
    DNAVector & d = (*this)[i];
    int n = d.size() / 3;
    if (n * 3 != d.size())
      cout << "WARNING: size of sequence " << (*this)[i].getName() << " is not a multiple of 3 (protein conversion: " << d.size() << " )" << endl;

    char * p = new char[n+1];

    for (j=0; j<n; j++) {
      p[j] = trans.GetCodon(&d[j*3]);
    }
    p[j] = 0;
    //cout << "Set " << p << endl;
    if (!b) 
      d.SetFromProteins(p);
    else
      d.SetFromBases(p);
    delete [] p;
  }

}


void vecDNAVector::ReadV(const string & file) 
{
  int ver = 1;
  int i, j;
  CMReadFileStream s;
  s.Open(file.c_str());
  
  s.Read(ver);
  int new_size = 0;

  s.Read(new_size);
  resize(new_size);

  for (i=0; i<size(); i++) {
    CMString n;
    s.Read(n);
    DNAVector & d = (*this)[i];
    d.setName((const char *)n);
    m_name2index[d.getName()] = i;
    int len = 0;
    s.Read(len);
    d.resize(len);
    for (j=0; j<len; j++)
      s.Read(d[j]);
  }

  s.Close();
}

void vecDNAVector::WriteV(const string & file) const
{
  int ver = 1;
  int i, j;
  CMWriteFileStream s;
  s.Open(file.c_str());
  
  s.Write(ver);
  s.Write(size());

  for (i=0; i<size(); i++) {
    CMString n = (*this)[i].getName().c_str();
    s.Write(n);
    const DNAVector & d = (*this)[i];
    int len = d.size();
    s.Write(len);
    for (j=0; j<len; j++)
      s.Write(d[j]);
  }

  s.Close();
}



void vecDNAVector::Write(const string & fileName, bool bSkipEmpty) const
{
  FILE * p = fopen(fileName.c_str(), "w");
  int i;
  for (i=0; i<size(); i++) {
    if (bSkipEmpty && (*this)[i].size() == 0)
      continue;
    //cout << "size: " << (*this)[i].size() << " skip=" << bSkipEmpty << endl;
    if ((*this)[i].getName() == "") 
      fprintf(p, ">Sequence_%d\n", i);
    else
      fprintf(p, "%s\n", (*this)[i].getName().c_str());
    (*this)[i].Write(p);
    
  }
  fclose(p);
}

void vecDNAVector::WriteQuals(const string & fileName) const
{
  FILE * p = fopen(fileName.c_str(), "w");
  int i;
  for (i=0; i<size(); i++) {
    fprintf(p, "%s\n", (*this)[i].getName().c_str());
    (*this)[i].WriteQual(p);
  }
  fclose(p);
}



void vecDNAVector::ReadQuals(const string & fileName)
{
  if (m_data.size() == 0) {
    cout << "vecDNAVector ERROR: you need to load the bases first!" << endl;
    return;
  }
    
  FlatFileParser parser;
  
  parser.Open(fileName);
  
  int i = 0;
  int k = 0;
  int j;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    const char * p = parser.AsString(0).c_str();
    if (p[0] == '>') {

      string tmpName = parser.AsString(0);
      for (int x=1; x<parser.GetItemCount(); x++) {
	tmpName += "_";
	tmpName += parser.AsString(x);
      }

      if (tmpName != (*this)[i].getName()) {
	cout << "vecDNAVector ERROR: qual file is out of sync with fasta file!" << endl;	
	return;
      }
      k = 0;
      ++i;
      continue;
    }
    for (j=0; j<parser.GetItemCount(); j++) {
      (*this)[i-1].SetQual(k, parser.AsInt(j));
      k++;
    }
  }
}


void vecDNAVector::ReadQ(const string & fileName) // Reads a fastq file
{
  FlatFileParser parser;  
  parser.Open(fileName);
  cout << "Reading FASTQ format: " << fileName << endl;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    string n = ">" + parser.Line();
    parser.ParseLine();
    const string & s = parser.Line();
    DNAVector d;
    d.SetFromBases(s);
    // Not really efficient...
    //cout << n << " " << s << endl;
    push_back(d, n);

    parser.ParseLine();
    parser.ParseLine();
  }
  //cout << "Sequences: " << size() << endl;
  setupMap();
 
}

// Half-way efficient implementation...?
void vecDNAVector::Read(const string & fileName, bool bProteins, bool shortName, bool allUpper, bool bAppend)
{
  StringParser p;
  p.SetLine(fileName, ",");
  // p.ParseLine();
  bool bLocalAppend = bAppend; 
  for (int i=0; i<p.GetItemCount(); i++) {
    if (bLocalAppend)
      cout << "Appending " << p.AsString(i) << endl;
    ReadOne(p.AsString(i), bProteins, shortName, allUpper, bLocalAppend);
    bLocalAppend = true;
  }
}

void vecDNAVector::ReadOne(const string & fileName, bool bProteins, bool shortName, bool allUpper, bool bAppend)
{
	FlatFileParser parser;

	parser.Open(fileName);

	//	int k = 0;
	int i;

	if (!bAppend)
	  m_data.clear();

	int counter = m_data.size();
	int k = counter;
	//cout << "Counter: " << counter << endl;
	// reserve some space
	int chunk = 20000;

	// 200 Megs?
	int bigChunk = 200000000;

	//if (chunk > m_data.size())
	//m_data.resize(chunk);

	DNAVector * pVec = NULL;
	DNAVector tmpVec;

	int localcounter = 0;
	int j = 0;
	while (parser.ParseLine()) {
		if (parser.GetItemCount() == 0)
			continue;
		const char * p = parser.AsString(0).c_str();
		if (localcounter == 0 && p[0] == '@') { // It's a fastq file!!!
		  if (k == 0)
		    m_data.resize(0);
		  ReadQ(fileName);
		  return;
		}
		if (p[0] == '>') {
		  counter++;
		  localcounter++;
			if (pVec != NULL) {
				pVec->SetToSubOf(tmpVec, 0, j);
				if (allUpper)
					pVec->ToUpper();
				j = 0;
				if (bProteins)
					pVec->Proteinize();
			}


			string tmpName = parser.AsString(0);
			if ( !shortName )  {
				for (int x=1; x<parser.GetItemCount(); x++) {
					tmpName += "_";
					tmpName += parser.AsString(x);
				}
			}
			//cout << "In: " << tmpName << endl;

			if (k >= m_data.size())
				m_data.resize(k + chunk);

			pVec = &(*this)[k];
			pVec->setName(tmpName);
			k++;
			continue;
		}
		int n = strlen(p);

		//cout << "Parsing: " << p << endl;
		for (i=0; i<n; i++) {
			if (j >= tmpVec.size())
				tmpVec.resize(j + bigChunk);
			tmpVec[j] = p[i];
			j++;
		}

	}

	if (pVec != NULL) {
		pVec->SetToSubOf(tmpVec, 0, j);
		if (allUpper)
			pVec->ToUpper();
		if (bProteins)
			pVec->Proteinize();
	}

	//cout << "Read sequences: " << k << endl;
	m_data.resize(k);
	setupMap();
}

void vecDNAVector::Read(const string & fileName, std::vector<string> & names)
{
	Read(fileName);
	int i;
	names.resize(m_data.size());
	for (i=0; i<m_data.size(); i++) {
		const char * p = (*this)[i].getName().c_str();
		//if (strlen(p) > 0)
		names[i] = &p[1];
	}
}

void vecDNAVector::ReadPartial(const string & fileName, unsigned int firstToRead, unsigned int numToRead, bool bProteins, bool shortName, bool allUpper)
{
	if(numToRead == 0)
		return;

	FlatFileParser parser;
	parser.Open(fileName);

	clear();

	// reserve some space
	const int chunk = 20000;

	// 200 Megs?
	const int bigChunk = 200000000;

	m_data.resize(chunk);

	DNAVector tmpVec;

	unsigned int currSequenceLength = 0;
	int currIndex = -1;
	unsigned int numRead = 0;

	string currName;
	//When this terminates currName will hold the name of the first sequence to read
	while(currIndex < (int)firstToRead && parser.ParseLine()) {
		if (parser.AsString(0)[0] == '>') {
			currIndex++;
			currName = parser.AsString(0);
			if ( !shortName )  {
				for (int x=1; x<parser.GetItemCount(); x++) {
					currName += "_";
					currName += parser.AsString(x);
				}
			}
		}
	}

	while (numRead < numToRead && parser.ParseLine()) {
		if (parser.GetItemCount() == 0)
			continue;
		const char * p = parser.AsString(0).c_str();
		if (p[0] == '>') {
			(*this)[numRead].setName(currName);
			(*this)[numRead].SetToSubOf(tmpVec, 0, currSequenceLength);
			currSequenceLength = 0;

			if (allUpper)
				(*this)[numRead].ToUpper();
			if (bProteins)
				(*this)[numRead].Proteinize();
			numRead++;

			currName = parser.AsString(0);
			if ( !shortName )  {
				for (int x=1; x<parser.GetItemCount(); x++) {
					currName += "_";
					currName += parser.AsString(x);
				}
			}


			if (numRead >= m_data.size())
				m_data.resize(numRead + chunk);
		}
		else {
			int n = strlen(p);

			for (int i=0; i<n; i++) {
				if (currSequenceLength >= (unsigned int)tmpVec.size())
					tmpVec.resize(currSequenceLength + bigChunk);
				tmpVec[currSequenceLength] = p[i];
				currSequenceLength++;
			}
		}
	}

	if (numRead < numToRead) {
		(*this)[numRead].setName(currName);
		(*this)[numRead].SetToSubOf(tmpVec, 0, currSequenceLength);
		currSequenceLength = 0;

		if (allUpper)
			(*this)[numRead].ToUpper();
		if (bProteins)
			(*this)[numRead].Proteinize();

		numRead++;
	}

	m_data.resize(numRead);
	setupMap();
}

void vecDNAVector::ReadMultiple(const vector<string> &fileNames, bool bProteins, bool shortName, bool allUpper) {
	m_data.clear();

	int k = 0;
	int i;
	int j = 0;
	DNAVector * pVec = NULL;
	DNAVector tmpVec;

	// reserve some space
	int chunk = 20000;

	// 200 Megs?
	int bigChunk = 200000000;

	m_data.resize(chunk);

	for(unsigned int currFile = 0; currFile < fileNames.size(); currFile++) {
		FlatFileParser parser;
		parser.Open(fileNames[currFile]);

		while (parser.ParseLine()) {
			if (parser.GetItemCount() == 0)
				continue;
			const char * p = parser.AsString(0).c_str();
			if (p[0] == '>') {
				if (pVec != NULL) {
					pVec->SetToSubOf(tmpVec, 0, j);
					if (allUpper)
						pVec->ToUpper();
					j = 0;
					if (bProteins)
						pVec->Proteinize();
				}


				string tmpName = parser.AsString(0);
				if ( !shortName )  {
					for (int x=1; x<parser.GetItemCount(); x++) {
						tmpName += "_";
						tmpName += parser.AsString(x);
					}
				}
				//cout << "In: " << tmpName << endl;

				if (k >= m_data.size())
					m_data.resize(k + chunk);

				pVec = &(*this)[k];
				pVec->setName(tmpName);
				k++;
				continue;
			}
			int n = strlen(p);

			//cout << "Parsing: " << p << endl;
			for (i=0; i<n; i++) {
				if (j >= tmpVec.size())
					tmpVec.resize(j + bigChunk);
				tmpVec[j] = p[i];
				j++;
			}

		}

		if (pVec != NULL) {
			pVec->SetToSubOf(tmpVec, 0, j);
			if (allUpper)
				pVec->ToUpper();
			j = 0;
			if (bProteins)
				pVec->Proteinize();

			pVec = NULL;
		}
	}

	//cout << "Read sequences: " << k << endl;
	m_data.resize(k);
	setupMap();
}

void vecDNAVector::ReverseComplement() {
  for (int i=0; i < m_data.size(); i++)
    (*this)[i].ReverseComplement();
}

void vecDNAVector::UniqueSort() {
	map<DNAVector, string> tempNameMap;
	for(int i = 0; i < m_data.size(); i++)
		tempNameMap[(*this)[i]] = (*this)[i].getName();

  m_name2index.clear();
  //::UniqueSort(m_data);
  std::sort(m_data.begin(), m_data.end()); 
  auto last = std::unique(m_data.begin(), m_data.end());
  m_data.erase(last, m_data.end());

  for(int i = 0; i < m_data.size(); i++) {
  	m_name2index[tempNameMap[(*this)[i]]] = i;
  	tempNameMap.erase((*this)[i]);
  }
  for(map<DNAVector, string>::iterator currRemoved = tempNameMap.begin(); currRemoved != tempNameMap.end(); currRemoved++)
  	invalidateReferences(currRemoved->second);
}

void vecDNAVector::Sort() {
  map<DNAVector, string> tempNameMap;
  for(int i = 0; i < m_data.size(); i++)
    tempNameMap[(*this)[i]] = (*this)[i].getName();

  m_name2index.clear();
  std::sort(m_data.begin(),m_data.end());
  for(int i = 0; i < m_data.size(); i++) {
    m_name2index[tempNameMap[(*this)[i]]] = i;
    tempNameMap.erase((*this)[i]);
  }
  //for(map<DNAVector, string>::iterator currRemoved = tempNameMap.begin(); currRemoved != tempNameMap.end(); currRemoved++)
  //invalidateReferences(currRemoved->second);
}

bool vecDNAVector::SetSequence(const Coordinate& coords, DNAVector& resultSeq) const {
  if(!HasChromosome(coords.getChr())) { 
    //FILE_LOG(logWARNING) << "Check Genome data! - Chromosome: " << coords.getChr() 
    //<< " was not found in the given fasta file";
    return false; 
  }   
  bool set = resultSeq.SetToSubOf(this->operator()(coords.getChr()), coords.getStart(), coords.getStop()-coords.getStart()+1);
  if (coords.isReversed())
    resultSeq.ReverseComplement();
  return set;
}

void vecDNAVector::addReference(ReferenceListener *newReference) const {
	  references[newReference->name()].push_back(newReference);
}

void vecDNAVector::removeReference(ReferenceListener *toRemove) const {
	map<string, vector<ReferenceListener *> >::iterator containingVector = references.find(toRemove->name());
	if(containingVector != references.end()) {
		for(vector<ReferenceListener *>::iterator currReference = containingVector->second.begin(); currReference != containingVector->second.end(); currReference++) {
			if(*currReference == toRemove) {
				containingVector->second.erase(currReference);
				if(containingVector->second.empty())
					references.erase(toRemove->name());
				break;
			}
		}
	}
}

void vecDNAVector::invalidateReferences(string toInvalidate) {
	  string realName = toInvalidate;
	map<string, vector<ReferenceListener *> >::iterator containingVector = references.find(realName);
	if(containingVector == references.end()) {
		  string::size_type foundIndex = toInvalidate.find('>');
		  if(foundIndex == string::npos)
			  realName = ">" + toInvalidate;
		  else if(foundIndex == 0)
			  realName = toInvalidate.substr(1);

		  containingVector = references.find(realName);
	}

	if(containingVector != references.end())
		for(vector<ReferenceListener *>::iterator currReference = containingVector->second.begin(); currReference != containingVector->second.end(); currReference++)
			(*currReference)->invalidate();

	references.erase(realName);
}

void vecDNAVector::setupMap() {
	if ( m_data.empty() )
		return;

	for ( int i=0; i<(int) m_data.size(); ++i )
		m_name2index.insert(make_pair((*this)[i].getName(),i));
}

vecDNAVector::ReferenceListener::~ReferenceListener() {}

vecDNAVector::DNAVectorRef::DNAVectorRef(vecDNAVector * owner, string DNAName) : owner(owner), DNAName(DNAName) {
	owner->addReference(this);
}

vecDNAVector::DNAVectorRef::DNAVectorRef(const DNAVectorRef &toCopy) : ReferenceListener() {
	  owner = toCopy.owner;
	  DNAName = toCopy.DNAName;
	  if(owner != NULL)
		  owner->addReference(this);
}

vecDNAVector::DNAVectorRef::DNAVectorRef() {
	  owner = NULL;
	  DNAName = "";
}

vecDNAVector::DNAVectorRef::~DNAVectorRef() {
	if(owner != NULL)
		owner->removeReference(this);
}

vecDNAVector::DNAVectorRef &vecDNAVector::DNAVectorRef::operator = (const DNAVectorRef &other) {
	if(this != &other) {
		if(owner != NULL)
			owner->removeReference(this);

		owner = other.owner;
		DNAName = other.DNAName;
		if(owner != NULL)
			owner->addReference(this);
	}
	return *this;
}

bool vecDNAVector::DNAVectorRef::operator == (const DNAVectorRef & other) const {
	string strippedName = DNAName;
	if(strippedName.size() > 0 && strippedName[0] == '>')
		strippedName = strippedName.substr(1);

	string otherStrippedName = other.name();
	if(otherStrippedName.size() > 0 && otherStrippedName[0] == '>')
		otherStrippedName = otherStrippedName.substr(1);

	return owner == other.owner && strippedName == otherStrippedName;
}

bool vecDNAVector::DNAVectorRef::operator == (const const_DNAVectorRef & other) const {
	string strippedName = DNAName;
	if(strippedName.size() > 0 && strippedName[0] == '>')
		strippedName = strippedName.substr(1);

	string otherStrippedName = other.name();
	if(otherStrippedName.size() > 0 && otherStrippedName[0] == '>')
		otherStrippedName = otherStrippedName.substr(1);

	return owner == other.owner && strippedName == otherStrippedName;
}

bool vecDNAVector::DNAVectorRef::operator != (const DNAVectorRef & other) const {
	return !(*this == other);
}

bool vecDNAVector::DNAVectorRef::operator != (const const_DNAVectorRef & other) const {
	return !(*this == other);
}

DNAVector &vecDNAVector::DNAVectorRef::operator *() {
	  if(owner == NULL) {
		  cout << "FATAL ERROR: Tried to dereference an invalid DNAVectorRef for DNAVector " << DNAName << "." << endl;
		  exit(-1);
	  }
	  else
		  return (*owner)(DNAName);
}

DNAVector *vecDNAVector::DNAVectorRef::operator ->() {
	  if(owner == NULL) {
		  cout << "FATAL ERROR: Tried to dereference an invalid DNAVectorRef for DNAVector " << DNAName << "." << endl;
		  exit(-1);
	  }
	  else
		  return &(*owner)(DNAName);
}

string vecDNAVector::DNAVectorRef::name() const {
	return DNAName;
}

bool vecDNAVector::DNAVectorRef::isInvalid() const {
	return owner == NULL;
}

void vecDNAVector::DNAVectorRef::invalidate() {
	owner = NULL;
}

vecDNAVector::const_DNAVectorRef::const_DNAVectorRef(const vecDNAVector * owner, string DNAName) : owner(owner), DNAName(DNAName) {
	owner->addReference(this);
}

vecDNAVector::const_DNAVectorRef::const_DNAVectorRef(const const_DNAVectorRef &toCopy) : ReferenceListener() {
	owner = toCopy.owner;
	DNAName = toCopy.DNAName;
	if(owner != NULL)
		owner->addReference(this);
}

vecDNAVector::const_DNAVectorRef::const_DNAVectorRef(const DNAVectorRef &toCopy) : ReferenceListener() {
	owner = toCopy.owner;
	DNAName = toCopy.DNAName;
	if(owner != NULL)
		owner->addReference(this);
}

vecDNAVector::const_DNAVectorRef::const_DNAVectorRef() {
	  owner = NULL;
	  DNAName = "";
}

vecDNAVector::const_DNAVectorRef::~const_DNAVectorRef() {
	if(owner != NULL)
		owner->removeReference(this);
}

vecDNAVector::const_DNAVectorRef &vecDNAVector::const_DNAVectorRef::operator = (const const_DNAVectorRef &other) {
	if(this != &other) {
		if(owner != NULL)
			owner->removeReference(this);

		owner = other.owner;
		DNAName = other.DNAName;
		if(owner != NULL)
			owner->addReference(this);
	}
	return *this;
}

vecDNAVector::const_DNAVectorRef &vecDNAVector::const_DNAVectorRef::operator = (const DNAVectorRef &other) {
	if(owner != NULL)
		owner->removeReference(this);

	owner = other.owner;
	DNAName = other.DNAName;
	if(owner != NULL)
		owner->addReference(this);
	return *this;
}

bool vecDNAVector::const_DNAVectorRef::operator == (const DNAVectorRef & other) const {
	string strippedName = DNAName;
	if(strippedName.size() > 0 && strippedName[0] == '>')
		strippedName = strippedName.substr(1);

	string otherStrippedName = other.name();
	if(otherStrippedName.size() > 0 && otherStrippedName[0] == '>')
		otherStrippedName = otherStrippedName.substr(1);

	return owner == other.owner && strippedName == otherStrippedName;
}

bool vecDNAVector::const_DNAVectorRef::operator == (const const_DNAVectorRef & other) const {
	string strippedName = DNAName;
	if(strippedName.size() > 0 && strippedName[0] == '>')
		strippedName = strippedName.substr(1);

	string otherStrippedName = other.name();
	if(otherStrippedName.size() > 0 && otherStrippedName[0] == '>')
		otherStrippedName = otherStrippedName.substr(1);

	return owner == other.owner && strippedName == otherStrippedName;
}

bool vecDNAVector::const_DNAVectorRef::operator != (const DNAVectorRef & other) const {
	return !(*this == other);
}

bool vecDNAVector::const_DNAVectorRef::operator != (const const_DNAVectorRef & other) const {
	return !(*this == other);
}

const DNAVector &vecDNAVector::const_DNAVectorRef::operator *() const {
	  if(owner == NULL) {
		  cout << "FATAL ERROR: Tried to dereference an invalid const_DNAVectorRef for DNAVector " << DNAName << "." << endl;
		  exit(-1);
	  }
	  else
		  return (*owner)(DNAName);
}

const DNAVector *vecDNAVector::const_DNAVectorRef::operator ->() const {
	  if(owner == NULL) {
		  cout << "FATAL ERROR: Tried to dereference an invalid const_DNAVectorRef for DNAVector " << DNAName << "." << endl;
		  exit(-1);
	  }

	  return &(*owner)(DNAName);
}

string vecDNAVector::const_DNAVectorRef::name() const {
	return DNAName;
}

bool vecDNAVector::const_DNAVectorRef::isInvalid() const {
	return owner == NULL;
}

void vecDNAVector::const_DNAVectorRef::invalidate() {
	owner = NULL;
}


int CountNs(const DNAVector & d)
{
  int n = 0;
  for (int i=0; i<d.size(); i++) {
    if (d[i] == 'N' || d[i] == 'X')
      n++;
  }
  return n;
}

bool IsHomopolymer(const DNAVector &d, double threshold)
{
  map<char,int> basecounts;
  for (int i=0; i<d.size(); ++i)
  {
    basecounts[d[i]]++;
  }

  map<char,int>::iterator mIter = basecounts.begin();
  unsigned int max(basecounts.count(mIter->first));
  for (; mIter != basecounts.end(); ++mIter)
    if (basecounts.count(mIter->first) > max)
      max = basecounts.count(mIter->first);
  
  return ((double) max/d.size() >= threshold);
}




//=================================================================

void vecDNAVectorStream::ReadStream(const string & fileName)
{
  if (m_pParser != NULL)
    delete m_pParser;

  m_pParser = new FlatFileParser;
  
  m_pParser->Open(fileName);
}


const DNAVector & vecDNAVectorStream::Next()
{

  m_seq.resize(0);
  while (m_pParser->ParseLine()) {
    if (m_pParser->GetItemCount() == 0)
      continue;
    const char * p = m_pParser->AsString(0).c_str();
    if (p[0] == '>') {
      break;
    }
    DNAVector tmp;
    tmp.SetFromBases(m_pParser->Line());
    m_seq += tmp;
  }
  return m_seq;
}

vecDNAVectorStream::~vecDNAVectorStream() 
{
  if (m_pParser != NULL)
    delete m_pParser;
}
