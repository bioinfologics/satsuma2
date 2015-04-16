#ifndef MULTIPROTEIN
#define MULTIPROTEIN

#include "analysis/DNAVector.h"

class CMReadFileStream;
class CMWriteFileStream;

typedef svec<double> SeqFeature;



class MultiProtein
{
 public:

  MultiProtein() {
    m_feat.resize(21);
  }
  MultiProtein(const DNAVector & d) {
    m_feat.resize(21);
    AddSequence(d);
  }

  const DNAVector & Sequence() const {return m_seq;}
  
  virtual void AddSequence(const DNAVector & d) {
    int i, j;
    if (m_feat[0].isize() > 0 && d.isize() != m_feat[0].isize()) {
      cout << "ERROR: MultiProtein sizes don't match!!" << endl;
      throw;
    }
    if (m_feat[0].isize() == 0) {
      for (i=0; i<m_feat.isize(); i++)
	m_feat[i].resize(d.isize(), 0.);    
      m_seq = d;
    }

    for (j=0; j<21; j++) {
      SeqFeature & s = m_feat[j];
      char letter = AminoAcidLetter(j);
      for (i=0; i<s.isize(); i++) {
	if (d[i] == letter)
	  s[i] += 1.;
      }

    }
    for (i=0; i<d.isize(); i++) {
      if (m_seq[i] != d[i])
	m_seq[i] = '?';
    }

  }

  void Read(CMReadFileStream & f);
  void Write(CMWriteFileStream & f);


  const char & operator [] (int i) const {return m_seq[i];}
  int size() const {return m_feat[0].isize();}
  int isize() const {return m_feat[0].isize();}

  double Distance(const MultiProtein & seq, int x, int y) const {
    int i;
    //cout << "x=" << x << " y=" << y << endl;
    for (i=0; i<m_feat.isize(); i++) {
      const SeqFeature & f = m_feat[i];
      const SeqFeature & e = seq.m_feat[i];
      if (f[x] > 0.0 && e[y] > 0.)
	return 0.;
    }
    return 1.;
  }

 
  
 protected:
  DNAVector m_seq;
  svec<SeqFeature> m_feat;
};



class MultiRNA : public MultiProtein
{
 public:
  MultiRNA() {
    m_feat.resize(4);
  }
  MultiRNA(const DNAVector & d) {
    m_feat.resize(4);
    AddSequence(d);
  }
 
  
  virtual void AddSequence(const DNAVector & d) {
    int i, j;
    if (m_feat[0].isize() > 0 && d.isize() != m_feat[0].isize()) {
      cout << "ERROR: MultiRNA sizes don't match!!" << endl;
      throw;
    }
    if (m_feat[0].isize() == 0) {
      for (i=0; i<m_feat.isize(); i++)
	m_feat[i].resize(d.isize(), 0.);    
      m_seq = d;
    }

    for (j=0; j<m_feat.isize(); j++) {
      SeqFeature & s = m_feat[j];
      char letter = NucLetter(j);
      for (i=0; i<s.isize(); i++) {
	if (d[i] == letter)
	  s[i] += 1.;
      }

    }
    for (i=0; i<d.isize(); i++) {
      if (m_seq[i] != d[i])
	m_seq[i] = '?';
    }

  }


};




#endif //SEQWITHFEATURES
