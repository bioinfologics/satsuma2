#ifndef MULTIPROTEIN
#define MULTIPROTEIN

#include "analysis/DNAVector.h"

class CMReadFileStream;
class CMWriteFileStream;

typedef std::vector<double> SeqFeature;



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
    if (m_feat[0].size() > 0 && d.size() != m_feat[0].size()) {
      cout << "ERROR: MultiProtein sizes don't match!!" << endl;
      throw;
    }
    if (m_feat[0].size() == 0) {
      for (i=0; i<m_feat.size(); i++)
	m_feat[i].resize(d.size(), 0.);    
      m_seq = d;
    }

    for (j=0; j<21; j++) {
      SeqFeature & s = m_feat[j];
      char letter = AminoAcidLetter(j);
      for (i=0; i<s.size(); i++) {
	if (d[i] == letter)
	  s[i] += 1.;
      }

    }
    for (i=0; i<d.size(); i++) {
      if (m_seq[i] != d[i])
	m_seq[i] = '?';
    }

  }

  void Read(CMReadFileStream & f);
  void Write(CMWriteFileStream & f);


  const char & operator [] (int i) const {return m_seq[i];}
  int size() const {return m_feat[0].size();}

  double Distance(const MultiProtein & seq, int x, int y) const {
    int i;
    //cout << "x=" << x << " y=" << y << endl;
    for (i=0; i<m_feat.size(); i++) {
      const SeqFeature & f = m_feat[i];
      const SeqFeature & e = seq.m_feat[i];
      if (f[x] > 0.0 && e[y] > 0.)
	return 0.;
    }
    return 1.;
  }

 
  
 protected:
  DNAVector m_seq;
  std::vector<SeqFeature> m_feat;
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
    if (m_feat[0].size() > 0 && d.size() != m_feat[0].size()) {
      cout << "ERROR: MultiRNA sizes don't match!!" << endl;
      throw;
    }
    if (m_feat[0].size() == 0) {
      for (i=0; i<m_feat.size(); i++)
	m_feat[i].resize(d.size(), 0.);    
      m_seq = d;
    }

    for (j=0; j<m_feat.size(); j++) {
      SeqFeature & s = m_feat[j];
      char letter = NucLetter(j);
      for (i=0; i<s.size(); i++) {
	if (d[i] == letter)
	  s[i] += 1.;
      }

    }
    for (i=0; i<d.size(); i++) {
      if (m_seq[i] != d[i])
	m_seq[i] = '?';
    }

  }


};




#endif //SEQWITHFEATURES
