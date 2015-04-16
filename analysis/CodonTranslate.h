
#ifndef CODONTRANSLATE_H
#define CODONTRANSLATE_H

#include "base/SVector.h"
#include <string>

inline char as_i(char b) {
  char i = -1;

  switch(b) {
  case 'A':
    i = 0;
    break;
  case 'C':
    i = 1;
    break;
  case 'G':
    i = 2;
    break;
  case 'T':
    i = 3;
    break;
  default:   
    break;
  }
  return i;
}


class CodonMatrix
{
 public:
  CodonMatrix() {
    m_row.resize(21, 0);
  }
  void Set(int i, double v) {
    m_row[i] = v;
  }

  double Get(int i) const {return m_row[i];}
  int Size() const {return m_row.isize();}
 private:
  svec<double> m_row;
};


class CodonTranslator
{
 public:
  CodonTranslator() {
    m_index.resize(64, 0);
    m_codons.resize(21, 0);
    
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

    m_matrix.resize(21);
    
    SetFullMatrix();

  }

  int GetCodonIndexByIndex(char one, char two, char three) {
    int i = GetIndex(one, two, three);
    return m_index[i];
  }
  
  char GetCodonByIndex(char one, char two, char three) {
    int i = GetIndex(one, two, three);
    return m_codons[m_index[i]];
  }

  int GetCodonIndex(char one, char two, char three) {
    //cout << "one=" << one << " idx=" << (int)as_i(one) << endl;
    int i = GetIndex(as_i(one), as_i(two), as_i(three));
    //cout << "Index int=" << i << " m_index=" << m_index[i] << endl;
    return m_index[i];
  }
  
  char GetCodon(char one, char two, char three) {
    int i = GetIndex(as_i(one), as_i(two), as_i(three));
    return m_codons[m_index[i]];
  }

  const CodonMatrix & GetMatrix(int idx) const {return m_matrix[idx];}

 private:
  svec<char> m_codons;
  svec<int> m_index;
  svec<CodonMatrix> m_matrix;

  int GetIndex(char one, char two, char three) {
    
    int i = one;
    i *= 4;
    i += (int)two;
    i *= 4;
    i += (int)three;
    return i;
  }

  void Set(const char * codon, const char * three, int index) {
    
    int i = GetIndex(as_i(three[0]), as_i(three[1]), as_i(three[2]));
    m_index[i] = index;
    //cout << "i=" << i << " index=" << index << endl;
    m_codons[index] = codon[0];
  }

  void SetFullMatrix();

  void SetMatrix(char a, char b, double v) {
    int i;
    int idx1 = Find(a);
    int idx2 = Find(b);
    if (idx1 == -1 || idx2 == -1)
      return;

    CodonMatrix & m = m_matrix[idx1];
    m.Set(idx2, v);
  }

  int Find(char a) {
    for (int i=0; i<m_codons.isize(); i++) {
      if (m_codons[i] == a) {
	return i;
      }
    }
    return -1;
  }
};





#endif //CODONTRANSLATE_H

