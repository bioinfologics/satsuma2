#ifndef DNAMATRIX_H
#define DNAMATRIX_H

#include "DNAVector.h"



char DNA_Letter(double a, double c, double g, double t) 
{
  if (a > c && a > g && a > t)
    return 'A';
  if (c > g && c > t)
    return 'C';
  if (g > t)
    return 'G';
 
  return 'T';
  //return 'N';
}

class DNAMatrix
{
 public:
  DNAMatrix() {}
  
  void SetDNA(const DNAVector & d) {
    m_data.resize(d.isize()*4, 0.);

    int i;
    for (i=0; i<d.isize(); i++) {
      m_data[4*i]   = DNA_A(d[i]);
      m_data[4*i+1] = DNA_C(d[i]);
      m_data[4*i+2] = DNA_G(d[i]);
      m_data[4*i+3] = DNA_T(d[i]);
    }
  }
  
  void SetData(const svec<double> & d) {
    m_vec.resize(d.isize()/4);
    m_data = d;
    int i;
    for (i=0; i<d.isize()/4; i++) {
      char l = DNA_Letter(m_data[4*i], 
			  m_data[4*i+1],
			  m_data[4*i+2],
			  m_data[4*i+3]);
      m_vec[i] = l;
    }
  }


  void Print() const {
    int i;
    for (i=0; i<m_vec.isize(); i++) {
      cout << m_vec[i];
    }
    cout << endl;
    for (i=0; i<m_vec.isize(); i++) {
      cout << m_vec[i] << "\t" << m_data[4*i] << "\t";
      cout << m_data[4*i+1] << "\t" << m_data[4*i+2] << "\t" << m_data[4*i+3] << endl;
    }
  }

  svec<double> & DataDirect() {return m_data;}
  const svec<double> & Data() const {return m_data;}
  const DNAVector & DNA() const {return m_vec;}

 private:
  DNAVector m_vec;
  svec<double> m_data;
};


#endif 

