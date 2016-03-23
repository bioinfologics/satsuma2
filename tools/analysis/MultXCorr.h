#ifndef MULTXCORR_H
#define MULTXCORR_H

#include "analysis/CrossCorr.h"

class MultiSizeXCorr
{
public:
  MultiSizeXCorr() {
    m_xc.resize(32);
  }

  int Size(int sizeA, int sizeB) {
    //cout << "Sizes: " << sizeA << " " << sizeB << endl;
    
    if (sizeB > sizeA)
      sizeA = sizeB;

    int b = 1;
    while (b < sizeA) {
      b *= 2;
    }
    // Take care of 0 padding to avoid wrap-around.
    b *= 2;

    if (b < 512)
      b = 512;
    //cout << "Full size=" << b << endl;
    return b;     
  }

  
  void CrossCorrelate(svec<float> & out, const CCSignal & one, const CCSignal & two) {
    int i = Index(one.GetFullSize());
    if (i >= m_xc.isize()) {
      cout << "ERROR! FFT size not implemented: " << one.GetFullSize() << endl;
      throw;
      return;
    }

    //cout << "Start XC, i=" << i << endl;
    m_xc[i].CrossCorrelate(out, one, two);
    //cout << "done." << endl;
  }



private:
  int Index(int size) {
    int i = 0;
    while (size > 2) {
      i++;
      size /= 2;
    }
    return i;
  }

  svec<CrossCorrelation> m_xc;

};



#endif 

