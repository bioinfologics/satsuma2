//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif

#include "aligns/KmerAlignCore.h"
#include <string>

#include "base/CommandLineParser.h"
#include "analysis/CrossCorr.h"
#include "analysis/XCorrDynProg.h"
#include "util/SysTime.h"
#include "analysis/MakeSixFrames.h"
#include "analysis/MultiProtein.h"



class ProteinXCorr
{
public:
  ProteinXCorr() {
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



int main( int argc, char** argv )
{


  commandArg<string> aStringCmmd("-i","fasta file");
 

  commandLineParser P(argc,argv);
  P.SetDescription("Satsuma-based (cross-correlation) protein alignment tool.");
  P.registerArg(aStringCmmd);


  P.parse();

  string REF = P.GetStringValueFor(aStringCmmd);


  vecDNAVector dna;
 
  dna.Read(REF);

  int i, j;

  ProteinXCorr xc;
 

  int default_xc_size = 4096;


  for (i=0; i<dna.isize(); i++) {
    const DNAVector & t = dna[i];
    DNAVector d = t;
    d.ReverseComplement();

    int size = xc.Size(t.isize(), d.isize());
    svec<float> signal;

    
    cout << "Set up RNA signals" << endl;
    CCSignal target, query;      
    target.SetSequence(t, size);
    query.SetSequence(d, size);


    //const svec<float> & a = target.Get(0);
    //for (j=0; j<a.isize(); j++) {
    //  cout << "A: " << a[j] << " " << t[j] << endl;
    // }



    cout << "Start XC from here." << endl;
    xc.CrossCorrelate(signal, target, query);
    cout << "Signal size (RNA) =" << signal.isize() << " q=" << query.GetFullSize() << " t=" << target.GetFullSize() << endl;
      
    SignalFilter filter;

    //for (j=0; j<signal.isize(); j++)
    //cout << "xcsig " << signal[i] << endl;


 
    filter.Do(signal);
      
    XCDynProg dp;
    MultiRNA tt(t);
    MultiRNA dd(d);


    dp.SetExpect(0.35);
    double score = dp.Align(tt, dd, filter);
      
    if (score < 0.5)
      cout << "Summary " << dna.NameClean(i) << " score: " << score << endl << endl << endl;

    
  }


  return 0;
}
