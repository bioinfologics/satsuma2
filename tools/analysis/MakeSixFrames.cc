
#include "analysis/MakeSixFrames.h"


int Mult3(int l) {
  return 3 * (l/3);
}

void MakeSixFrames::AllSixFrames(vecDNAVector & out, const vecDNAVector & in)
{
  out.resize(in.isize()*6);

  int i;
  int k = 0;
  for (i=0; i<in.isize(); i++) {
    //cout << i << endl;

    DNAVector d = in[i];
    int n = d.isize();
    

    out[k].SetToSubOf(d, 0, Mult3(n));
    out.SetName(k, in.Name(i) + " frame:0+");
    k++;
    out[k].SetToSubOf(d, 1, Mult3(n-1));
    out.SetName(k, in.Name(i) + " frame:1+");
    k++;
    out[k].SetToSubOf(d, 2, Mult3(n-2));
    out.SetName(k, in.Name(i) + " frame:2+");
    k++;

    d.ReverseComplement();
  
    out[k].SetToSubOf(d, 0, Mult3(n));
    out.SetName(k, in.Name(i) + " frame:0-");
    k++;
    out[k].SetToSubOf(d, 1, Mult3(n-1));
    out.SetName(k, in.Name(i) + " frame:1-");
    k++;
    out[k].SetToSubOf(d, 2, Mult3(n-2));
    out.SetName(k, in.Name(i) + " frame:2-");
    k++;
  }
}


void MakeSixFrames::AllSixFrames(vecDNAVector & inout)
{
  vecDNAVector tmp;
  AllSixFrames(tmp, inout);
  //cout << "Copy." << endl;
  inout = tmp;
  //cout << "Done." << endl;
}
