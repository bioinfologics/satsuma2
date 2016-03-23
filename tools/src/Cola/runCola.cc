#include <string>

#include "base/CommandLineParser.h"
#include "src/Cola/Cola.h"


int main(int argc,char** argv)
{
  commandArg<string> aStringCmd("-q","query sequence");
  commandArg<string> bStringCmd("-t","target sequence");
  commandArg<int>    alignerTypeCmd("-a","Aligner type - Choose 1 : NSGA , 2 : NS , 3 : SWGA, 4 : SW ");
  commandArg<int>    gapOpenCmd("-o", "Aligner gap open penalty", false);
  commandArg<int>    mismatchCmd("-m", "Aligner Mismatch penalty", false);
  commandArg<int>    gapExtCmd("-e", "Aligner gap extension penalty", false);
  commandArg<bool>   selfCmd("-all","all alignments", false);
  commandArg<double> maxPCmd("-p","Maximum acceptable P-value", 1.0);
  commandArg<double> minIdentCmd("-i","Minium acceptable identity", 0.0);
  commandArg<int>    bandedCmd("-b", "The bandwidth for banded mode, default is for unbanded", -1);

  commandLineParser P(argc,argv);
  P.SetDescription("Aligns two fasta files using Cola");
  P.registerArg(aStringCmd);
  P.registerArg(bStringCmd);
  P.registerArg(alignerTypeCmd);
  P.registerArg(gapOpenCmd);
  P.registerArg(mismatchCmd);
  P.registerArg(gapExtCmd);
  P.registerArg(selfCmd);
  P.registerArg(maxPCmd);
  P.registerArg(minIdentCmd);
  P.registerArg(bandedCmd);

  P.parse();

  string      aString     = P.GetStringValueFor(aStringCmd);
  string      bString     = P.GetStringValueFor(bStringCmd);
  AlignerType aType       = AlignerType(P.GetIntValueFor(alignerTypeCmd));
  int         gapOpenPen  = P.GetIntValueFor(gapOpenCmd);
  int         mismatchPen = P.GetIntValueFor(mismatchCmd);
  int         gapExtPen   = P.GetIntValueFor(gapExtCmd);
  bool        bAll        = P.GetBoolValueFor(selfCmd);
  double      maxP        = P.GetDoubleValueFor(maxPCmd);
  double      minIdent    = P.GetDoubleValueFor(minIdentCmd);
  int         banded      = P.GetIntValueFor(bandedCmd);

  vecDNAVector query, target;

  query.Read(aString);
  target.Read(bString);
  query.RemoveGaps();
  target.RemoveGaps();
  int i, j;
  
  for (i=0; i<target.isize(); i++) {
    for (j=0; j<query.isize(); j++) {
      if (bAll || i==j) {
        Cola cola1 = Cola();
        if(gapOpenPen && mismatchPen && gapExtPen) {
          cola1.createAlignment(target[i], query[j],
             AlignerParams(banded, aType, -gapOpenPen, -mismatchPen, -gapExtPen));
        } else { // If params are not given, use default mode
          cola1.createAlignment(target[i], query[j], AlignerParams(banded, aType));
        }
        Alignment cAlgn = cola1.getAlignment();
        if(cAlgn.calcPVal()<=maxP && cAlgn.getIdentityScore()>=minIdent) {
          cout << target.Name(i) << " vs " << query.Name(j) << endl;
          cAlgn.print(0,1,cout,100);
        } else {
          cout<<"No Alignment at given significance threshold"<<endl;
        }
      }
    }
  }

  return 0;

}

