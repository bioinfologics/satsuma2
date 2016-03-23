#include <string>
#include <ctime>

#include "extern/easel/esl_gumbel.h"

#include "base/CommandLineParser.h"
#include "src/Cola/Cola.h"

int main(int argc,char** argv)
{
  commandArg<int>   alignerTypeCmd("-a","Aligner type - Choose 1 : NSGA , 2 : NS , 3 : SWGA, 4 : SW ");
  commandArg<int>   gapOpenCmd("-o", "Aligner gap open penalty");
  commandArg<int>   mismatchCmd("-m", "Aligner Mismatch penalty");
  commandArg<int>   gapExtCmd("-e", "Aligner gap extension penalty");

  commandLineParser P(argc,argv);
  P.SetDescription("Estimate EVD parameters for the given cola aligner/parameters");
  P.registerArg(alignerTypeCmd);
  P.registerArg(gapOpenCmd);
  P.registerArg(mismatchCmd);
  P.registerArg(gapExtCmd);

  P.parse();

  AlignerType aType       = AlignerType(P.GetIntValueFor(alignerTypeCmd));
  int         gapOpenPen  = P.GetIntValueFor(gapOpenCmd);
  int         mismatchPen = P.GetIntValueFor(mismatchCmd);
  int         gapExtPen   = P.GetIntValueFor(gapExtCmd);

  vecDNAVector query, target;

  string file1 = "data/FDRTestFW1000.fasta";
  string file2 = "data/FDRTestRC1000.fasta";
  query.Read(file1);
  target.Read(file2);
  query.RemoveGaps();
  target.RemoveGaps();
  double* samples = new double[100000];
  int numOfSamples = 0;
  int i, j;
  srand((unsigned)time(0));
  Cola cola1 = Cola();
  for (i=0; i<target.isize(); i++) {
    for (j=0; j<query.isize(); j++) { //TODO
      if (i!=j && rand()%10000==1) {
        DNAVector t, q;  
        t.SetToSubOf(target[i], 0, rand()%800+200);
        q.SetToSubOf(query[j], 0, rand()%800+200);
        AlignmentCola algn = cola1.createAlignment(t, q,
           AlignerParams(-1, aType, -gapOpenPen,
             -mismatchPen, -gapExtPen));
        algn.print(0,1,cout,150);
        cout<<i<<", "<<samples[numOfSamples-1]<<endl;
        samples[numOfSamples++] = algn.getSWScore();
           }
    }
  }
  double ret_mu, ret_lambda;
  esl_gumbel_FitComplete(samples, numOfSamples, &ret_mu, &ret_lambda);
  delete[] samples;

  cout << "Mu: "<<ret_mu<<" Lambda: "<<ret_lambda<<endl;

  return 0;
}

