#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"




int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<int> nCmmd("-n","how many");
  commandLineParser P(argc,argv);
  P.SetDescription("Splits a Satsuma file into equal chunks.");
  P.registerArg(fileCmmd);
  P.registerArg(nCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  
  int n = P.GetIntValueFor(nCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);


  svec<string> all;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    all.push_back(parser.Line());
  }
  int cc = all.isize() / n;
  
  int k = 0;
  FILE * pOut = NULL;
  
  for (int i=0; i<all.isize(); i++) {
    if ((i % cc) == 0) {
      char num[256];
      sprintf(num, "%d", k);
      k++;
      string out = fileName;
      out += ".";
      if (strlen(num) == 1)
	out += "0";
      out += num;
      if (pOut != NULL)
	fclose(pOut);
      pOut = fopen(out.c_str(), "w");
      cout << "ColaAlignSatsuma -s " << out << " -t target -q query -o cola/ColaAligns" << num << ".out &" << endl;
    }
    fprintf(pOut, "%s\n", all[i].c_str()); 
  }
  fclose(pOut);
  return 0;
}
