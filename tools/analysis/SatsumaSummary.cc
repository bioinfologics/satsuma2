#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Median length and identity.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  svec<double> ident;
  svec<int> len;

  int sum = 0;

  while (parser.ParseLine()) {
    len.push_back(parser.AsInt(2)-parser.AsInt(1));
    ident.push_back(parser.AsFloat(6));
    sum += parser.AsInt(2) - parser.AsInt(1);
  }

  cout << "Total sequence:  " << sum << endl;
  cout << "Median length:   " << len[len.isize()/2] << endl;
  cout << "Median identity: " << ident[ident.isize()/2] << endl;

  return 0;
}
