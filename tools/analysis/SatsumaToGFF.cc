#include <string>
#include "../base/CommandLineParser.h"

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input satsuma file");
  commandArg<bool> tCmmd("-t","output coordinates based on target", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Generates a GFF3 file from a satsuma file.");
  P.registerArg(fileCmmd);
  P.registerArg(tCmmd);

  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  bool bTarget = P.GetBoolValueFor(tCmmd);

  FlatFileParser parser;
  
  parser.Open(fileName);

  int count = 1;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (bTarget) {
      // target based
      cout << parser.AsString(0) << "\tSatsuma\tblock\t" << parser.AsInt(1) << "\t" << parser.AsInt(2) << "\t\t" << parser.AsString(7) << "\t.\t";
      cout << "ID=sblock_" << to_string(count) << ";Name=sblock_" << to_string(count) << ";Target=" << parser.AsString(3) << endl;
    } else {
      // query based
      cout << parser.AsString(3) << "\tSatsuma\tblock\t" << parser.AsInt(4) << "\t" << parser.AsInt(5) << "\t\t" << parser.AsString(7) << "\t.\t";
      cout << "ID=sblock_" << to_string(count) << ";Name=sblock_" << to_string(count) << ";Target=" << parser.AsString(0) << endl;
    }

    count++;

  }
  
  return 0;
}
