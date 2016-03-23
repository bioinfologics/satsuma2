#include <string>
#include "../base/CommandLineParser.h"
#include "../base/FileParser.h"
#include "DNAVector.h"


void Print(const string &s) {
  int n = strlen(s.c_str());
  for (int i=0; i<n; i++) {
    if (s[i] != '-')
      cout << s[i];
  }
  cout << endl;
}

int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","satsuma file");
  commandArg<string> fastaCmmd("-f","fasta file to get sequence from");
  commandArg<bool> tCmmd("-t","output target sequence", false);
  commandArg<int> slackCmmd("-s","flanks on each side", 0);
  commandLineParser P(argc,argv);
  P.SetDescription("Generates a FASTA file from a satsuma summary file.");
  P.registerArg(fileCmmd);
  P.registerArg(tCmmd);
  P.registerArg(fastaCmmd);
  P.registerArg(slackCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  bool bTarget = P.GetBoolValueFor(tCmmd);
  int slack = P.GetIntValueFor(slackCmmd);
  
  vecDNAVector dna;
  dna.Read(fastaName);

  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    string chr;
    int start = 0;
    int stop = 0;
    bool bRC = false;
    if (bTarget) {
      chr = parser.AsString(0);
      start = parser.AsInt(1);
      stop = parser.AsInt(2);
    } else {
      chr = parser.AsString(3);
      start = parser.AsInt(4);
      stop = parser.AsInt(5);
      if (parser.AsString(7) == "-")
	bRC = true;
    }
    const DNAVector & d = dna(chr);
    start -= slack;
    if (start < 0)
      start = 0;
    stop += slack;
    if (stop >= d.isize())
      stop = d.isize();

    cout << ">" << chr << ":" << start << "-" << stop << endl;
    int k = 0;

    DNAVector small;
    small.SetToSubOf(d, start, stop-start);
    if (bRC)
      small.ReverseComplement();
    for (int i=0; i<small.isize(); i++) {
      cout << small[i];
      if (k > 0 && k % 80 == 0)
	cout << endl;
      k++;
    }
    cout << endl;
  }
  return 0;
}
