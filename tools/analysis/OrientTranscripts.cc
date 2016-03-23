#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "analysis/DNAVector.h"



int main( int argc, char** argv )
{

  commandArg<string> fileCmmd("-i","input file");
  commandArg<string> fastaCmmd("-f","fasta file");
  commandArg<string> outCmmd("-o","out fasta file");
  commandLineParser P(argc,argv);
  P.SetDescription("Orients sequences according to alignments.");
  P.registerArg(fileCmmd);
  P.registerArg(fastaCmmd);
  P.registerArg(outCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  string fastaName = P.GetStringValueFor(fastaCmmd);
  string outName = P.GetStringValueFor(outCmmd);
  vecDNAVector dna, out;
  dna.Read(fastaName);


  //comment. ???
  FlatFileParser parser;
  
  parser.Open(fileName);

  svec<string> name;
  svec<string> ori;
  svec<int> start;
  svec<int> stop;

  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    if (parser.AsFloat(6) < 0.7)
      continue;
    if (parser.AsInt(2) - parser.AsInt(1) < 200)
      continue;

    name.push_back(parser.AsString(3));
    ori.push_back(parser.AsString(7));
    start.push_back(parser.AsInt(4));
    stop.push_back(parser.AsInt(5));

  }

  int i, j;
  for (i=0; i<dna.isize(); i++) {
    string o;
    string n = dna.NameClean(i);
    //bool b = true;

    int fromF = 0; 
    int toF = 0;
    int fromR = 0; 
    int toR = 0;
    

    for (j=0; j<name.isize(); j++) {
      if (name[j] == n) {
	//if (dna.Name(i)== ">sample22_9")
	//cout << "===========> " << ori[j] << endl;
	/*if (o == "") {
	  o = ori[j];
	} else {
	  if (o != ori[j])
	    b = false;
	    }*/
	if (ori[j] == "+") {
	  if (stop[j] - start[j] > toF - fromF) {
	    fromF = start[j];
	    toF = stop[j];
	  } 
	} else {
	  if (stop[j] - start[j] > toR - fromR) {
	    fromR = start[j];
	    toR = stop[j];
	  } 
	}

      }
    }
    cout << "+ " << fromF << " " << toF << " - " << fromR << " " << toR << endl;
    if (toF - fromF < 200 && toR - fromR < 200) {
      cout << "Skipping: " << dna.Name(i) << endl;
      continue;
    }
    DNAVector tmp;
    if (toR - fromR  > toF - fromF) {
      //dna[i].ReverseComplement();
      if (toF - fromF > 100) {
	tmp.SetToSubOf(dna[i], fromR, toR-fromR);
      } else {
	tmp = dna[i];
      }
      tmp.ReverseComplement();
      cout << "Reversing " << dna.Name(i) << endl;
    } else {
      if (toR - fromR > 100) {
	tmp.SetToSubOf(dna[i], fromF, toF-fromF);
      } else {
 	tmp = dna[i];
      }
      cout << "Keeping " << dna.Name(i) << endl;
    }
    out.push_back(tmp, dna.Name(i));

  }


  out.Write(outName);

  return 0;
}
