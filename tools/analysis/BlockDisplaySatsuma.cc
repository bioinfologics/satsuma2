
#include "../base/FileParser.h"
#include "../base/CommandLineParser.h"

#include "../base/SVector.h"
#include "DNAVector.h"



class SatsumaCoords
{
public:
  SatsumaCoords() {
    m_ident = 0.;
    m_startTarget = 0;
    m_startQuery = 0;
    m_endTarget = 0;
    m_endQuery = 0;

  }

  void Set(const string & target,
	   const string & query,
	   const string & ori,
	   int startTarget,
	   int endTarget,
	   int startQuery,
	   int endQuery,
	   double ident,
	   int scale) {

    m_target = target;
    m_query = query;
    m_ori = ori;
    m_ident = ident;
    m_startTarget = startTarget / scale;
    m_startQuery = startQuery / scale;
    m_endTarget = endTarget / scale;
    m_endQuery = endQuery / scale;    
  }


  bool operator < (const SatsumaCoords & c) const {
    if (m_target != c.m_target)
      return (m_target < c.m_target);
    if (m_ori != c.m_ori)
      return (m_ori < c.m_ori);
    if (m_query != c.m_query)
      return (m_query < c.m_query);
  
    return (m_startTarget < c.m_startTarget);

  }


  const string & Target() const {return m_target;}
  const string & Query() const {return m_query;}
  const string & Ori()  const {return m_ori;}
  double Identity()    const {return m_ident;}
  int StartTarget()   const {return m_startTarget;}
  int StartQuery()   const {return m_startQuery;}
  int EndTarget()   const {return m_endTarget;}
  int EndQuery()   const {return m_endQuery;}
 

  bool IsCont(const SatsumaCoords & s) const {
    if (s.Target() != m_target)
      return false;
    if (s.Query() != m_query)
      return false;
    if (s.Ori() != m_ori)
      return false;

    if (m_ori == "+") {
      if (s.StartQuery() > m_startQuery || s.EndQuery() > m_startQuery) {
	int diff = s.StartQuery() - m_startQuery - (s.StartTarget() - m_startTarget);

	if (diff < 100000 && diff > -100000)
	  return true;
	else
	  return false;
      }
    } else {

      if (s.StartQuery() < m_startQuery || s.EndQuery() < m_startQuery) {
	int diff = s.StartQuery() - m_startQuery - (s.StartTarget() - m_startTarget);
	
	if (diff < 100000 && diff > -100000)
 	  return true;
	else
	  return false;
      }
      
    }
    return false;
  }

  void Print() {
    /*cout << m_target << "\t" << m_startTarget << "\t" << m_endTarget << "\t"  << "+1\t" << m_ident << "\t";
    cout << m_query << "\t" << m_startQuery << "\t" << m_endQuery;
    if (m_ori == "+")
      cout << "\t+1\t";
    else
      cout << "\t-1\t";
      cout << m_ident << endl;*/


    cout << m_startTarget << "\t" << m_endTarget << "\t"  << "+1\t";
    cout << "\t" << m_startQuery << "\t" << m_endQuery;
    if (m_ori == "+")
      cout << "\t+1\t";
    else
      cout << "\t-1\t";
    cout << m_ident << endl;
  }
  
private:
  string m_target;
  string m_query;
  string m_ori;
  double m_ident;
  int m_startTarget;
  int m_startQuery;
  int m_endTarget;
  int m_endQuery;
  
 
};




int main( int argc, char** argv )
{
  commandArg<string> aStringCmmd("-i","satsuma summary file");
  commandArg<string> tCmmd("-t","target fasta file");
  commandArg<string> qCmmd("-q","query fasta file");
  commandArg<int> minCmmd("-min","minimum block size", 3);
  commandArg<int> minScaff("-s","minimum scaffold size", 100000);
  commandArg<bool> transCmmd("-transpose","switch query and target", false);
  commandLineParser P(argc,argv);
  P.SetDescription("Takes a satsuma summary file and writes displayable blocks.");
  P.registerArg(aStringCmmd);
  P.registerArg(tCmmd);
  P.registerArg(qCmmd);
  P.registerArg(minCmmd);
  P.registerArg(minScaff);
  P.registerArg(transCmmd);
  //P.registerArg(bStringCmmd);

  P.parse();

  string in = P.GetStringValueFor(aStringCmmd);  
  int min_block = P.GetIntValueFor(minCmmd);  
  int min_scaff = P.GetIntValueFor(minScaff);  
  bool transpose = P.GetBoolValueFor(transCmmd);  
  //String bString = P.GetStringValueFor(bStringCmmd);  
  
  //comment. ???
  FlatFileParser parser;
  
  parser.Open(in);

  int scale = 10;
  vector<SatsumaCoords> coords;

  vecDNAVector target, query;
  target.Read(P.GetStringValueFor(tCmmd));
  query.Read(P.GetStringValueFor(qCmmd));

  
  int i;
  cout << "genome\tGenome1" << endl;
  int n = 0;
  for (i=0; i<target.isize(); i++) {
    if (target[i].isize() > min_scaff) {
      n++;
    }
  }
  cout << "chromosomes\t" << n << endl;
  for (i=0; i<target.isize(); i++) {
    if (target[i].isize() > min_scaff) {
      cout << target.NameClean(i) << "\t" << target[i].isize()/scale << endl;
    }
  }
  cout << endl;
  
  cout << "genome\tGenome2" << endl;
  n = 0;
  for (i=0; i<query.isize(); i++) {
    if (query[i].isize() > min_scaff) {
      n++;
    }
  }
  cout << "chromosomes\t" << n << endl;
  for (i=0; i<query.isize(); i++) {
    if (query[i].isize() > min_scaff) {
      cout << query.NameClean(i) << "Q\t" << query[i].isize()/scale << endl;
    }
  }
  cout << endl;
 
  int k = 0;
  while (parser.ParseLine()) {
    if (parser.GetItemCount() == 0)
      continue;
    
    if (k >= (int) coords.size()) {
      coords.resize(k + 100000);
    }


    if (transpose) {
      coords[k].Set(parser.AsString(3),
		    parser.AsString(0),
		    parser.AsString(7),
		    parser.AsInt(4),
		    parser.AsInt(5),
		    parser.AsInt(1),
		    parser.AsInt(2),
		    parser.AsFloat(6),
		    scale);
    } else {

      coords[k].Set(parser.AsString(0),
		    parser.AsString(3),
		    parser.AsString(7),
		    parser.AsInt(1),
		    parser.AsInt(2),
		    parser.AsInt(4),
		    parser.AsInt(5),
		    parser.AsFloat(6),
		    scale);
    }

    k++;
  }

  coords.resize(k);


  //cout << "Before sort: " << endl;
  if (transpose)
    sort(coords.begin(),coords.end());
  //cout << "after: " << coords.isize() << endl; 


  
  int j = 0;
  int l;
  n = 0;
  for (i=0; i<k; i++) {

    const SatsumaCoords & first = coords[i];
    double ident = first.Identity();
    for (j=i+1; j<k; j++) {
      const SatsumaCoords & next = coords[j];
      if (j+1 == k || !first.IsCont(coords[j])) {
	break;
      }
      ident += coords[j].Identity();
    }


    if (j-i >= min_block) {

      int l1 = (target(coords[i].Target())).isize();
      int l2 = (query(coords[i].Query())).isize();

      if (l1 > min_scaff && l2 > min_scaff) {
	cout << /*"block:" << n << "\t" <<*/ coords[i].Target() << "\t" << coords[i].StartTarget();
	cout << "\t" << coords[j-1].EndTarget() << "\t" << "+1\t" /*<< ident / (double)(j-i) << "\t"*/; 
	cout << coords[i].Query() << "Q\t";
	
	if (coords[i].Ori() == "+")
	  cout << coords[i].StartQuery() << "\t" << coords[j-1].EndQuery() << "\t+1\t";
	else
	  cout << coords[j-1].StartQuery() << "\t" << coords[i].EndQuery() << "\t-1\t";
	cout << ident / (double)(j-i) << "\t" << j-i << endl;
	
	n++;
	for (l=i; l<j; l++) {
	  coords[l].Print();
	}
	cout << endl;
      }
    }

    i = j-1;
    
  }




  return 0;
}
