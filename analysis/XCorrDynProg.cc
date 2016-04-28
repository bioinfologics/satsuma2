//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif

#include "XCorrDynProg.h"

MatchDynProg::MatchDynProg(int targetStart, int targetLen) 
{
  m_data.resize(targetLen);
  m_targetStart = targetStart;
}


double MatchDynProg::PrettyPrint(const std::vector<int> & index, int start, const DNAVector & target, const DNAVector & query, int len)
{
  string t, q, m, amb;

  int i;
  double matches = 0;

  int last = -1;

  double coveredPositions = 0;
  for (i=0; i<index.size(); i++) {
    int k = start + i;
    if (index[i] == -1) {
      t += target[k];
      m += " ";
      q += "-";
      amb += target[k];
      continue;
    }
    
    if (last != -1) { 
      while (last < index[i]) {
	t += "-";
	m += " ";
	q += query[last];
	//amb += query[last];
	amb += "-";
	last++;
      }
    }
    last = index[i]+1;

    t += target[k];
    q += query[index[i]];
    coveredPositions++;

    string tmp;
    tmp += target[k];
    tmp += query[index[i]];
    amb += ResolveAmbiguous(tmp);

    double dist = DNA_Equal(target[k], query[index[i]]);
    matches += dist;
    
    string mI = " ";
    if (dist > 0.3)
      mI = ".";
    if (dist > 0.49)
      mI = "~";
    if (dist > 0.95)
      mI = "|";
			 
    m += mI;

    //if (target[k] == query[index[i]]) {
    //  m += "|";
    //  matches++;
    //} else {
    //  m += " ";
    //}
  }

  //cout << q << endl;
  //cout << m << endl;
  //cout << t << endl;

  std::vector<string> tt, qq, mm, ambamb;
  string tw, qw, mw, ambw;

  int length = strlen(m.c_str());
 
  if (coveredPositions < 0.5)
    coveredPositions = 1.;


  for (i=0; i<length; i++) {    
    if (i > 0 && i % 80 == 0) {
      tt.push_back(tw);
      qq.push_back(qw);
      mm.push_back(mw);
      ambamb.push_back(ambw);
      tw = "";
      qw = "";
      mw = "";
      ambw = "";
    }
    tw += t[i];
    qw += q[i];
    mw += m[i];
    ambw += amb[i];
  }
  tt.push_back(tw);
  qq.push_back(qw);
  mm.push_back(mw);
  ambamb.push_back(ambw);

  cout << "Identity (w/ indel count): " << 100. * matches/(double)len << " %" << endl;
  cout << "-------------------------------------------------------------------------------" << endl << endl;
  for (i=0; i<tt.size(); i++) {
    cout << qq[i] << endl;
    cout << mm[i] << endl;
    cout << tt[i] << endl;
    cout << endl;
    cout << ambamb[i] << " AMB" << endl;
    cout << "-------------------------------------------------------------------------------" << endl << endl;
  }

  return (double)matches/(double)len;
}

  
void MatchDynProg::AddMatch(const DNAVector & target, const DNAVector & query, const SingleMatch & m)
{
  //cout << "adding match from " << m.GetStartTarget() << " to " << m.GetStartTarget() + m.GetLength() << endl;
  //cout << "m_start= " << m_targetStart << endl;
  //cout << "tID=" << m.GetTargetID() << " qID=" << m.GetQueryID() << " startQ=" << m.GetStartQuery();
  //if (m.IsRC())
  //  cout << " RC" << endl;
  //else
  //  cout << " FW" << endl;



  // Add to the dyn prog structure
  int i;
  int k;
  int startTarget = m.GetStartTarget();
  int length = m.GetLength();
  if (startTarget < 0)
    startTarget = 0;
  if (startTarget + length >= target.size())
    length = target.size() - startTarget - 1;


  if (m.GetStartQuery() < 0)
    cout << "ERROR!" << endl;

  for (i=startTarget; i<startTarget + length; i++) {
    k = i-startTarget + m.GetStartQuery();
    double pen = 1.;


    
    // Scoring!!
    pen = 1. - DNA_Equal(target[i], query[k]);

    //if (target[i] == query[k])
    //pen = 0.;
 

    if (i - m_targetStart >= m_data.size()) {
      cout << "WARNING: adjusting dyn prog array size!" << endl;
      m_data.resize(i - m_targetStart + 1);
    }
    MDItemList & l = m_data[i - m_targetStart];
    l.Add(MDItem(k, pen));
  }
}

int MatchDynProg::Merge(std::vector<int> & index)
{
  int i, j, k;
  if (m_data.size() == 0)
    return -1;

  MDItemList & first = m_data[0];
  MDItemList & last = m_data[m_data.size()-1];
  // Reset scores
  for (i=0; i<first.GetCount(); i++) {
    first.Get(i).Merge(-1, 0, 0);
  }
  //cout << "Init." << endl;
  for (i=1; i<m_data.size(); i++) {
    MDItemList & one = m_data[i-1];
    MDItemList & two = m_data[i];

    //cout << "i=" << i << " items: " << two.GetCount() << endl;

    for (j=0; j<one.GetCount(); j++) {
      MDItem & from = one.Get(j);
      for (k=0; k<(int)two.GetCount(); k++) {
	MDItem & to = two.Get(k);
	int diff = to.GetIndex() - from.GetIndex() - 1;
	//cout << "to=" << to.GetIndex() << " from=" << from.GetIndex() << endl;
	double pen = 0.;
	if (diff != 0)
	  pen = 1.;
	double pen_plus = (double)diff;
	if (pen_plus < 0)
	  pen_plus = -pen_plus;
	pen += pen_plus / 2;
	to.Merge(j, from.GetScore(), pen);
	//cout << "i=" << i << " score=" << from.GetScore() << " pen=" << pen << endl;
      }
    }
  }
  
  // Now traceback:
  //cout << "Traceback." << endl;
  int best = -1;
  double bestScore = XCINFINITY;


  for (i=0; i<last.GetCount(); i++) {
    //cout << "score=" << last.Get(i).GetScore() << endl;
    //cout << i << endl;
    if (last.Get(i).GetScore() < bestScore) {
      best = i;
      bestScore = last.Get(i).GetScore();
    }
  }
  if (best < 0)
    return -1;
  index.clear();
  index.resize(m_data.size(), -1);
  //cout << "best: " << best << " " << m_data.size() << endl;
  index[m_data.size()-1] = last.Get(best).GetIndex();  
  int prev = last.Get(best).GetPrevIndex();

  
  //cout << "Clean." << endl;
  for (i=m_data.size()-2; i>=0; i--) {
    index[i] = m_data[i].Get(prev).GetIndex();
    prev = m_data[i].Get(prev).GetPrevIndex();
    //cout << "i=" << i << " index=" << index[i] << endl;
  }

  
  k = 0;
  for (i=0; i<index.size(); i++) {
    int curr = index[i];
    while (i+1<index.size() && index[i+1] <= curr) {     
      index[i+1] = -1;
      i++;
    }
  }

  //cout << "Done!" << endl;
  return index.size() + k;
}



//--------------------------------------------------------------------------------------------
double CDF(double x, double m, double s)
{
  
  double r = 0.5 * (1. + erf((x-m)/s/1.414213562));

  //cout << "x=" << x << " m=" << m << " s=" << s << " cdf=" << r << endl;
  return r;
}

double Sigma(double p, int N)
{
  return sqrt(p * (1. - p) * (double)N);
}
