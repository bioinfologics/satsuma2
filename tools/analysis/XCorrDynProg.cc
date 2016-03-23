//#ifndef FORCE_DEBUG
//#define NDEBUG
//#endif

#include "analysis/XCorrDynProg.h"
#include <math.h>
#include "analysis/Blosum.h"
#include "analysis/MultiProtein.h"

MatchDynProg::MatchDynProg(int targetStart, int targetLen) 
{
  m_data.resize(targetLen);
  m_targetStart = targetStart;
}


double MatchDynProg::PrettyPrint(const svec<int> & index, int start, const DNAVector & target, const DNAVector & query, int len)
{
  string t, q, m, amb;

  int i;
  double matches = 0;

  int last = -1;

  double coveredPositions = 0;
  for (i=0; i<index.isize(); i++) {
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

  svec<string> tt, qq, mm, ambamb;
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
  for (i=0; i<tt.isize(); i++) {
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
  if (startTarget + length >= target.isize())
    length = target.isize() - startTarget - 1;


  if (m.GetStartQuery() < 0)
    cout << "ERROR!" << endl;

  for (i=startTarget; i<startTarget + length; i++) {
    k = i-startTarget + m.GetStartQuery();
    double pen = 1.;


    
    // Scoring!!
    pen = 1. - DNA_Equal(target[i], query[k]);

    //if (target[i] == query[k])
    //pen = 0.;
 

    if (i - m_targetStart >= m_data.isize()) {
      cout << "WARNING: adjusting dyn prog array size!" << endl;
      m_data.resize(i - m_targetStart + 1);
    }
    MDItemList & l = m_data[i - m_targetStart];
    l.Add(MDItem(k, pen));
  }
}

int MatchDynProg::Merge(svec<int> & index)
{
  int i, j, k;
  if (m_data.isize() == 0)
    return -1;

  MDItemList & first = m_data[0];
  MDItemList & last = m_data[m_data.isize()-1];
  // Reset scores
  for (i=0; i<first.GetCount(); i++) {
    first.Get(i).Merge(-1, 0, 0);
  }
  //cout << "Init." << endl;
  for (i=1; i<m_data.isize(); i++) {
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
  index.resize(m_data.isize(), -1);
  //cout << "best: " << best << " " << m_data.isize() << endl;
  index[m_data.isize()-1] = last.Get(best).GetIndex();  
  int prev = last.Get(best).GetPrevIndex();

  
  //cout << "Clean." << endl;
  for (i=m_data.isize()-2; i>=0; i--) {
    index[i] = m_data[i].Get(prev).GetIndex();
    prev = m_data[i].Get(prev).GetPrevIndex();
    //cout << "i=" << i << " index=" << index[i] << endl;
  }

  
  k = 0;
  for (i=0; i<index.isize(); i++) {
    int curr = index[i];
    while (i+1<index.isize() && index[i+1] <= curr) {     
      index[i+1] = -1;
      i++;
    }
  }

  //cout << "Done!" << endl;
  return index.isize() + k;
}



//--------------------------------------------------------------------------------------------

double BlosumScore(char a, char b) 
{
  if (a == b)
    return 0.;
  if (Similarity(a, b) > 0)
    return 0.75;
  return 1.;
}

bool XCDynProg::FindBestBracket(int & last, int & m, int & first)
{
  int i, j;

  int tSize = m_matrix[0].Size();

  //cout << "Traceback, finding best bracket." << endl;
  //int last = -1;
  double bestTotal = 9999999999999999999.;

  svec<double> minScores;
  minScores.resize(tSize, bestTotal);

  for (i=0; i<tSize; i++) {
    for (j=0; j<m_matrix.isize(); j++) {
      if (m_matrix[j].Score(i) < minScores[i]) {
	//cout << "i=" << i << " j=" << j << " score=" << m_matrix[j].Score(i) << endl;
	minScores[i] = m_matrix[j].Score(i);
      }
    }
  }

  int k;

  int max = 0;

  last = 0;
  m = -1;
  first = 0;

  //double bestScore = 

  
  for (i=tSize-1; i>=1; i--) {    
    for (j=0; j<m_matrix.isize(); j++) {
      double score = m_matrix[j].Score(i);

      if (score == minScores[i] && m_matrix[j].Score(i-1) == score) {
	//cout << "First: " << i << endl;
	first = i;
      }

      //cout << "Score: " << score << " Prev=" << m_matrix[j].Score(i-1) << " min=" << minScores[j] << endl;
      //cout << "Consider: " << score << " m " << m << "minScore " << minScores[i] << " mat " << m_matrix[j].Score(i-1) << endl;
      if (m == -1 && score == minScores[i] && m_matrix[j].Score(i-1) == score) {
	//cout << "Last: " << i << endl;
	last = i;
	m = j;
      }
    }
  }

      /*
      if (score > bestTotal/2)
	continue;

      for (k=i-25; k>=0; k--) {
	//if (i == 719 && j == 0) {
	//  cout << "Score: " << i-k << "\t" << score - minScores[k] << endl;
	//}

	if (minScores[k] >= bestTotal/2)
	  continue;
	double diff = score - minScores[k];
	double s = NormScore(diff, i-k);
	//cout << "Normscore: " << s << " thresh=" << m_thresh << endl;
	if (s < bestTotal) {
	  bestTotal = s;

	//if (i-k > max) {
	  max = i-k;
	  last = i;
	  m = j;
	  first = k;
	  //cout << "BEST!" << endl;
	}
	if (s > m_thresh)
	  break;
	//}
      }
    }
    }*/


  // cout << "Best last node: " << last << " matrix: " << m << " first=" << first << endl;
  
  return false;
}


double ComputePrint(const string & qq, const string & tt, bool bPrint, int & len)
{
  //bPrint = true;

  int i, j;
  int l = strlen(qq.c_str());
  int k = 0;

  double matches = 0;
  do {
    //cout << endl;
    for (i=0; i<80; i++) {
      if (bPrint)
	cout << "-";
    }
    if (bPrint)
      cout << endl;

    for (i=k; i<l; i++) {
      if (i == k+80)
	break;
      if (bPrint)
	cout << qq[l-1-i];
    }
    if (bPrint)
      cout << endl;
    
    for (i=k; i<l; i++) {
      if (i == k+80)
	break;
      if (tt[l-1-i] == qq[l-1-i]) {
	if (bPrint)
	  cout << tt[l-1-i];
	matches += 1.;
      } else {
	if (Similarity(tt[l-1-i], qq[l-1-i]) > 0) {
	  if (bPrint)
	    cout << "+";
	  matches += 0.3;
	} else {
	  if (bPrint)
	    cout << " "; 
	}
      }
    }
    if (bPrint)
      cout << endl;
    
    for (i=k; i<l; i++) {
      if (i == k+80)
	break;
     if (bPrint)
       cout << tt[l-1-i];
    }
    if (bPrint)
      cout << endl;
    k += 80;
  } while (k < l);
  if (bPrint)
    cout << endl;

  len = l;

  if (l == 0)
    return 0.;

  if (matches < 20)
    return 0.;

  cout << "Matches; " << matches << " l: " << l << endl;
  return matches / (double)l;

}

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


double GetPValue(double ident,
		 int length,
		 double expect_ident)
{

      
  double p_match = expect_ident;
     
  //cout << "ident=" << ident << " p_match=" << p_match << endl;


  double s = Sigma(p_match, length);

  //cout << "Sigma=" << s << endl; 
      
  double m = p_match * (double)length;

  double x = (double)length * ident;
    
  //cout << "m=" << m << " x=" << x << endl;

  double cdf = CDF(m, x, s);

  cdf = erfc((x-m)/s/1.414213562)/2;


  //cout << "cdf=" << cdf << endl;

  double targetSize = 1;

  double expect = cdf * targetSize; // Twice the genome size if we do RC!!
  double p_val = -expm1(-expect);
     
  //cout << "expect=" << expect << " p-val=" << p_val << endl;
  return p_val;
}





double XCDynProg::Align(const MultiProtein & target, const MultiProtein &query, const SignalFilter & filter)
{
  int i, j, k;

  m_matrix.resize(filter.GetTopCount());


  //cout << "Setting up structures" << endl;
  for (i=0; i<filter.GetTopCount(); i++) {
    //cout << " -> " << i << " shift=" << filter.GetShift(i) << " " << m_matrix.isize() << endl;
    
    m_matrix[i].SetUp(query.Sequence(), filter.GetShift(i), target.isize());
  }

  //cout << "Dynprog." << endl;
  for (i=1; i<target.isize(); i++) {
    for (j=0; j<m_matrix.isize(); j++) {
      if (m_matrix[j].Letter(i-1) == 0)
	continue;
      for (k=0; k<m_matrix.isize(); k++) {

	//?????????????????????????
	//if (m_matrix[k].Letter(i-1) == 0)
	//continue;
	
	int jump = m_matrix[j].Shift() - m_matrix[k].Shift();

	double trans = jump;
	if (trans < 0)
	  trans = -trans;
	double pen = 2.;


	int posK = i+m_matrix[k].Shift();

	//if (m_matrix[k].Letter(posK) == 0)
	//continue;
	if (m_matrix[k].Letter(i) == 0)
	  continue;


	//pen *= BlosumScore(target[i], m_matrix[k].Letter(i));
	pen *= target.Distance(query, i, posK);
      

	double score = pen + trans + m_matrix[j].Score(i-1);
	
       
	//if (target[i] == query[posK]) {
	//  cout << "Compare i=" << i << " " << target[i] << " j=" << j <<  " k=" << k << " jump=" << jump << " " << query[posK] << endl;
	//  cout << "Score=" << score << " curr=" << m_matrix[k].Score(i) << endl;
	//}


	if (score < m_matrix[k].Score(i)) {
	  //cout << "Transition: " << j << " -> " << k << " i=" << i << endl;
	  m_matrix[k].SetScore(i, score);
	  m_matrix[k].SetBack(i, j);
	}
	//if (score < m_matrix[k].Score(jump)) {
	//  m_matrix[k].SetScore(jump, score);
	//  m_matrix[k].SetBack(jump, j);
	//}
      }
    }
    
  }

  /*
  for (i=0; i<target.isize(); i++) {
    cout << "i=" << i << "\t";
    for (j=0; j<5; j++) {
      cout << m_matrix[j].Score(i) << " (" << m_matrix[j].Back(i) << ")\t";
    }
    cout << endl;
  }
  */

  /*
  cout << "Traceback." << endl;
  // Traceback...
  int last = -1;
  double bestTotal = 9999999999999999999.;
  for (i=target.isize()-1; i>=0; i--) {
    
    for (j=0; j<m_matrix.isize(); j++) {
      if (m_matrix[j].Score(i) < bestTotal) {
	bestTotal = m_matrix[j].Score(i);
	last = j;
      }
    }
    if (last != -1)
      break;
  }
  cout << "Best last node: " << i << " matrix: " << last << " score=" << bestTotal << endl;
  */
  //cout << "Best last node: " << i << " matrix: " << last << " score=" << bestTotal << endl;

  //cout << "BRACKET!" << endl;
  int last, mIndex, first;
  FindBestBracket(last, mIndex, first);
  //cout << "Best last node: " << i << " matrix: " << last << " index: " << mIndex << endl;

  //cout << "Done." << endl;

  string tt;
  string qq;

  i = last;

  //mIndex = 1;

  int lastQuery = 0;
  int firstQuery = 99999999;

  do {
    //cout << "i=" << i << " mIndex=" << mIndex << endl;
    const XCDynProgLine & l = m_matrix[mIndex];
    if (i < first)
      break;
    int xx = i;
    int back = l.Back(i);
    
    if (back != mIndex) {
      //cout << "Jump, i=" << i << " back=" << back << endl;
    }
    if (back < 0)
      break;
    int diff = l.Shift() - m_matrix[back].Shift();
    if (back != mIndex) {
      //cout << "Found shift: " << diff << " i=" << i << endl;
      if (diff < 0) {
	for (j=0; j<-diff; j++) {
	  if (i < 0)
	    break;
	  //qq += l.Letter(i-j);
	  //tt += "-";
	  tt += target[i];
	  qq += "-";
	  i--;
	}
      }
    }
    if (i < 0)
      break;
    tt += target[i];
    qq += l.Letter(xx);
    int qPos = xx + l.Shift();
    if (qPos > lastQuery)
      lastQuery = qPos;
    if (qPos < firstQuery)
      firstQuery = qPos;

    if (back != mIndex) {
      if (diff > 0) {
	//const XCDynProgLine & bb = m_matrix[back];    
	for (j=0; j<diff; j++) {
	  if (xx-j-1 < 0)
	    break;
	  qq += l.Letter(xx-j-1);
	  tt += "-";
	  //cout << "  insert i=" << xx-j-1 << " letter=" << l.Letter(xx-j-1) << endl;
	  //i--;
	  //tt += target[i];
	  //qq += "-";
	}
      }
    } else {
      //tt += target[i];
      //qq += l.Letter(xx);      
    }


    //cout << "i=" << xx << " letter=" << l.Letter(xx) << "  " << mIndex << endl;
    


    mIndex = back;
    i--;

  } while (mIndex != -1);


  

  //cout << "Printing." << endl;
  int len;
  double ident = ComputePrint(qq, tt, false, len);
  double p_value = GetPValue(ident, len, m_expect);
  if (strlen(qq.c_str()) < 30 || p_value > 0.00001) {
  //if (strlen(qq.c_str()) < 30 || ident < 0.3) {
    //cout << "No alignment." << endl;
    return 1.;
  }
  
  cout << "Target: " << first << "-" << last << " Query: " << firstQuery << "-" << lastQuery << endl;

  ComputePrint(qq, tt, true, len);

  cout << "Identity: " << 100. * ident << " % p=" << p_value << endl<< endl;
  return p_value;
}
