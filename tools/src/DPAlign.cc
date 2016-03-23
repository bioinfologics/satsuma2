#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "src/DPAlign.h"
#include "analysis/AlignProbability.h"
#include <iostream>
#include <fstream>
#include <string.h>

void FullAlignment::WriteBinary(ofstream &oStrm)
{
  int n = m_nameA.size();
  oStrm.write((char *) &n, sizeof(int));
  oStrm.write(m_nameA.c_str(), n);
    
  n=m_nameB.size();
  oStrm.write((char *) &n, sizeof(int));
  oStrm.write(m_nameB.c_str(), n);

  n=m_indexA.size();

  oStrm.write((char *) &n, sizeof(int));
  oStrm.write((char *) &m_direction[0], n*sizeof(int));
  oStrm.write((char *) &m_question[0], n*sizeof(int));
  oStrm.write((char *) &m_indexA[0], n*sizeof(int));
  oStrm.write((char *) &m_indexB[0], n*sizeof(int));

  n = m_offsetA;
  oStrm.write((char *) &n, sizeof(int)); 
  n = m_offsetB;
  oStrm.write((char *) &n, sizeof(int)); 

}


void FullAlignment::ReadBinary(ifstream &iStrm)
{
  
  int n;
  iStrm.read((char *) &n, sizeof(int));
 
  vector<char> buff(n);
  iStrm.read(&buff[0], n);

  m_nameA  = string(&buff[0],n);

  iStrm.read((char *) &n, sizeof(int));

//   buff.clear();
//   buff.resize(n);
  vector<char> buff2(n);
  iStrm.read(&buff2[0],n);

  m_nameB = string(&buff2[0],n);

  iStrm.read((char *) &n, sizeof(int));

  m_direction.clear();
  m_question.clear();
  m_indexA.clear();
  m_indexB.clear();

  m_direction.resize(n);
  m_question.resize(n);
  m_indexA.resize(n);
  m_indexB.resize(n);


  iStrm.read((char *) &m_direction[0], n*sizeof(int));
  iStrm.read((char *) &m_question[0], n*sizeof(int));
  iStrm.read((char *) &m_indexA[0], n*sizeof(int));
  iStrm.read((char *) &m_indexB[0], n*sizeof(int));

  iStrm.read((char *) &n, sizeof(int));
  m_offsetA = n;

  iStrm.read((char *) &n, sizeof(int));
  m_offsetB = n;

}


double FullAlignment::Identity(const DNAVector & a, 
			       const DNAVector & b,
			       bool rcExplicit,
			       int offA,
			       int offB,
			       int printOff)
{
  int i, j;

  //int n = 80;
  int totalMatch = 0; 
  string seqA, seqB, match, dir;
  //cout << "Size=" << Size() << endl;

  int k = printOff;
  int ret = 0;
  int totalSize = a.isize();


  for (i=0; i<Size(); i++) {
    //cout << "*** i=" << i << " idxB=" << IndexB(i) << endl;
    int d = Direction(i);
    if (rcExplicit) {
      dir += "<";
    } else {
      if (d == -1)
	dir += "<";
      if (d == 1)
	dir += ">";
      if (d == 0)
	dir += " ";
    }
    char lettA = '-';
    char lettB = '-';
    if (IndexA(i) != -1)
      lettA = a[IndexA(i)+offA];

    if (IndexB(i) != -1) {
      if (d == 1)
	lettB = b[IndexB(i)+offB];
      else
	lettB = GetRC(b[IndexB(i)+offB]);
    }
    if (lettA == lettB) {
      match += "|";
      totalMatch++;
    } else {
      match += " ";
    }
    /*
    if (Quest(i) > 20) {
      if (lettA != '-')
	lettA = (char)tolower(lettA);
      if (lettB != '-')
	lettB = (char)tolower(lettB);
	}*/


    seqA += lettA;
    seqB += lettB;

    /*
    if ((k+1) % n == 0 || i+1 == Size()) {
      for (j=0; j<n; j++) {
	cout << "-";
      } 
      cout << endl;

      cout << dir << endl;
      cout << seqB << endl;
      cout << match << endl;
      cout << seqA << endl;
      ret = strlen(seqA.c_str());
      if (ret == n)
	ret = 0;

      dir = "";
      seqA = "";
      seqB = "";
      match = "";

      }*/
    k++;

  }
  return (double)totalMatch / (double)totalSize;
}
int FullAlignment::PrettyPrint(const DNAVector & a, 
			       const DNAVector & b,
			       bool rcExplicit,
			       int offA,
			       int offB,
			       int printOff)
{
  int i, j;

  int n = 80;
  int totalMatch = 0; 
  string seqA, seqB, match, dir;
  //cout << "Size=" << Size() << endl;

  int k = printOff;
  int ret = 0;

  for (i=0; i<Size(); i++) {
    //cout << "*** i=" << i << " idxB=" << IndexB(i) << endl;
    int d = Direction(i);
    if (rcExplicit) {
      dir += "<";
    } else {
      if (d == -1)
	dir += "<";
      if (d == 1)
	dir += ">";
      if (d == 0)
	dir += " ";
    }
    char lettA = '-';
    char lettB = '-';
    if (IndexA(i) != -1)
      lettA = a[IndexA(i)+offA];

    if (IndexB(i) != -1) {
      if (d == 1)
	lettB = b[IndexB(i)+offB];
      else
	lettB = GetRC(b[IndexB(i)+offB]);
    }
    if (lettA == lettB) {
      match += "|";
      totalMatch++;
    } else {
      match += " ";
    }
    /*
    if (Quest(i) > 20) {
      if (lettA != '-')
	lettA = (char)tolower(lettA);
      if (lettB != '-')
	lettB = (char)tolower(lettB);
	}*/


    seqA += lettA;
    seqB += lettB;

    if ((k+1) % n == 0 || i+1 == Size()) {
      for (j=0; j<n; j++) {
	cout << "-";
      } 
      cout << endl;

      cout << dir << endl;
      cout << seqB << endl;
      cout << match << endl;
      cout << seqA << endl;
      ret = strlen(seqA.c_str());
      if (ret == n)
	ret = 0;

      dir = "";
      seqA = "";
      seqB = "";
      match = "";

    }
    k++;

  }
  return totalMatch;
}

//=================================================
void AlignsPrinter::PrettyPrint(const FullAlignment & align,
				const DNAVector & a, 
				const DNAVector & b,
				int offA,
				int offB,
				bool bNew)
{
  int i, j;

  int n = 80;

  //string seqA, seqB, match, dir;
  //cout << "Size=" << Size() << endl;

  //int k = m_k;
  int ret = 0;

  if (bNew && strlen(dir.c_str()) > 0) {
    for (j=0; j<n; j++) {
      cout << "-";
    } 
    cout << endl;
    
    cout << dir << endl;
    cout << seqB << endl;
    cout << match << endl;
    cout << seqA << endl;
    
    dir = "";
    seqA = "";
    seqB = "";
    match = "";
    k = 0;

  }


  for (i=0; i<align.Size(); i++) {
    //cout << "*** i=" << i << " idxB=" << IndexB(i) << endl;
    int d = align.Direction(i);
    if (d == -1)
      dir += "<";
    if (d == 1)
      dir += ">";
    if (d == 0)
      dir += " ";

    char lettA = '-';
    char lettB = '-';
    if (align.IndexA(i) != -1)
      lettA = a[align.IndexA(i)+offA];

    if (align.IndexB(i) != -1) {
      if (d == 1)
	lettB = b[align.IndexB(i)+offB];
      else
	lettB = GetRC(b[align.IndexB(i)+offB]);
    }
    if (lettA == lettB)
      match += "|";
    else
      match += " ";

    /*
    if (Quest(i) > 20) {
      if (lettA != '-')
	lettA = (char)tolower(lettA);
      if (lettB != '-')
	lettB = (char)tolower(lettB);
	}*/


    seqA += lettA;
    seqB += lettB;

    if ((k+1) % n == 0 || i+1 == align.Size()) {
      ret = strlen(seqA.c_str());
      if (ret == n) {
	ret = 0;
	for (j=0; j<n; j++) {
	  cout << "-";
	} 
	cout << endl;
	
	cout << dir << endl;
	cout << seqB << endl;
	cout << match << endl;
	cout << seqA << endl;

	dir = "";
	seqA = "";
	seqB = "";
	match = "";
      }

    }
    k++;

  }
  k = ret;
}




//=======================================================================
void DPAligner::SetUpColumn(int x_off, int len)
{
  m_x_off = x_off;
  m_columns.resize(len);
}

void DPAligner::SetUpRow(int x, int y_off, int len)
{
  int i = x - m_x_off;
  m_columns[i].resize(len);
  m_columns[i].Set(y_off);
  //if (i == 0)
  //m_columns[i].SetScore(0.);
  if (i == 0) {
    m_columns[i].SetScore(0., 0);
  }
}

void DPAligner::SetRect(int x_len, int y_len)
{
  SetUpColumn(0, x_len);
  int i;
  for (i=0; i<x_len; i++) {
    SetUpRow(i, 0, y_len);
  }

}


void DPAligner::SetReward(const DNAVector & colBases, 
			  const DNAVector & rowBases) {

  int i, j, k;
  if (m_columns.isize() != colBases.isize()) {
    cout << "ERROR: SetReward does not work!!" << endl;
  }

  int win = 50;

  for (i=0; i<m_columns.isize(); i++) {
    const DPAlignNodeColumn & c = m_columns[i];
    if (c.isize() != rowBases.isize()) {
      cout << "ERROR: SetReward will not work!!" << endl;
    }
    
    for (j=c.isize()-1; j>=0; j--) {
      // Slow, slow!!
      int n = i;
      double count = 0.;
      for (k=j; k>=0; k--) {
	//cout << "k=" << k << " n=" << n << endl;
	count += DNA_Equal(colBases[n], rowBases[k]);
       
	n--;
	if (n < 0)
	  break;
	if (i-n == win)
	  break;
      }
      if (n >=0 && k >= 0 && count > (double)win - 10.) {
	n++;
	DPAlignNodeColumn & cc = m_columns[n];
	//cout << "Setting reward to i=" << n << " j=" << k << endl;
	cc[k].SetReward(-count);
      }
    }
  } 
  
}


void DPAligner::Align(FullAlignment & align,
		      const DNAVector & colBases, 
		      const DNAVector & rowBases)
{
  int i, j, k, l;

  int loopsize = 0;
  int hiY = 0;
  
  for (i=0; i<m_columns.isize(); i++) {
    const DPAlignNodeColumn & c = m_columns[i];
    int hi = c.GetY() + c.isize();
    if (hi > loopsize) {
      hiY = hi;
    }
  }
  //cout << "Maximum y=" << hiY << " Maximum x=" << m_columns.isize() << endl;

  //if (m_columns.isize() > loopsize)
  loopsize = hiY + m_columns.isize();

  //  if (m_bUseRewards)
  //SetReward(colBases, rowBases);
  

  //cout << "loopsize=" << loopsize << endl;

  for (l=0; l<loopsize; l++) {
    //cout << "processing l=" << l << endl;

    // Find the best score first
    double minScore = 999999999.;
    int bestI = -1;
    int bestJ = -1;
    i = l;
    for (j=0; j <loopsize ; j++) {      
      if (i >= m_columns.isize() - 1) {
	i--;
	continue;
      }
      if (i < 0)
	break;

      if (j >= hiY)
	continue;

      //cout << "c0 i=" << i << endl;
      const DPAlignNodeColumn & c = m_columns[i];      
      
      double sFW = c[j].GetScore();
      //cout << "c0 done" << endl;

      if (sFW < minScore) {
	minScore = sFW;
	bestI = i;
	bestJ = j;
      }

      i--;
    }

    //cout << "l=" << l << " -> best score=" << minScore << " i=" << bestI << " j=" << bestJ << endl;
    // Now do the update
    i = l;

    //const DPAlignNodeColumn & ccc = m_columns[1696];
    //cout << "Check score=" << ccc[2689].GetScore() << " i=" << ccc[2689].GetBackCol();
    //cout << " j=" << ccc[2689].GetBackRow() << endl;

    for (j=0; j <loopsize ; j++) {
      
      if (i >= m_columns.isize() - 1) {
	i--;
	continue;
      }
      if (i < 0)
	break;

      if (j >= hiY)
	continue;

      //cout << "Processing i=" << i << " j=" << j << endl;

      //if (l >=2653)
      //cout << "c1" << endl;
      const DPAlignNodeColumn & c = m_columns[i];      
      int y = c.GetY();    
      DPAlignNodeColumn & c_next = m_columns[i+1];
      //if (l >=2653)
      //cout << "c1 done" << endl;
      
      //double localScoreFW = 1. - 3. * DNA_Equal(colBases[m_x_off + i], rowBases[y + j]);
      double localScoreFW = 1. - DNA_Equal(colBases[m_x_off + i], rowBases[y + j]);
      
      double sFW = c[j].GetScore();

      
      int col = i;
      int row = j;


      /*
      if (i>1688 && i <1900 && j > 2681 && j < 2893
	  && j - i == 993) {
	cout << "---> i=" << i << " j=" << j << " score=" << sFW << " back i=" << c[j].GetBackCol();
	cout << " j=" << c[j].GetBackRow() << endl;
	}*/
 
      //double tbScore = TraceBackScore(row, col, colBases, rowBases);
      //sFW -= tbScore;
      //cout << "i=" << i << " j=" << j << endl;
      
      //cout << "Local fw=" << localScoreFW << " rc=" << localScoreRC << " acc fw=" << sFW << " rc=" << sRC << endl;
      
      
      // Forward
      c_next.Minimize(j + y+1, sFW + localScoreFW, i, j);
      
      //cout << "Score diff=" << sFW - minScore << " skip=" << GetSkipLen(sFW - minScore) << endl;

      // Horizontal
      int to = i + GetSkipLen(sFW - minScore); // This might not be right!!

      if (to >= m_columns.isize())
	to = m_columns.isize();
      
      for (k=i+2; k<to; k++) {
	//if (l >=2653)
	//cout << "c2" << endl;
 	DPAlignNodeColumn & c_hori = m_columns[k];
	c_hori.Minimize(j + m_columns[k].GetY() + 1, sFW + localScoreFW + SkipPen(k, i), i, j);
	//if (l >=2653)
	//cout << "c2 done" << endl;
      }
      
      //Vertical
      //int from = y + j - m_maxSkip;
      to = y + j + GetSkipLen(sFW - minScore) ; 
      
      /*
	for (k=from; k<=y+j-2; k++) {
	c_next.MinimizeRC(k, sRC + localScoreRC + SkipPen(j, k), i, j);
	}*/
      
      for (k=y+j+2; k<=to; k++) {
	/*
	if (i+1>1688 && i+1 <1900 && k > 2681 && k < 2893
	    && k - (i+1) == 993) {
	  cout << "+++ Seeding " << i+1 << "/" << k << " from " << i <<"/" << j << " score=";
	  cout << sFW + localScoreFW << " skip=" <<  SkipPen(j, k) << endl; 
	}
	*/
	c_next.Minimize(k, sFW + localScoreFW + SkipPen(j, k), i, j);	
      }
  
      i--;
    }
            
  }
  TraceBack(align);

}

double DPAligner::TraceBackScore(int row, int col, 			
				 const DNAVector & colBases, 
				 const DNAVector & rowBases)
{
  int maxBack = 35;
  double sum = 0.;
  int bestRow = row;
  int bestCol = col;
  double score = 0.;
  while (bestRow != -1) {
    const DPAlignNodeColumn & item = m_columns[bestCol];
    int y = item.GetY();    
    int br = bestRow;
    int bc = bestCol;
    double localScore = DNA_Equal(colBases[m_x_off + bestCol], rowBases[y + bestRow]);
    score += localScore;
    sum++;
    bestCol = item[br].GetBackCol();  
    bestRow = item[br].GetBackRow();
    double dist = (bc - bestCol);
    dist += (br - bestRow);
    dist /= 2;
    
    sum += dist;
  }
  if (sum < 1.)
    return 0.;
  double ident = score / sum;
  double s = GetRawScore(ident, (int)(sum + 0.5), (int)(sum + 0.5), 0.5, 0.5);
  
  //cout << "Traceback score - ident=" << ident << " len=" << sum << " score=" <<  s << endl;
  return 2 * s;
}

void DPAligner::TraceBack(FullAlignment & align)
{
  int i, j;
 
  const DPAlignNodeColumn & last = m_columns[m_columns.isize()-1];
 
  double bestScore = 999999999999.;
  int bestDir = 0;
  int bestIdx = -1;
  for (i=0; i<last.isize(); i++) {
    //==================================================================
    if (i < last.isize() - 1)
      continue;
    //==================================================================

    //cout << "Check # " << i << " Score=" << last[i].GetScore() << endl;
    if (last[i].GetScore() < bestScore) {
      bestScore = last[i].GetScore();
      bestIdx = i;
      bestDir = last[i].IsBestRC();
    }
  }
  //cout << "TRACEBACK: Best index=" << bestIdx << " score=" << bestScore << endl;

  int bestRow = bestIdx;
  int bestCol = m_columns.isize()-1;

  svec<TraceBackItem> trace;

  int len = 0;

  while (bestRow != -1) {
 
    const DPAlignNodeColumn & item = m_columns[bestCol];
    int br = bestRow;

    //cout << "br=" << br << endl;
    bestScore = item[br].GetScore();
    int q = 0;

    //cout << "Pushing " << bestRow << "\t" <<  bestCol  << "\t" <<  item[br].IsBestRC() << "\t" << bestScore << endl;
    trace.push_back(TraceBackItem(bestRow, bestCol, bestDir, bestScore, q));
    bestCol = item[br].GetBackCol();
    //bestDir = item[br].IsRC();
    bestRow = item[br].GetBackRow();
    
    
    if (bestRow >= br) {
      cout << "ERROR: br=" << br << " last br=" << bestRow << endl;
    } 
    //cout << "Win score: " << bestScore << " frameFW: " << m_bestFrameScoreFW[br] << " frameRC: " << m_bestFrameScoreRC[br] << endl;
  }
  //cout << "# of items: " << trace.isize() << endl;

  //int lastRow = m_columns[0].GetY();
  //int lastCol = m_x_off;
  int lastDir = 0;

  //int indexRow = m_columns[0].GetY();
  int lastCol = m_x_off;
  int lastRow = m_columns[0].GetY();

  for (i=trace.isize()-1; i>=0; i--) {
    int row = trace[i].Row();
    int col = trace[i].Col();
    int dir = trace[i].Dir();
    int q = trace[i].Quest();
    const DPAlignNodeColumn & item = m_columns[col];

    row += item.GetY();
    col += m_x_off;

    //cout << " c=" << col << " r=" << row << " dir=" << dir << endl;


    // Vertical indel
    for (j=lastCol+1; j<col; j++) {
      align.Push(j, -1, dir, q);
      //cout << "Pushing v indel (col)" << endl;
    }
    
    // Horizontal indel
    if (lastDir == dir) {
      if (dir == 1) {
	for (j=lastRow+1; j<row; j++) {
	  align.Push(-1, j, dir, q);
	  //cout << "Pushing h indel (row fw)" << endl;
	}
      } else {

	for (j=row-1; j>lastRow; j--) {
	  align.Push(-1, j, dir, q);
	  //cout << "Pushing h indel (row rc)" << endl;
	}

      }
    }
    align.Push(col, row, dir, q);

    lastCol = col;
    lastRow = row;
    lastDir = dir;

  }


  //cout << "Align length: " << align.Size() << endl;
  
}

