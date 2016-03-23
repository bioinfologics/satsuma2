#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "analysis/GridSearch.h"



inline long long big_random( )
{    
  return ( (long long)( rand( ) ) << 31 ) ^ (long long)( rand( ) );    
}

GridSearch::GridSearch(int size, int blocks)
{
  m_size = size;
  
  //m_blocks = 24;
  m_blocks = blocks;
}



void GridSearch::SetUp(const vecDNAVector & target, 
		       const vecDNAVector & query)
{
  int i, j;

  ChunkManager tChunk(m_size, m_size/4);
  ChunkManager qChunk(m_size, 0);
  m_targetSeq.resize(target.isize());
  m_querySeq.resize(query.isize());


  // Target first
  svec<string> names;
  names.resize(target.isize());
  for (i=0; i<target.isize(); i++) {
    const char * p = target.Name(i).c_str();
    names[i] = &p[1];
    m_targetSeq[i].SetName(names[i], target[i].isize());
  }

  vecDNAVector tmp;
  tChunk.ChunkIt(tmp, m_targetChunks, target, names);
  m_allNs.resize(tmp.isize(), false);
  // Remember where the bad chunks are!!

  for (i=0; i<tmp.isize(); i++) {
    if (tmp[i].isize() == 0) {
      //cout << "Dismissing, all N's" << endl;
      m_allNs[i] = true;
    }
  }


  // Query
  names.resize(query.isize());
  for (i=0; i<query.isize(); i++) {
    const char * p = query.Name(i).c_str();
    names[i] = &p[1];
    m_querySeq[i].SetName(names[i], query[i].isize());
  }

  qChunk.ChunkIt(tmp, m_queryChunks, query, names);

  // Find first and last blocks...

  j = 0;
  int first = -1;
  int last = 0;
  for (i=0; i<m_targetSeq.isize(); i++) {
    const string & n = m_targetSeq[i].Name();
    first = -1;
    last = 0;
    for (j=0; j<m_targetChunks.isize(); j++) {
      if (m_targetChunks[j].GetName() == n) {
      //	m_targetSeq[i].SetBlocks(first, j-1);
      //	break;
      //}
	if (first == -1)
	  first = j;
	last = j;
      }
    }
    m_targetSeq[i].SetBlocks(first, last);
    //j--;
  }
   
  j = 0;
  first = -1;
  last = 0;
  for (i=0; i<m_querySeq.isize(); i++) {
    const string & n = m_querySeq[i].Name();
    first = -1;
    last = 0;
    for (j=0; j<m_queryChunks.isize(); j++) {
      if (m_queryChunks[j].GetName() == n) {
	if (first == -1)
	  first = j;
	last = j;
      }
    }
    m_querySeq[i].SetBlocks(first, last);

    //    for (; j<m_queryChunks.isize(); j++) {
    //if (m_queryChunks[j].GetName() != n || j+1 == m_queryChunks.isize()) {
    //m_querySeq[i].SetBlocks(first, j-1);
    //break;
    //}
    //if (first == -1)
    //first = j;
    //}
    //j--;
  }


  m_matrix.Setup(m_targetChunks.isize(), m_queryChunks.isize(), m_blocks, m_blocks);

  
  m_vertical.resize(m_queryChunks.isize(), 0);
  m_horizontal.resize(m_targetChunks.isize(), 0);

  m_repTrack.Setup(target.isize(), query.isize());

  for (i=0; i<target.isize(); i++)
    m_repTrack.SetTargetSize(i, target[i].isize());
  for (i=0; i<query.isize(); i++)
    m_repTrack.SetQuerySize(i, query[i].isize());

}

void GridSearch::SetUsed(int targetID,
			 int startTarget,
			 int endTarget,
			 int queryID,
			 int startQuery,
			 int endQuery)
{
  int targetBlock, queryBlock;
  int x, y;
  CoordsToBlocks(targetBlock,
		 queryBlock,
		 targetID,
		 startTarget,
		 queryID,
		 startQuery);

  m_matrix.BlockByAbsolute(x, y, targetBlock, queryBlock);

  int i, j;
  for (i=x-1; i<=x+1; i++) {
    for (j=y-1; j<=y+1; j++) {
      if (i < 0 || j < 0)
	continue;
      if (i >= m_matrix.TargetBlocks() || j >= m_matrix.QueryBlocks())
	continue;
      
      m_matrix.Set(i, j, SEARCH_DONE);
    }
  }
}

void GridSearch::ConsiderTargets(int targetID,
				 int startTarget,
				 int endTarget,
				 int queryID,
				 int startQuery,
				 int endQuery,
				 double ident)
{
  int targetBlock, queryBlock;


  if (m_repTrack.IsRepeat(targetID, startTarget, endTarget, queryID, startQuery, endQuery)) {
    //cout << "Repeat????" << endl;
    //return;
  }

  // This is very stupid...
  m_repTrack.SetRepeat(targetID, startTarget, endTarget, queryID, startQuery, endQuery);


  CoordsToBlocks(targetBlock,
		 queryBlock,
		 targetID,
		 (startTarget + endTarget)/2,
		 queryID,
		 (startQuery + endQuery)/2);

  //cout << "Mapped to blocks t=" << targetBlock << " q=" << queryBlock << endl;

  int x, y;
  //m_matrix.ClearCounts();

  m_matrix.BlockByAbsolute(x, y, targetBlock, queryBlock);

  //cout << "Grid coords x=" << x << " y=" << y << " (size_x=" << m_matrix.TargetBlocks();
  //cout << " size_y=" << m_matrix.QueryBlocks() << ")" << endl;
  
  //m_matrix.Set(x, y, SEARCH_HIT);


  const SearchCoordinates & c = m_matrix.Coordinates(x, y);
  //cout << "Translate back to blocks: " << c.TargetFirst() << "\t" << c.TargetLast();
  //cout << "\t" << c.QueryFirst() << "\t" << c.QueryLast() << endl;


  m_vertical[queryBlock] += endTarget - startTarget;
  m_horizontal[targetBlock] += endTarget - startTarget;
     
  int i, j;
  int plus = 0;
  if (m_matrix.Get(x, y) != SEARCH_UNKNOWN)
    plus = 1;

  for (i=x-plus; i<=x+plus; i++) {
    for (j=y-plus; j<=y+plus; j++) {
      if (i < 0 || j < 0)
	continue;
      if (i >= m_matrix.TargetBlocks() || j >= m_matrix.QueryBlocks())
	continue;
      
      //cout << "Adding to element " << i << " " << j << ", is " << m_matrix.Print(i, j) << endl;
      int w = (int)(100. * ident * (double)(endTarget - startTarget));
      m_matrix.AddCount(i, j, w);

 
    }
  }
}


void GridSearch::CollectSeeds(svec<GridTarget> & targetsOut, int n)
{
  if (n == 0)
    return;

  cout << "Looking for additional seeds" << endl;

  svec<GridTarget> targets;

  int i;

  double sum = 0.;
  double num = 0.;
  for (i=0; i<m_vertical.isize(); i++) {
    if (m_vertical[i] > 0) {
      sum += (double)m_vertical[i];
      num += 1.;
    }
  }
  sum /= (num + 1.);
  cout << "Average vertical occupancy: " << sum << endl;


  int skipMax = 16;

  //cout << "Vertical size: " << m_vertical.isize() << endl;

  int j;

  int lineCount = 0;
  int lineCountBig = 0;

  for (i=0; i<m_horizontal.isize(); i++) {
    if (m_horizontal[i] > 0)
      continue;

    if (m_allNs[i])
      continue;
      

    for (j=i; j<m_horizontal.isize(); j++) {
      if (m_horizontal[j] > 0 || m_allNs[j])
	break;
    }

 
    int len = j-i-1;

    if (len < 4)
      continue;

    lineCount++;
    

    if (len >= 4)
      lineCountBig++;

    int x = i + len/2;

    i = j;

    for (int y=0; y<m_vertical.isize(); y++) {
      if (m_vertical[y] > 0)
	continue;

      int skip = 0;

      for (j=y; j<m_vertical.isize(); j++) {      
	if (m_vertical[j] > 0) {
	  skip++;
	  if (skip >= skipMax)
	    break;
	} else {
	  skip = 0;
	}
      }
      int y1 = y; 
      int y2 = j-1;

      if (y2-y1 < 2)
	continue;

      
      targets.push_back(GridTarget(x,
				   x,
				   y1,
				   y2,
				   x, 
				   y1, 
				   y2, 
				   len,
				   true));

      y = j;
    }
  }



  //cout << "Potential targets: " << targets.isize() << endl;

  Sort(targets);

  cout << "Grid search: potential lines to be farmed out: " << lineCount << " big ones: "<< lineCountBig << endl;
 
  if (n < targets.isize()) {
    int lastTarget = targets[n-1].TargetFirst();
    for (i=n; i<targets.isize(); i++) {
      if (targets[i].TargetFirst() != lastTarget) {
	break;
      }
    }
    targets.resize(i);
    cout << "Requested: " << n << " retrieved: " << i << endl;
  }

  for (i=0; i<targets.isize(); i++) {
    for (j=targets[i].TargetFirst(); j<=targets[i].TargetLast(); j++)
      m_horizontal[j] = 1;
    targetsOut.push_back(targets[i]);
  }


}




int GridSearch::CollectTargets(svec<GridTarget> & targets, int n, int nSeeds)
{
  int i, j;

  for (i=0; i<m_matrix.TargetBlocks(); i++) {
    for (j=0; j<m_matrix.QueryBlocks(); j++) {
      if (m_matrix.Get(i, j) != SEARCH_UNKNOWN)
	continue;
      

      if (m_matrix.GetCount(i, j) == 0)
	continue;

      
      const SearchCoordinates & c = m_matrix.Coordinates(i, j);
      
      targets.push_back(GridTarget(c.TargetFirst(),
				   c.TargetLast(),
				   c.QueryFirst(),
				   c.QueryLast(),
				   i, j,  m_matrix.GetCount(i, j)));
    }
  }


  Sort(targets);

  cout << "Potential targets: " << targets.isize() << endl;

  int pTar = targets.isize();

  if (n < targets.isize()) {
    cout << "Downsize to " << n << endl;
    targets.resize(n-nSeeds);
  }

  for (i=0; i<targets.isize(); i++) {
    m_matrix.Set(targets[i].X(), targets[i].Y(), SEARCH_SUBMITTED);
    if (targets[i].IsFast())
      cout << "ERROR: Target is FAST!" << endl;

    cout << "Mark used x=" << targets[i].X() << " y=" <<  targets[i].Y() << " count=" << targets[i].GetCount();
    cout << " [" << targets[i].TargetFirst() << "-" << targets[i].TargetLast() << ",";
    cout << targets[i].QueryFirst() << "-" << targets[i].QueryLast() << "]" << endl;
  }

  if (targets.isize() < n) {
    cout << "Adding in SEED targets..." << endl;
    CollectSeeds(targets, n-targets.isize());
  }

  //for (i=0; i<targets.isize(); i++) {
  //  if (targets[i].IsFast())
  //    cout << "ERROR(2): Target is FAST!" << endl;
  //}

  return pTar;
}


void GridSearch::CoordsToBlocks(int & targetBlock,
				int & queryBlock,
				int target,
				int startTarget,
				int query,
				int startQuery)
{
  //cout << "Mapping t=" << target << " q=" << query << endl;

  int bStartTarget = m_targetSeq[target].First();
  int bStartQuery = m_querySeq[query].First();

  targetBlock = queryBlock = -1;


  int off = startTarget/(m_size-m_size/4);
  targetBlock = m_targetSeq[target].First() + startTarget/(m_size-m_size/4);
  if (off > 0 && startTarget - (m_size-m_size/4)*off < m_size/4)
    targetBlock--;

  queryBlock = m_querySeq[query].First() + startQuery/m_size;

}

bool GridSearch::GetSelect(SeqChunkSelect & vertical, svec<int> & horizontal, int howMany)
{
  int i, j;
 
  vertical.Set(m_vertical.isize(), m_horizontal.isize());

  double sum = 0.;
  double n = 0.;
  for (i=0; i<m_vertical.isize(); i++) {
    if (m_vertical[i] > 0) {
      sum += (double)m_vertical[i];
      n += 1.;
    }
  }
  sum /= (n + 1.);
  cout << "Average vertical occupancy: " << sum << endl;

  double thresh = sum * 0.25;

  int count = 0;
  for (i=0; i<m_vertical.isize(); i++) {
    if ((double)m_vertical[i] < thresh) {
      vertical.SetQuery(i, 1);
      count++;
    }
  }

  cout << count << " chunks out of " << m_vertical.isize() << " are free." << endl;

  count = 0;
  for (i=0; i<m_horizontal.isize(); i++) {
    if (m_horizontal[i] == 0) 
      count++;
  }

  double prob = (double)howMany / (double)count;

  count = 0;
  for (i=0; i<m_horizontal.isize(); i++) {
    if (m_horizontal[i] == 0) {
      long long l = big_random( ) % 10000;
      double p = (double)l / 10000;

      if (p < prob) {
	horizontal.push_back(i);
	m_horizontal[i] = 1;
      }
    }
  }

  cout << "Selected " << horizontal.isize() << " horizontal regions" << endl;

  if (horizontal.isize() > 0)
    return true;
  else
    return false;
}
