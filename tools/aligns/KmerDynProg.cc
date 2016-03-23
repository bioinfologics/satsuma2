

#ifndef FORCE_DEBUG

#define NDEBUG

#endif

#include "aligns/KmerDynProg.h"


#define SILENT

int GetSepDevScore(int sep, int dev) 
{
  //cout << "sep=" << sep << " dev=" << dev << endl;
  if (dev == 0)
    return sep;
  
  if (sep < 0)
    sep = -sep;
  
  sep -= 4*dev;
  if (sep < 0)
    sep = 0;

  if (sep > 128)
    sep = 128;

  //cout << "SepDevScore=" << sep << endl;


  return sep;
}


string Date() {
  static string dummy;
  return dummy;
}

int MatchScore::GetScore(const KSAMatch & from,
			 const KSAMatch & to,
			 int sep,
			 int dev)
{

  if (to.GetRefID() == NULL_KMER_ID && from.GetRefID() == NULL_KMER_ID)
    return 2 * m_nullScore;

  if (to.GetRefID() == NULL_KMER_ID) {
    return m_globScore;
  }

  int fromRef = from.GetRefID();

  bool bRev = false;
  if (from.GetRefID() == NULL_KMER_ID) {
    fromRef = from.GetFeed();

    if (fromRef == NULL_KMER_ID) {
      return m_nullScore;
    }
    return 1;
  }
  

  if (from.IsRC() != to.IsRC())
    return  m_globScore;

  if (fromRef != to.GetRefID())
    return m_globScore;

  int score = 0;
  if (from.IsRC()) {
    score = from.GetRefStart() - to.GetRefStart();// - sep;
  } else {
    score = to.GetRefStart() - from.GetRefStart();// - sep;
  }
  if (score < 0)
    score = -score;

  score -= sep;

  score -= 4*dev;

  if (score < 0)
    score = 0;
  

  if (score > m_globScore)
    score = m_globScore;


#ifndef SILENT
  if (bRev) {
    *m_log << "Score: " << score << " from=" << from.GetRefID() << " to " << to.GetRefID() << " sep=" << sep << endl;
  }
#endif
  //if (score < 0)
  //score = 0;
  return score;
}



//#define PRINT_DETAILS





KmerSuperAligner::KmerSuperAligner( ostream *log )
{
  m_lastBestScore = 0;
  m_upCount = 0;
  m_kmerShift = 1;

  m_maxFreq = -1;

  m_matchScore.SetLog( log );

  m_log = log ? log : (ostream *)&cout;

  m_bUseCore = true;

  m_numKmers = 2;
  
  m_bUseCombs = false;
  m_bUseProteins = false;
  
  m_lookAheadFrames = 16;
  m_lookAheadStep = 2;
  m_lookAheadMaxFreq = 1000;
}

KmerSuperAligner::~KmerSuperAligner() 
{
}

void KmerSuperAligner::SetRefBases(const vecDNAVector & b)
{
  int i;
  m_core.SetNumTables(m_numKmers);
  if (m_bUseCombs) {
    m_core.SetTranslator(&m_combTrans);
  } else {
    if (m_bUseProteins) {
      m_core.SetTranslator(&m_proteinTrans);
    } else {
      m_core.SetTranslator(&m_numberTrans);
    }
  }
  
  
  cout << "Sequences: " << b.size() << endl;
  //---------------------------
  m_core.AddData(b, false, 1);
  //---------------------------
  m_core.SortAll();
}

void KmerSuperAligner::SetRefBases(const vecDNAVector & b, const vecNumVector & tags, int min)
{
  int i;
  m_core.SetNumTables(m_numKmers);
  if (m_bUseCombs) {
    m_core.SetTranslator(&m_combTrans);
  } else {
    if (m_bUseProteins) {
      m_core.SetTranslator(&m_proteinTrans);
    } else {
      m_core.SetTranslator(&m_numberTrans);
    }
  }
  
  
  cout << "Sequences: " << b.size() << endl;
  cout << "Using repeat mask!" << endl;
  //---------------------------
  m_core.AddData(b, tags, min, false, 1);
  //---------------------------
  m_core.SortAll();
 }




int KmerSuperAligner::MatchCloseExtensions(int fromIndex, int toIndex, int addScore)
{
  KSAItem & from = m_dynArray[fromIndex];
  KSAItem & to = m_dynArray[toIndex];
  int i, j, k;

  int dev = to.GetDev();

  int globDist = m_matchScore.GetGlobalDistanceLimit();
  
  int bestIndex = -1;
  int bestScore = 0x7FFFFFFF;

  j = 0;

  int sep = to.GetFromOnSuper() - from.GetFromOnSuper();


  int skipPenalty = 0;
  if (toIndex - fromIndex > 1) {
   int wordSize = KSA_K;
   
   if (m_bUseCore) 
     wordSize = m_core.GetWordSize();
  
   skipPenalty = toIndex - fromIndex - wordSize;
   if (skipPenalty < 1)
     skipPenalty = 1;
  }

  //cout << "From index=" << fromIndex << " toIndex=" << toIndex << endl;

  for (i=0; i<from.GetCount(); i++) {
    const KSAMatch & f = from.Get(i);

    //cout << "i=" << i << " score=" << f.GetScore() << " id=" << f.GetRefID() << " refStart=" << f.GetRefStart() << endl;

    if (f.GetScore() + addScore < bestScore) {
      bestIndex = i;
      bestScore = f.GetScore() + addScore;
    }
    //*m_log << "ff" << endl;
    // Let's fast forward here...
    for (; j<to.GetCount(); j++) {
      const KSAMatch & t = to.Get(j);
      if (t < f) {
	continue;
      }
      break;
    }
    //*m_log << "j=" << j << endl;
    
    k = j-1;

    //*m_log << "forward" << endl;
    // Run forward as much as we have to 
    for (; j<to.GetCount(); j++) {
      const KSAMatch & t = to.Get(j);

  
      int score = m_matchScore.GetScore(f, t, sep, dev) + addScore;
      if (score >= globDist)
	break;
      int transScore = score + f.GetScore() + skipPenalty;
      to.Update(j, transScore, i, fromIndex);
      m_upCount++;
   }
    //*m_log << "j=" << j << endl;
    j = k;

    //*m_log << "back" << endl;
    //...and backwards, and leave the index j right where it ends up...
    for (; j>=0; j--) {
      const KSAMatch & t = to.Get(j);
      int score = m_matchScore.GetScore(f, t, sep, dev) + addScore;
      if (score >= globDist)
	break;
      int transScore = score + f.GetScore() + skipPenalty;
      to.Update(j, transScore, i, fromIndex);      
      m_upCount++;
    }
    //*m_log << "j=" << j << endl;
    if (j < 0)
      j = 0;
  }

  return bestIndex;
}

int KmerSuperAligner::ExtendAll(int fromIndex, int bestIndex, int toIndex)
{
  KSAItem & from = m_dynArray[fromIndex];
  KSAItem & to = m_dynArray[toIndex];
  int j;

  int bestLastIndex = -1;
  const KSAMatch & f = from.Get(bestIndex);
  m_lastBestScore = 0x7FFFFFFF;
  int sep = to.GetFromOnSuper() - from.GetFromOnSuper();

  for (j=0; j<to.GetCount(); j++) {
    const KSAMatch & t = to.Get(j);
    int dev = to.GetDev();
    int score = m_matchScore.GetScore(f, t, sep, dev);
    int transScore = score + f.GetScore();
    bool b = to.Update(j, transScore, bestIndex, fromIndex);
    m_upCount++;

    // NULL kmer??
    if (b && t.GetRefID() == NULL_KMER_ID) {
      KSAMatch & two = to.GetIt(j);

      if (f.GetRefID() != NULL_KMER_ID)
	two.SetFeed(f.GetRefID());
      else
	two.SetFeed(f.GetFeed());
      two.Set(NULL_KMER_ID, 
	      f.GetRefStart() + 1, 
	      f.GetRefStop() + 1,
	      f.IsRC());
	  
    }

    if (f.GetScore() < m_lastBestScore) {
      m_lastBestScore = f.GetScore();
      bestLastIndex = j;
    }
    
  }
  return bestLastIndex;
}

void KmerSuperAligner::Align(SuperAlign & result, 
			     const vecDNAVector & contigBases,
			     svec<int> & contigStartInSuper,
			     svec<int> & gapDevInSuper,
			     int skipContig)
{
  int i, j, k;

  //cout << "Skip " << skipContig << endl;

  int totalCount = 0;
  int wordSize = KSA_K;

  if (m_bUseCore) {
    wordSize = m_core.GetWordSize();
    //cout << "wordSize=" << wordSize << endl;
  }


  for (i=0; i<(int)contigBases.size(); i++) {
    for (j=0; j<(int)contigBases[i].size()-wordSize; j+=m_kmerShift) {
      totalCount++;
    }
    totalCount++;
    //totalCount += (int)contigBases[i].size() - KSA_K;
  }

  //cout << "total bases: " << totalCount << ", contigs: " << contigBases.size() << endl;

  m_dynArray.clear();
  m_dynArray.resize(totalCount);
  
  /*
  k = 0;
  for (i=0; i<(int)contigBases.size(); i++) {
    int offset = contigStartInSuper[i];
    
    for (j=0; j<(int)contigBases[i].size()-wordSize; j+=m_kmerShift) {
      int dev = 0;
      if (i > 0 && j == 0)
	dev = gapDevInSuper[i-1];

      m_dynArray[k].Set(offset + j, j, i, dev);
      // Add the "null" kmer!!
      //m_dynArray[k].AddMatch(KSAMatch());
      k++;
    }
    }*/

  // Now, let's fill out the dynamic programming structures...


  DNAVector kmerbases;
  k = 0;

#ifndef SILENT
  *m_log << Date() << " Building dyn prog structure..." << endl;
#endif

  int plusminus = m_core.GetLookAhead();
 

  for (i=0; i<(int)contigBases.size(); i++) {
    const DNAVector & b = contigBases[i];
    int offset = contigStartInSuper[i];
    //cout << (int)b.size()-wordSize << endl;

    for (j=0; j<=(int)b.size()-wordSize; j+=m_kmerShift) {      
      if (j + wordSize > b.isize())
	break;
      //cout << "Adding j=" << j << endl;
      int dev = 0;
      if (i > 0 && j == 0)
	dev = gapDevInSuper[i-1];

      m_dynArray[k].Set(offset + j, j, i, dev);


      kmerbases.SetToSubOf(b, j, wordSize);

      //cout << b.isize() << " - " << kmerbases.isize() << endl;

      KSAItem & dyn = m_dynArray[k];


      svec<KmerAlignCoreRecordWithScore> matches;
      svec<KmerAlignCoreRecordWithScore> rc_matches;
      m_core.GetMatches(matches, kmerbases, 0);
      
      kmerbases.ReverseComplement();
      m_core.GetMatches(rc_matches, kmerbases, 0); 
 
      //cout << "RC Matches: " << rc_matches.isize() <<  " FW Matches: " << matches.isize() << endl;

      int skipSelfCount = 0;
      if (skipContig != -1) {
	if (m_maxFreq == -1 || matches.isize() <= m_maxFreq) {
	  for (int x=0; x<matches.isize(); x++) {
	    if (matches[x].GetContig() == skipContig && j == matches[x].GetPosition()) {
	      skipSelfCount++;
	    }
	  }
	}
      }

      

      int count = dyn.GetSize();
      int startCount = count;
      if (m_maxFreq == -1 || matches.isize() <= m_maxFreq)
	count += matches.isize() - skipSelfCount;
      
      
      if (m_maxFreq == -1 || rc_matches.isize() <= m_maxFreq) 
	count += rc_matches.isize();
      
      if (count == 0)
	continue;

      int checkCount = count;
      dyn.SetSize(count);
      
      count = startCount;
      
      if (m_maxFreq == -1 || matches.isize() <= m_maxFreq) {
	for (int x=0; x<matches.isize(); x++) {
	  if (matches[x].GetContig() == skipContig && j == matches[x].GetPosition()) {
	    //cout << "Skipping matches from contig " << skipContig << endl;
	    continue;
	  }
	  dyn.AddMatch(matches[x].GetContig(), matches[x].GetPosition(), matches[x].GetPosition() + wordSize - 1, 
		       false, count, matches[x].GetScore());
	  //dyn.AddMatch(matches[x].GetContig(), matches[x].GetPosition(), matches[x].GetPosition() + wordSize - 1, false);
	  count++;
	}
      }
      
      
      
      if (m_maxFreq == -1 || rc_matches.isize() <= m_maxFreq) {
	for (int x=0; x<rc_matches.isize(); x++) {
	  //cout << "Adding rc match from contig " << rc_matches[x].GetContig() << endl;
	  //if (rc_matches[x].GetContig() == skipContig)
	  //continue;
	  dyn.AddMatch(rc_matches[x].GetContig(), rc_matches[x].GetPosition(), rc_matches[x].GetPosition() + wordSize - 1, 
		       true, count, rc_matches[x].GetScore());
	  //dyn.AddMatch(rc_matches[x].GetContig(), rc_matches[x].GetPosition(), rc_matches[x].GetPosition() + wordSize - 1, true);
	  count++;
	}
      }
      
      if (count != checkCount)
	dyn.SetSize(count);
      
      
      k++;
    }
  }


  //cout << "Downsize array from " << m_dynArray.isize() << " to " << k << endl;
  m_dynArray.resize(k);
  if (k == 0)
    return;


#ifndef SILENT
  *m_log << Date() << " done!" << endl;

  *m_log << Date() << " Dynamic programming..." << endl;
#endif
  int bestScore = 0x7FFFFFFF;
  int lastBestIndex = -1;


  KSAItem & first = m_dynArray[0];

  for (i=0; i<first.GetCount(); i++) {
    if (first.Get(i).GetRefID() == -1) {
      first.Update(i, 100000, -1, 0);
    } else {
      first.Update(i, 0, -1, 0);
    }
  }

#ifndef SILENT
  *m_log << Date() << " Sorting structures..." << endl;
#endif

  for (i=0; i<m_dynArray.isize(); i++) {
    m_dynArray[i].SortMatches();
  }

#ifndef SILENT
  *m_log << Date() << " done!" << endl;
#endif


  int lookAheadFrames = m_lookAheadFrames * m_lookAheadStep;

  for (i=1; i<m_dynArray.isize(); i++) {

    //int posLast = m_dynArray[i-1].GetFromOnSuper();
    //int posNow = m_dynArray[i].GetFromOnSuper();
    //cout << "Transition to " << i << endl;
    lastBestIndex = MatchCloseExtensions(i-1, i, 0);
    
    for (j=i+wordSize+1; j<i+wordSize+1+lookAheadFrames; j+= m_lookAheadStep) {
      if (j >= m_dynArray.isize())
	break;
      int pen = j-i;
    
      int targetSize = m_dynArray[j].GetSize();
      if (targetSize < m_lookAheadMaxFreq)
	MatchCloseExtensions(i-1, j, pen);
    }

    lastBestIndex = ExtendAll(i-1, lastBestIndex, i);

    bestScore = m_lastBestScore;
    //if (i > 41212 && i < 41236)
    //cout << "i=" << i << " bestScore=" << bestScore << endl;
  }


#ifndef SILENT
  *m_log << Date() << " Best score=" << bestScore << endl;
  *m_log << Date() << " Traceback (backwards!!!)..." << endl;
  *m_log << Date() << " Hypotheses updated: " << m_upCount << endl;
#endif

  SuperAlign tmpMatches;

  int theBestScore = m_dynArray[m_dynArray.isize()-1].Get(0).GetScore();

  KSAMatch lastBest;

  //for (i=m_dynArray.isize()-1; i>=0; i--) {
  //  if (m_dynArray[i].GetSize() > 1)
  //    break;
  //}
  int frameIndex = m_dynArray.isize()-1;


  const KSAItem & lastFrame = m_dynArray[frameIndex];
  for (i=0; i<lastFrame.GetSize(); i++) {
    const KSAMatch & cand  =  m_dynArray[frameIndex].Get(i); 
    if (cand.GetScore() <= theBestScore) {
      theBestScore = cand.GetScore();
      lastBestIndex = i;
      lastBest = cand;
    }
  }

  //cout << "lastBestIndex= " << lastBestIndex << " score=" << lastBest.GetScore() << endl;





  /*
  for (i=0; i<m_dynArray.isize(); i++) {
    const KSAItem & frame = m_dynArray[i];
    cout << i <<  endl;
    for (int y=0; y<frame.GetSize(); y++) {
      cout << "  #=" << y << " score=" << frame.Get(y).GetScore() << endl;
    }
    }*/
  

  while (frameIndex >= 0 && lastBestIndex != -1) {
    const KSAItem & frame = m_dynArray[frameIndex];
    const KSAMatch & best = frame.Get(lastBestIndex);
    
 
    //cout << frameIndex << "  " << lastBestIndex <<  " score=" << best.GetScore() << endl;

    int startContig = frame.GetContigIndex();
    
    int startOnContig = frame.GetFromOnContig();
    int endOnContig = startOnContig + wordSize - 1;

    int startOnSuper = frame.GetFromOnSuper();
    int startOnRef = best.GetRefStart();
    int endOnRef = startOnRef + wordSize - 1;

    int len = endOnContig - startOnContig;
    SuperMatch match;
    match.Set(startContig,
	      startOnContig,
	      endOnContig,
	      best.GetRefID(),
	      best.IsRC(),
	      startOnRef,
	      endOnRef);

    match.SetScore(0);

    tmpMatches.AddMatch(match);

    lastBestIndex = best.GetPrevIndex();
    frameIndex = best.GetPrevFrame();
    
  }

  //cout << "done!" << endl;

  for (i=tmpMatches.GetMatches().isize()-1; i>=0; i--) {
    const SuperMatch & match = tmpMatches.GetMatch(i);
    for (j=i-1; j>=0; j--) {
      const SuperMatch & match2 = tmpMatches.GetMatch(j+1);
      const SuperMatch & match3 = tmpMatches.GetMatch(j);
      if (!match2.IsContiguous(match3))
	break;
    }
    i = j+1;

    const SuperMatch & match4 = tmpMatches.GetMatch(i);

    
    SuperMatch theMatch;

    int startOnRef = match.GetRefStart();
    int endOnRef = match4.GetRefEnd();

    if (match.GetRC()) {
      startOnRef = match4.GetRefStart();
      endOnRef = match.GetRefEnd();
    }

    theMatch.Set(match.GetContig(),
		 match.GetFirstBase(),
		 match4.GetLastBase(),
		 match.GetRefID(),
		 match.GetRC(),
		 startOnRef,
		 endOnRef);
    if (match.GetRefID() != -1)
      result.AddMatch(theMatch);
  }


  //cout << "Leaving dyn prog." << endl;

      //--------------------------------------------------------
#if 0 

  //for (i=m_dynArray.isize()-2; i>=0; i--) {
  i = lastBest.GetPrevFrame();
  int shift = m_dynArray.isize() - 1 - i;

  bool bFirst = true;
  //cout << "Initial shift=" << shift << endl;
  while (i >= 0 && lastBestIndex != -1) {
    const KSAItem & to = m_dynArray[i];
    const KSAMatch & best = to.Get(lastBestIndex);

    //cout << "i=" << i << " lastBestIndex=" << lastBestIndex << " prevFrame=" << best.GetPrevFrame() << endl;

    //*m_log << "i=" << i << endl;
    //*m_log << "best id=" << best.GetRefID() << " index=" << lastBestIndex << " out of " << to.GetCount() << endl;
    
    //if (lastBest == -1)
    //cout << "OUCH!" << endl;

    if (i>0 && best.IsContiguous(lastBest, shift) 
	&& to.GetContigIndex () == m_dynArray[i+1].GetContigIndex ()) {
      lastBest = best;
      lastBestIndex = best.GetPrevIndex();
      shift = i - best.GetPrevFrame();
      i = best.GetPrevFrame();
      continue;
    }

    //*m_log << "Adding block!" << endl;
    //*m_log << "last id=" << lastBest.GetRefID() << endl;
 

    //const KSAItem & lastGood = m_dynArray[i+1];
    const KSAItem & lastGood = m_dynArray[i];
    const KSAItem & firstGood = m_dynArray[firstMatchIndex];
    int startContig = lastGood.GetContigIndex();
    
    int startOnContig = lastGood.GetFromOnContig();
    int startOnSuper = lastGood.GetFromOnSuper();

    //int startOnRef = startOnSuper;
    //int startOnRef = lastBest.GetRefStart();

    int startOnRef = best.GetRefStart();
    if (bFirst) {
      startOnRef = lastBest.GetRefStart();
      bFirst = false;
    }

    int endOnContig = firstGood.GetFromOnContig() + wordSize - 1;
    int len = endOnContig - startOnContig;

    //    cout << "startOnContig=" << startOnContig << " startOnSuper=" << startOnSuper << " startOnRef=" << startOnRef << endl;

    if (lastBest.GetRefID() == NULL_KMER_ID) {
      len -= 2 * wordSize;
    }
    //if (best.GetRefID() == NULL_KMER_ID) {
    //  len -= 2 * wordSize;
    //}
    int endOnSuper = startOnSuper + len;
    int endOnRef = startOnRef + len;
    if (lastBest.IsRC()) {
      endOnRef = lastBest.GetRefStart() + wordSize - 1;
      startOnRef = endOnRef - len;
    }

    if (lastBest.GetRefID() != NULL_KMER_ID) {

#ifdef PRINT_DETAILS
      *m_log << "c" << startContig << "[ " << startOnContig << " - " << endOnContig << " ]";
      *m_log << " w/ " << lastBest.GetRefID();
      if (lastBest.GetRefID() != NULL_KMER_ID)
	*m_log << " [ " << startOnRef << " - " << endOnRef << " ]";
      if (lastBest.IsRC())
	*m_log << " rc";
      else
	*m_log << " fw";
      
      *m_log << " len=" << len;
      *m_log << " SCORE=" << lastBest.GetScore() << endl;
#endif
      //*m_log << firstGood.GetRefID() << endl;


      SuperMatch match;
      match.Set(startContig,
		startOnContig,
		endOnContig,
		lastBest.GetRefID(),
		lastBest.IsRC(),
		startOnRef,
		endOnRef);

      match.SetScore(lastBest.GetScore());

      tmpMatches.AddMatch(match);

    }



    //*m_log << "Match contig " << to.GetContigIndex() << " w/ ref " << best.GetRefID() << " @ " << best.GetRefStart();
    //if (best.IsRC())
    //  *m_log << "rc";
    //else
    //  *m_log << "fw";
    //*m_log << endl;
    
    firstMatchIndex = i; 
    lastBest = best;
    lastBestIndex = best.GetPrevIndex();
    shift = i - best.GetPrevFrame();

    if (i == 0)
      break;
    i = best.GetPrevFrame();
  }

  for (i=tmpMatches.GetMatches().isize()-1; i>=0; i--) {
    result.AddMatch(tmpMatches.GetMatches()[i]);
  }

#endif



#ifndef SILENT
  *m_log << Date() << " done!" << endl << endl;
#endif

}




 







