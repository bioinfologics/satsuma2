#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "KmerAlignCore.h"

template<class DataType>
KmerAlignCore<DataType>::KmerAlignCore() 
{
  m_numTables = 2;    
  m_pTrans = NULL;
  m_lookAhead = 0;
  m_lookAheadMaxFreq = 50000;

  m_max12 = 0x7FFFFFFF;
}

template<class DataType>
void KmerAlignCore<DataType>::AddData(const vecDNAVector & bases, bool filterLC, int stepSize)
{
  vecNumVector dummy;
  dummy.resize(bases.size());
  AddData(bases, dummy, 1, filterLC, stepSize);
}

template<class DataType>
void KmerAlignCore<DataType>::AddData(const vecDNAVector & bases,
   const vecNumVector & tags, int min, bool filterLC, int stepSize)
{
  // First, count k-mers
  int i, j, k;
  int size = m_pTrans->GetSize();

  svec<int> counts;
  counts.resize(m_pTrans->GetBoundValue(), 0);
  svec< svec<int> > baseConversions;
  baseConversions.resize(bases.isize());

  cout << "Counting k-mers..." << endl;
  for (j=0; j<(int)bases.size(); j++) {
    const DNAVector & b = bases[j];
    const NumVector & t = tags[j];

    if (j % 10000 == 0)
      cout << flush << "\rContigs: " << j;
    baseConversions[j].resize(b.size()-size+1, -1);
    for(k=0; k <= (int)b.size()-size; k+=stepSize) {
      if (!IsRepeat(t, k, size, min) && (!filterLC || (filterLC  && !IsLowComplexity(b, k, size)))) {
	int n = m_pTrans->BasesToNumber(b, k);
	if (n >= 0) { counts[n]++; }
        baseConversions[j][k] = n;
      } 
    }
  }
  cout << "done, re-sizing arrays..." << endl;
  for (j=0; j<counts.isize(); j++) {
    m_table[j].Resize(counts[j]);
    counts[j] = 0;    
  }

  cout << "done, assigning k-mers..." << endl;
  for (j=0; j<(int)bases.size(); j++) {
    const DNAVector & b = bases[j];
    const NumVector & t = tags[j];

    if (j % 10000 == 0)
      cout << flush << "\rContigs: " << j;

    for(k=0; k <= (int)b.size()-size; k+=stepSize) {
      int n = baseConversions[j][k];
      if (n >= 0) {
        m_table[n].Add(j, k, counts[n]);
        counts[n]++;
      }
    }
  }
}

template<class DataType>
void KmerAlignCore<DataType>::AddData(const DNAVector & b, int contig, int offset, bool bSort)
{
  int i, j, k;

  int size = m_pTrans->GetSize();
  k = 0;

  while (k < (int)b.size()-size) {
    KmerAlignCoreRecordStoreTable<DataType> & t = m_table;
    int n = m_pTrans->BasesToNumber(b, k);
    if (n >= 0) {
      KmerAlignCoreRecordStore<DataType> & s   = t[n];
      s.Add(contig, k+offset);
    }
    k++;      
  }
  if (bSort)
    SortAll();
}

template<class DataType>
void KmerAlignCore<DataType>::SortAll()
{
  int i, j;
  KmerAlignCoreRecordStoreTable<DataType> & t = m_table;
  for (j=0; j<t.GetSize(); j++) {
    KmerAlignCoreRecordStore<DataType> & s    = t[j];
    s.Sort();    
  }
}

template<class DataType>
const svec<DataType>& KmerAlignCore<DataType>::GetMatchesDirectly(const DNAVector & b, int start)
{
  int size = m_pTrans->GetSize();
  int i, j;

  int n = m_pTrans->BasesToNumber(b, start);

  if (n < 0) {
    static svec<DataType> dummy;
    return dummy;
  }

  KmerAlignCoreRecordStore<DataType> & s = m_table[n];   
  return s.GetData();
}


template<class DataType>
bool KmerAlignCore<DataType>::GetMatches(svec<DataType>& matches, const DNAVector & b, int start)
{
  int size = m_pTrans->GetSize();
  int i, j, k, l;
  k = start;

  if (start + m_numTables * size > (int)b.size()) {
    cout << "Error: sequence length=" << b.size() << " and k-kmer end is " << start + m_numTables * size << endl; 
    return false;
  } 


  //cout << "Start position: " << start << endl;

  KmerAlignCoreRecordStoreTable<DataType> hits;

  hits.SetSize(m_numTables);

  for (i=0; i<m_numTables; i++) {
    KmerAlignCoreRecordStoreTable<DataType>& t      = m_table;
    KmerAlignCoreRecordStore<DataType> & r = hits[i];

    if (m_lookAhead == 0) {
      int n = m_pTrans->BasesToNumber(b, k + i * size);

      if (n >= 0) {
	KmerAlignCoreRecordStore<DataType> & s = t[n];   
	if (s.GetNumRecords() > m_max12) {
	  matches.clear();
	  return false;
	}
	
	// Pre-size it!
	r.Resize(s.GetNumRecords());
	int count = 0;
	//cout << "size=" << s.GetNumRecords() << endl;
	
	for (j=0; j<s.GetNumRecords(); j++) {
	  //cout << "Adding single match, pos=" << s.GetRecord(j).GetPosition() << endl;
	  r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size, count);
	  //cout << "Done" << endl;
	  count++;
	  //r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size);
	}
      }
    } else {
 
      int penalty = 0;
      for (l=0; l<m_lookAhead; l++) {
	if (k + i * size + l + size > (int)b.size())
	  break;
	//cout << "n=" << k + i * size + l << endl;
	int n = m_pTrans->BasesToNumber(b, k + i * size + l);

	if (n > 0) {
	  KmerAlignCoreRecordStore<DataType> & s = t[n];   
	  
	  int count = r.GetSize();
	  if (l > 0)
	    penalty = 1;
	  
	  if (l > 0 && count > m_lookAheadMaxFreq) {
	    //cout << "***** Count: " << count << endl;
	    continue;
	  }
	  r.Resize(count + s.GetNumRecords());
	  for (j=0; j<s.GetNumRecords(); j++) {
	    //cout << "Adding single match, pos=" << s.GetRecord(j).GetPosition() << endl;
	    r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size - l, count, penalty);
	    count++;
	    //r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size - l);
	  }
	}
	//cout << "Before sort: " << r.GetSize() << endl;
	r.UniqueSort();      
	//cout << "After sort:  " << r.GetSize() << endl;
	//for (l=0; l<r.GetSize(); l++) {
	//const KmerAlignCoreRecord & record = r.Get(l);
	//cout << "   c" << record.GetContig() << " @" << record.GetPosition() << endl;
	//}
      }
    }
  }

  matches.clear();
  svec<DataType> tmp;

  if (hits.GetSize() == 0)
    return false;

  //cout << "Merging, hits size=" << hits.GetSize() << endl;
  tmp = hits[0].GetData();
  for (i=1; i<m_numTables; i++) {
    MergeSortFilter(matches, tmp, hits[i].GetData());
    tmp = matches;
    //cout << "Matches remaining: " << matches.isize() << endl;
  }

  if (m_numTables == 1)
    matches = tmp;
  //cout << "All done!" << endl;

  return (matches.isize() > 0);
}

template<class DataType>
void KmerAlignCore<DataType>::MergeSortFilter(svec<DataType> & result,
				    const svec<DataType> & one,
				    const svec<DataType> & two) {
  // Stupid for now...
  int i;
  
  result.clear();

  //cout << "one=" << one.isize() << "  two=" << two.isize() << endl;
  if (one.isize() == 0 || two.isize() == 0)
    return;

  svec<DataType> tmp;
  result.resize(one.isize() + two.isize());
  tmp.resize(one.isize() + two.isize());

  int k = 0;
  
  int x = 0;
  int y = 0;
  for (x=0; x<one.isize(); x++) {
    while (y<two.isize() && two[y] < one[x]) {
      y++;
    }
    if (x >= one.isize() || y >= two.isize())
      break;
    if (one[x] == two[y]) {
      result[k] = one[x];
      k++;
    }
  }
  result.resize(k);
}

template<class DataType>
bool KmerAlignCore<DataType>::IsRepeat(const NumVector & t, int i, int size, int min)
{
  if (t.size() == 0)
    return false;

  int r = 0; 
  int n = 0; 
  for (int j=i; j<i+size; j+=4) {
    if (t[j] >= min)
      r++;
    n++;
  }
  if (r >= n / 2) {
    return true;
  } else {
    return false;
  }
}

template<class DataType>
bool KmerAlignCore<DataType>::IsLowComplexity(const DNAVector & seq, int startIdx, int numOfBases)
{
  if (seq.size()-startIdx < numOfBases) { numOfBases = seq.size()-startIdx; }
  for (int i=startIdx; i<startIdx+numOfBases-1; i++) {
    if (seq[i] != seq[i+1]) { return false; }
  }
  return true;
}
 
//Explicit Template Instantiation
template class KmerAlignCore<KmerAlignCoreRecord>;
template class KmerAlignCore<KmerAlignCoreRecordWithScore>;



