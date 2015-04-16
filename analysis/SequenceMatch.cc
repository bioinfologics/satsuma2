#include "analysis/SequenceMatch.h"
#include "util/mutil.h"

#include "base/FileParser.h"



void SingleMatch::Read(CMReadFileStream & f, int version)
{
  f.Read(m_targetID);
  f.Read(m_queryID);
  f.Read(m_queryLen);
  f.Read(m_startTarget);
  f.Read(m_startQuery);
  f.Read(m_length);
  int n;
  f.Read(n);
  if (n == 1)
    m_bRC = true;
  else
    m_bRC = false;
  f.Read(m_matches);
  f.Read(m_prob);
  if (version >= 3)
    f.Read(m_ident);
}

void SingleMatch::ReadAppend(CMReadFileStream & f, int version, int app)
{
  f.Read(m_targetID);
  f.Read(m_queryID);
  m_queryID += app;
  f.Read(m_queryLen);
  f.Read(m_startTarget);
  f.Read(m_startQuery);
  f.Read(m_length);
  int n;
  f.Read(n);
  if (n == 1)
    m_bRC = true;
  else
    m_bRC = false;
  f.Read(m_matches);
  f.Read(m_prob);
  if (version >= 3)
    f.Read(m_ident);
}

void SingleMatch::Write(CMWriteFileStream & f)
{
  f.Write(m_targetID);
  f.Write(m_queryID);
  f.Write(m_queryLen);
  f.Write(m_startTarget);
  f.Write(m_startQuery);
  f.Write(m_length);
  int n = 0;
  if (m_bRC)
    n = 1;
  f.Write(n);
  f.Write(m_matches);
  f.Write(m_prob);
  f.Write(m_ident);
}


//================================================================
void MultiMatches::Merge(const MultiMatches & in)
{
  int i;

  cout << "Start merging." << endl;
  if (m_targetNames.isize() == 0) {
    SetCounts(in.GetTargetCount(), in.GetQueryCount());
    for (i=0; i<in.GetTargetCount(); i++) {
      SetTargetSize(i, in.GetTargetSize(i));
      SetTargetName(i, in.GetTargetName(i));
    }
  }
  for (i=0; i<in.GetQueryCount(); i++) {
    SetQuerySize(i, in.GetQuerySize(i));
    SetQueryName(i, in.GetQueryName(i));
  }
  
  cout << "Before: " << m_count << endl;
  m_matches.resize(m_count + in.GetMatchCount());
  cout << "After: " << m_count + in.GetMatchCount() << endl;
  
  for (i=0; i<in.GetMatchCount(); i++) {
    m_matches[i+m_count] = in.GetMatch(i);
  }

  m_count = m_matches.isize();

}

int MultiMatches::Update(svec<string> & names, svec<int> & size, const string & n, int s)
{
  int i;
  for (i=0; i<names.isize(); i++) {
    if (names[i] == n) {
      if (s+1 > size[i])
	size[i] = s+1;
      return i;
    }
  }
  names.push_back(n);
  size.push_back(s);

  return names.isize()-1;
}

void MultiMatches::Read(const string & file)
{  
  m_targetNames.clear();
  m_queryNames.clear();
  m_matches.clear();
  m_count = 0;

  MergeRead(file);
}

//====================================================
void MultiMatches::ReadSummary(const string & file)
{
  int i;
  
  FlatFileParser parser;
  
  parser.Open(file);
  
  string lastT, lastQ;
  int k = 0;
  int tID = -1;
  int qID = -1;
  while (parser.ParseLine()) {
    if (k >= m_matches.isize())
      m_matches.resize(k + 500000);
    
    const string & t = parser.AsString(0);
    const string & q = parser.AsString(3);
    if (lastT != t)
      tID = Update(m_targetNames, m_targetSize, t, parser.AsInt(2));
    if (lastQ != q) 
      qID = Update(m_queryNames, m_querySize, q, parser.AsInt(5));
    SingleMatch & m = m_matches[k];
    
    m.SetQueryTargetID(qID, tID, 0);
    bool bRC = false;
    //if (parser.AsString(7) == "-")
    //bRC = true;

    m.SetPos(parser.AsInt(4), 
	     parser.AsInt(1), 
	     parser.AsInt(2)-parser.AsInt(1), 
	     bRC);
    m.SetIdentity(parser.AsFloat(6));    

 
    k++;
    m.SetProbability(1.);
    if (k % 10000 == 0) {
      cout << "Processed: " << k << endl;
    }

    lastT = t;
    lastQ = q;
    
  }
  cout << "Total matches: " << k << endl;
  cout << "Target seq's:  " << m_targetNames.isize() << endl;
  cout << "Query seq's:   " << m_queryNames.isize() << endl;
  m_matches.resize(k);
  m_count = m_matches.isize();
}

void MultiMatches::MergeRead(const string & file)
{
  CMReadFileStream f;
  
  f.Open(file.c_str());
  int ver = -1;

  f.Read(ver);

  cout << "File version " << ver << endl;
  if (ver == -1) {
    cout << "I/O ERROR: file is corrupt!!" << endl;
    return;
  }

  if (ver > m_currVer) {
    cout << "Incorrect binary version, trying ASCII..." << endl;
    f.Close();
    ReadSummary(file);
    return;
  }

  int i;
  int n;

  f.Read(n);
  m_targetNames.resize(n);
  for (i=0; i<m_targetNames.isize(); i++) {
    CMString s;
    f.Read(s);
    m_targetNames[i] = (const char*)s;
  }
  f.Read(n);
  m_queryNames.resize(n);
  for (i=0; i<m_queryNames.isize(); i++) {
    CMString s;
    f.Read(s);
    m_queryNames[i] = (const char*)s;
  }

  int start = m_count;

  int len = 0;

  f.Read(len);
  //cout << "Count= " << m_count << endl;

  //int start = m_matches.isize();
  m_matches.resize(len + m_count);
  m_count = m_matches.isize();
  
  for (i=start; i<m_count; i++)
    m_matches[i].Read(f, ver);


  m_targetSize.clear();
  m_querySize.clear();

  m_targetSize.resize(m_targetNames.isize(), 0);
  m_querySize.resize(m_queryNames.isize(), 0);

  if (ver > 1) {
    for (i=0; i<m_targetSize.isize(); i++)
      f.Read(m_targetSize[i]);
    for (i=0; i<m_querySize.isize(); i++)
      f.Read(m_querySize[i]);
  }

  f.Close();
}

void MultiMatches::MergeReadAppend(const string & file)
{
  CMReadFileStream f;
  
  f.Open(file.c_str());
  int ver = 1;

  f.Read(ver);

  cout << "File version " << ver << endl;

  int i;
  int n;

  //m_matches.resize(m_count);
  
  f.Read(n);
  
  cout << "Target size=" << m_targetNames.isize() << ", appending " << n << endl;
  int before = m_targetNames.isize();
  m_targetNames.resize(n+before);
  for (i=before; i<m_targetNames.isize(); i++) {
    CMString s;
    f.Read(s);
    m_targetNames[i] = (const char*)s;
  }

  f.Read(n);
  int curr = m_queryNames.isize();
  m_queryNames.resize(n);
  cout << "Query size=" << curr << endl;

  
  for (i=0; i<m_queryNames.isize(); i++) {
    CMString s;
    f.Read(s);
    m_queryNames[i] = (const char*)s;
  }

  f.Read(m_count);
  //cout << "Count= " << m_count << endl;

  int start = m_matches.isize();
  m_matches.resize(m_matches.isize() + m_count);
  m_count = m_matches.isize();
  
  for (i=start; i<m_count; i++)
    m_matches[i].ReadAppend(f, ver, curr);


  m_targetSize.clear();
  

  //m_querySize.clear();

  m_targetSize.resize(m_targetNames.isize(), 0);
  m_querySize.resize(m_queryNames.isize(), 0);

  if (ver > 1) {
    for (i=before; i<m_targetSize.isize(); i++)
      f.Read(m_targetSize[i]);
    for (i=curr; i<m_querySize.isize(); i++)
      f.Read(m_querySize[i]);
  }

  f.Close();
}


void MultiMatches::Write(const string file)
{
  CMWriteFileStream f;
  
  f.Open(file.c_str());
  int ver = m_currVer;

  f.Write(ver);

  int i;
  
  f.Write(m_targetNames.isize());
  for (i=0; i<m_targetNames.isize(); i++) {
    CMString s = m_targetNames[i].c_str();
    f.Write(s);
  }
  f.Write(m_queryNames.isize());
  for (i=0; i<m_queryNames.isize(); i++) {
    CMString s = m_queryNames[i].c_str();
    f.Write(s);
  }

  //cout << "Count= " << m_count << endl;
  f.Write(m_count);

  for (i=0; i<m_count; i++)
    m_matches[i].Write(f);

  for (i=0; i<m_targetSize.isize(); i++)
    f.Write(m_targetSize[i]);
  for (i=0; i<m_querySize.isize(); i++)
    f.Write(m_querySize[i]);

  

  f.Close();

}

void MultiMatches::SetTargetSize(int i, int size)
{
  if (m_targetSize.isize() == 0)
    m_targetSize.resize(m_targetNames.isize(), 0);
  m_targetSize[i] = size;
}

void MultiMatches::SetQuerySize(int i, int size)
{
  if (m_querySize.isize() == 0)
    m_querySize.resize(m_queryNames.isize(), 0);
  m_querySize[i] = size;
}

void MultiMatches::SetCounts(int target, int query)
{
  m_targetNames.resize(target);
  m_queryNames.resize(query);
  m_targetSize.resize(target, 0);
  m_querySize.resize(query, 0);
}

void MultiMatches::SetTargetName(int i, const string & n)
{
  m_targetNames[i] = n;
}

void MultiMatches::SetQueryName(int i, const string & n)
{
  m_queryNames[i] = n;
}


void MultiMatches::SetNames(const svec<string> & qNames, const svec<string> & tNames) {
  m_targetNames.resize(tNames.size());
  m_queryNames.resize(qNames.size());
  int i;
  for (i=0; i<m_targetNames.isize(); i++) {
    m_targetNames[i] = tNames[i];
  }
  for (i=0; i<m_queryNames.isize(); i++) {
    m_queryNames[i] = qNames[i];
  }
}

bool Laps(int a, int b)
{
  int c = a-b;
  if (c < 0)
    c = -c;
  if (c < 4)
    return true;
  else
    return false;
}

void MultiMatches::Collapse()
{
  int i;
  svec<SingleMatch> tmp;
  tmp.resize(m_matches.isize());
  int k = 0;

  SingleMatch n = m_matches[0];
  cout << "Matches before collapse: " << m_count << endl;
  for (i=1; i<m_count; i++) {
    const SingleMatch & s1 = m_matches[i-1];
    const SingleMatch & s2 = m_matches[i];

    //cout << "start=" << s2.GetStartTarget() << "len=" << s2.GetLength() << endl;

    if (s1.GetTargetID() == s2.GetTargetID() &&
	s1.IsRC() == s2.IsRC() &&
	Laps(s1.GetStartTarget(), s2.GetStartTarget()) &&
	Laps(s1.GetStartQuery(), s2.GetStartQuery())) {
      if (s2.GetStartQuery() < n.GetStartQuery()) {
	//cout << "BEFORE: start=" << n.GetStartQuery() << " len=" << n.GetLength() << endl;
	//cout << "THIS:   start=" << s2.GetStartQuery() << " len=" << s2.GetLength() << endl;
	//cout << "AND:    start=" << s1.GetStartQuery() << " len=" << s1.GetLength() << endl;
	n.SetPos(s2.GetStartQuery(), n.GetStartTarget(), 
		 n.GetStartQuery() + n.GetLength() - s2.GetStartQuery(),
		 n.IsRC());
 	//cout << "AFTER: start=" << n.GetStartQuery() << " len=" << n.GetLength() << endl;
     } else {
	//cout << "BEFORE: start=" << n.GetStartQuery() << " len=" << n.GetLength() << endl;
	//cout << "THIS:   start=" << s2.GetStartQuery() << " len=" << s2.GetLength() << endl;
	n.SetPos(n.GetStartQuery(), n.GetStartTarget(), 
		 s2.GetStartQuery() + s2.GetLength() - n.GetStartQuery(),
       		 n.IsRC());
	//cout << "AFTER: start=" << n.GetStartQuery() << " len=" << n.GetLength() << endl;
     }

    } else {	
      //if (s1.IsRC())
      //cout << "New one." << endl;
      tmp[k] = n;
      k++;
      n = s2;
    }
  }

  m_count = k;
  tmp.resize(k);
  m_matches = tmp;
  cout << "Matches after collapse:  " << m_count << endl;


}
