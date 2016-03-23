#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include <fstream>
#include <iomanip>
#include <cmath>
#include "Alignment.h"

//=====================================================================
void Alignment::getSeqAuxInfo(int& tOrigOffset, int& qOrigOffset, char& tSeqStrand, char& qSeqStrand) {
  tOrigOffset      = targetOrigOffset;              
  qOrigOffset      = queryOrigOffset;          
  tSeqStrand       = targetSeqStrand;        
  qSeqStrand       = querySeqStrand;              
} 

void Alignment::setSeqAuxInfo(int tOrigOffset, int qOrigOffset, char tSeqStrand, char qSeqStrand) {
  targetOrigOffset = tOrigOffset;              
  queryOrigOffset  = qOrigOffset;          
  targetSeqStrand  = tSeqStrand;        
  querySeqStrand   = qSeqStrand;              
} 

void Alignment::setSeqAuxInfo(int tOrigOffset, int qOrigOffset, bool tSeqStrand, bool qSeqStrand) {
  targetOrigOffset = tOrigOffset;              
  queryOrigOffset  = qOrigOffset;          
  targetSeqStrand  = (tSeqStrand?'+':'-');   
  querySeqStrand   = (qSeqStrand?'+':'-');              
}

void Alignment::print(int outputType, double pValLimit, ostream& sout, int screenWidth, bool withInfo) const{
  switch(outputType) {
    case 0: printFull(pValLimit, sout, screenWidth, withInfo); break;
    case 1: printInfoCSV(sout); break;
    case 2: printMFAFormat(pValLimit, sout, screenWidth); 
  }
}

int CountValid(const string & s)
{
  int k = 0;
  for (unsigned int i=0; i<s.length(); i++) {
    if (s[i] != '-')
      k++;
  }
  return k;
} 

string Alignment::toString(int screenWidth, bool withInfo) const{
  stringstream sout;
  if(withInfo) {
    sout << "**********************************************"              << endl
         << "Target sequence size:             " << getTargetLength()     << endl
         << "Query sequence size:              " << getQueryLength()      << endl
         << "Target offset:                    " << getTargetOffset()     << endl
         << "Query offset:                     " << getQueryOffset()      << endl 
         << "Target aligned basepairs:         " << getTargetBaseAligned()<< endl
         << "Query aligned basepairs:          " << getQueryBaseAligned() << endl
         << "Raw Score:                        " << getRawScore()         << endl 
         << "Identity score:                   " << calcIdentityScore()   << endl
         << "Total Edit Count                  " << calcEditCount()       << endl
         << "Mean Contiguity length            " << calcMeanContig()      << endl
         << "Mod-Smith-waterman score:         " << calcModSWScore()      << endl
         << "Significance P-value:             " << calcPVal()            << endl
         << "***********************************************"             << endl;
  } 
  // Let's not count gaps.
  int countQ = getQueryOffset();
  int countT = getTargetOffset();
  for( int i=0; i<=(int)(matchesStr.size()/screenWidth); i++) {
    string query = queryStr.substr(i*screenWidth,screenWidth); 
    sout << "Query: " << setw(6) << countQ  << " "
         << query << " " 
         << countQ + CountValid(query)-1 << endl
         << setw(14) << " " << matchesStr.substr(i*screenWidth,screenWidth) << endl
         << "Sbjct: " << setw(6) << countT << " " 
         << targetStr.substr(i*screenWidth,screenWidth) << " " 
         << countT + CountValid(targetStr.substr(i*screenWidth,screenWidth)) - 1 << endl
         << endl << endl << endl;
    countQ += CountValid(query);
    countT += CountValid(targetStr.substr(i*screenWidth,screenWidth));
  }
  return sout.str();
}

void Alignment::printFull(double pValLimit, ostream& sout, int screenWidth, bool withInfo) const{
  if(calcPVal()>pValLimit) { return; } // Significance is less than min required
  sout << toString(screenWidth, withInfo);
}

void Alignment::printMFAFormat(double pValLimit, ostream& sout, int screenWidth) const {

  int countQ = getQueryOffset()  + queryOrigOffset;
  int countT = getTargetOffset() + targetOrigOffset;
  if (targetSeqStrand == '-') {
    countT = targetOrigOffset + 
      getTargetLength() - 
      getTargetBaseAligned() - 
      getTargetOffset();
  }


  //sout << "DEBUG " << getTargetOffset() << " " << targetOrigOffset << " " << countT << " " << queryOrigOffset << endl;

  for( int i=0; i<=(int)(matchesStr.size()/screenWidth); i++) {
    string query = queryStr.substr(i*screenWidth,screenWidth); 
    int setwSize = max(getQuerySeq().Name().size(), getTargetSeq().Name().size()); 
    sout << std::left << setw(setwSize) << getQuerySeq().Name() << " " <<  querySeqStrand << " "  << std::right << setw(10) << countQ  << " "
         << query << " " 
         << countQ + CountValid(query)-1 << endl
         << std::left << setw(setwSize+14) << " " << matchesStr.substr(i*screenWidth,screenWidth) << endl
         << std::left << setw(setwSize) << getTargetSeq().Name() << " " <<  targetSeqStrand  << " " << std::right  << setw(10) << countT << " " 
         << targetStr.substr(i*screenWidth,screenWidth) << " " 
         << countT + CountValid(targetStr.substr(i*screenWidth,screenWidth)) - 1 << endl
         << endl << endl << endl;
    countQ += CountValid(query);
    countT += CountValid(targetStr.substr(i*screenWidth,screenWidth));
  }
} 

void Alignment::printXMLFormat(double pValLimit, ostream& sout, int screenWidth) const {

  sout <<"<Hsp>" << endl
       << "\t" << "<Hsp_num>"         << ""                   << "</Hsp_num>"         << endl
       << "\t" << "<Hsp_bit-score>"   << getSWScore()         << "</Hsp_bit-score>"   << endl
       << "\t" << "<Hsp_score>"       << getRawScore()        << "</Hsp_score>"       << endl
       << "\t" << "<Hsp_evalue>"      << getEValue()          << "</Hsp_evalue>"      << endl
       << "\t" << "<Hsp_query-from>"  << getQueryOffset()     << "</Hsp_query-from>" << endl
       << "\t" << "<Hsp_query-to>"    << getQueryOffset() + getQueryBaseAligned()  << "</Hsp_query-to>" << endl
       << "\t" << "<Hsp_hit-from>"    << getTargetOffset()    << "</Hsp_hit-from>"    << endl
       << "\t" << "<Hsp_hit-to>"      << getTargetOffset() + getTargetBaseAligned()  << "</Hsp_hit-to>" << endl
       << "\t" << "<Hsp_query-frame>" << ""                   << "</Hsp_query-frame>" << endl
       << "\t" << "<Hsp_hit-frame>"   << ""                   << "</Hsp_hit-frame>"   << endl
       << "\t" << "<Hsp_identity>"    << getIdentityScore()   << "</Hsp_identity>"    << endl
       << "\t" << "<Hsp_positive>"    << ""                   <<"</Hsp_positive>"     << endl
       << "\t" << "<Hsp_align-len>"   << getAlignmentLen()    <<"</Hsp_align-len>"    << endl
       << "\t" << "<Hsp_qseq>"        << getQueryString()     <<"</Hsp_qseq>"         << endl
       << "\t" << "<Hsp_hseq>"        << getTargetString()    <<"</Hsp_hseq>"         << endl
       << "\t" << "<Hsp_midline>"     << getMatchString()     <<"</Hsp_midline>"      << endl
       << "</Hsp>" << endl;

}

 
void Alignment::printInfoCSV( ostream& sout) const {
  sout<< getTargetOffset()     <<","
      << getQueryOffset()      <<","
      << getTargetBaseAligned()<<","
      << getQueryBaseAligned() <<","
      << calcEditCount()       <<","
      << calcMeanContig()      <<","
      << getRawScore()         <<","
      << calcIdentityScore()   <<","
      << getSWScore()          <<","
      << calcPVal()            <<","
      << getRuntime()          <<","
      << getRuntimeCoef()      <<endl;
}

char Alignment::getQueryAlignCharForTarget(int targetIndex) const { 
  int qIdx = targetIdxsInQuery[targetIndex];
  if( qIdx >= 0 ) { 
    return querySeq[qIdx]; 
  } else {
    return GAP_CHAR; 
  }
}

int Alignment::getQueryAlignIndexForTarget(int targetIndex) const { 
  return targetIdxsInQuery[targetIndex];
}

char Alignment::getTargetAlignCharForQuery(int queryIndex) const { 
  int tIdx = queryIdxsInTarget[queryIndex];
  if( tIdx >= 0 ) { 
    return targetSeq[tIdx]; 
  } else {
    return GAP_CHAR; 
  }
}

int Alignment::getTargetAlignIndexForQuery(int targetIndex) const { 
  return targetIdxsInQuery[targetIndex];
}

double Alignment::calcMeanContig() const {
  int lenCount = 0;
  int totCount = 0;
  int totSum   = 0;
  int matchesLen = matchesStr.size();
  for(int i=0; i<matchesLen; i++) {
    if(matchesStr[i]==MATCH_CHAR) { lenCount++; }
    else { 
      totSum += lenCount; 
      if(lenCount != 0 ) { totCount++; }   
      lenCount  = 0; //reset contig count  
    }
  }  
  totSum += lenCount; 
  if(lenCount !=0 || totCount==0 ) { totCount++; }   
  return double(totSum)/totCount;
}

int Alignment::calcEditCount() const {
  int numOfEdits   = 0;  
  int gapCnt       = 0;
  // Note that one indel (regardless of its length) is considered as one edit
  int matchesLen = matchesStr.size();
  for(int i=0; i<matchesLen; i++) {
    if(targetStr[i]==GAP_CHAR || queryStr[i]==GAP_CHAR) { gapCnt++;}
    else {
      if(gapCnt != 0) { 
        numOfEdits++; // gap ended
        gapCnt = 0;   //reset 
      }
      if(targetStr[i] != queryStr[i]) { numOfEdits++; }  //Mismatch
    }
  }
  return numOfEdits;
}

int Alignment::calcModSWScore() const {
  int totScore = 0;  
  int gapCnt   = 0;
  int matchesLen = matchesStr.size();
  for(int i=0; i<matchesLen; i++) {
    if(targetStr[i]==GAP_CHAR || queryStr[i]==GAP_CHAR) { gapCnt++;}
    else {
      if(gapCnt != 0) { 
        int maxPenalty = min(gapCnt, 10); // Limit gap penalty to 10
        totScore -= maxPenalty; 
        gapCnt = 0;   //reset - gap ended
      }
      if(targetStr[i] != queryStr[i]) { totScore--; }  //Mismatch
      else { totScore++; } //Match
    }
  }
  return totScore;
}

double Alignment::calcPVal() const {
  double lambda = 0.51;
  double mu     = 15.0;
  double x      = calcModSWScore();
  double p      = 1 - exp(-exp(-lambda*(x-mu)));
  return p;
}

void Alignment::resetContainers() {
  targetStr  = "";
  queryStr   = "";
  matchesStr = "";
  targetIdxsInQuery.clear();
  queryIdxsInTarget.clear();
}

//=====================================================================
void AlignmentInfo::resetCalculations(int tos, int qos) {
  tBaseAligned       = 0;
  qBaseAligned       = 0;
  baseMatched        = 0;
  alignmentLen       = 0;
  smithWatermanScore = 0;

  // Set offset values
  tOffset            = tos;
  qOffset            = qos;  
}

double AlignmentInfo::calcIdentity() const {
  return double(baseMatched)/getMaxBaseAligned();
}

