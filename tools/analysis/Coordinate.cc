#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "analysis/Coordinate.h"



//======================================================
void Coordinate::set(const string& ch, bool ori,
         int str, int stp) {
  setChr(ch);
  setOrient(ori);
  setStart(str);
  setStop(stp);
}

void Coordinate::set(const string& ch, char ori,
         int str, int stp) {
  setChr(ch);
  setOrient(ori);
  setStart(str);
  setStop(stp);
}

bool Coordinate::contains(const Coordinate& other) const {
  return( (other.getChr() == getChr())
          && ( (other.getStart()>=getStart())&& (other.getStop()<=getStop())
          )
        );
}

bool Coordinate::isSameCoords(const Coordinate& other) const {
  return( (other.getChr() == getChr())
          && ( (other.getStart()==getStart())&& (other.getStop()==getStop())
          )
        );
}

bool Coordinate::hasOverlap(const Coordinate& other) const {
 //For any overlap either the start or stop of one of the coords has to lie within the other coord
  return( (other.getChr() == getChr())
          && !( (other.getStart()>getStop())|| (other.getStop()<getStart())
          )
        );
}

int Coordinate::findOverlapCnt(const Coordinate& other) const {
  if(!hasOverlap(other) || !isSameOrient(other)) { return .0; }
  int overlap   = min(getStop(), other.getStop()) - max(getStart(), other.getStart()) + 1;
  return overlap; 
}

string Coordinate::toString(char sep) const {
  stringstream outStream;
  outStream << getChr() << sep << getStart() << sep << getStop() << sep << getOrient(); 
  return outStream.str();
}

string Coordinate::toString_noOrient(char sep) const {
  stringstream outStream;
  outStream << getChr() << sep << getStart() << sep << getStop(); 
  return outStream.str();
}

void Coordinate::setFromTarget(const AlignmentBlock& b) {
  set(b.getTargetChrom(), b.getOrient(), b.getTargetStart(), b.getTargetStop());
}

//======================================================
