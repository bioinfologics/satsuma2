#ifndef _COORDINATE_H_
#define _COORDINATE_H_

#include <map>
#include <string>
#include <sstream>
#include "base/SVector.h"
#include "src/AlignmentBlock.h"
//#include "base/FileParser.h"
//#include "src/AnnotationQuery/NCList.h"

//======================================================
/** Generic DNA Sequence Coordinates */
class Coordinate {
public:
  Coordinate():chr(""), orient(false), start(-1), stop(-1) {}

  Coordinate(string ch, bool ori, int str, int stp)
          :chr(ch), orient(ori), start(str), stop(stp) {}

  bool operator<(const Coordinate & c) const {
    if (getChr() != c.getChr())      { return (getChr() < c.getChr()); }
    return (getStart() < c.getStart()); 
  }
  bool operator==(const Coordinate & c) const {
    return ( isSameCoords(c) && isSameOrient(c) );
  }

  const string & getChr() const   { return chr;                     }
  char getOrient() const          { return orient?'+':'-';          }
  int getStart() const            { return start;                   }
  int getStop() const             { return stop;                    }
  void setChr(const string& ch)   { chr    = ch;                    }
  void setOrient(bool ori)        { orient = ori;                   }
  void setOrient(char ori)        { orient = (ori=='-')?false:true; }
  void setStart(int str)          { start  = str;                   }
  void setStop(int stp)           { stop   = stp;                   }
  // Backward compatiblity with GenomeCoord class
  void SetChr(char ch, char ori) {
    chr = ch;
    orient = ori;
  }
  bool isReversed() const { return !orient; }
  void set(const string& ch, bool ori, int str, int stp);
  void set(const string& ch, char ori, int str, int stp); 

  /** Return the length of the coordinate by subtracting start from stop (absolute value) */
  int findLength() const { return (abs(getStop()-getStart())+1); }
  /** Returns true if both coords have the same orientation  */
  bool isSameOrient(const Coordinate& other) const {return (other.getOrient() == getOrient()); }
  /** Returns true if the two are on the same chromosome & share at least one base */
  bool hasOverlap(const Coordinate& other) const;  
  /** Returns true if this coord contains other */
  bool contains(const Coordinate& other) const;    
  /** Returns true if start and stop and chromosome of both coords are the same */
  bool isSameCoords(const Coordinate& other) const; 
  /** Returns the percentage of basepairs overlapping over total length of item */
  int findOverlapCnt(const Coordinate& other) const;  

  /** 
   * Creates a string containing the details of the object
   * A separator character is provided to use for separating the fields
   */
  string toString(char sep) const;
  string toString_noOrient(char sep) const; 
  void print()const { cout << toString('\t') << endl; }
  void setFromTarget(const AlignmentBlock & b);

private:
  string chr;    /// Chromosome
  bool orient;   /// Orientation (true for forward (+) and false otherwise (-))
  int  start;    /// Start point
  int  stop;     /// Stop point 
};

#endif //_COORDINATE_H_
