#ifndef _ALIGNMENTCOORDINATES_H_
#define _ALIGNMENTCOORDINATES_H_

#include <sstream>
#include "../base/FileParser.h"

//==================================================================
/**
 * The AlignmentBlock class represents an allignment through 
 * its coordinates from the query and target sequences.
 * A parsed line  containing such coordinates in the blast/satsuma.... 
 * format can be used to construct an object. There are also
 * utility functions such as the merging of two blocks, etc. 
 **/
class AlignmentBlock 
{
  public:
    AlignmentBlock():targetStart(-1), targetStop(-1), queryStart(-1),
        queryStop(-1), orientation(' '), targetChrom(""),
        queryChrom("") {}

    /** Constructor that accepts a parser object and inserts fields */
    AlignmentBlock(FlatFileParser& parser) { parse(parser, false); }

    //TODO fix comment - adobted from SatsumaBlcok
    bool operator < (const AlignmentBlock& b) const {
      if (getTargetChrom() != b.getTargetChrom())
        return (getTargetChrom() < b.getTargetChrom());
      return (getTargetStart() < b.getTargetStart());
    }
    void set(const string& tChr, int tStart, int tStop) {
      targetChrom = tChr;
      targetStart = tStart;
      targetStop = tStop;
    }

    string toString() const; 
    void print () const; 

    /** Use parser object to read next line and parse into the object fields */
    bool parse(FlatFileParser& parser, bool flip);

    /** Merge two alignment blocks if they are within a given boundary and can be merged*/
    bool merge(const AlignmentBlock& toMerge, int mergeBoundary ); 

    /** 
     * Make sure the alignment blocks are from the same chromosomes and orientation
     * Used to make checks before merge and also in other cases for compatability
     */
    bool isCompatible(const AlignmentBlock& otherBlock)const; 

    bool isReversed() const { return (orientation=='-'); }

    int getTargetStart() const            { return targetStart; }
    int getTargetStop() const             { return targetStop;  }
    int getQueryStart() const             { return queryStart;  }
    int getQueryStop() const              { return queryStop;   } 
    char getOrient() const                { return orientation; }
    const string & getTargetChrom() const { return targetChrom; }
    const string & getQueryChrom() const  { return queryChrom;  }

    void setTargetStart(int tstr)         { targetStart = tstr; }
    void setTargetStop(int tstp)          { targetStop  = tstp; }
    void setQueryStart(int qstr)          { queryStart  = qstr; }
    void setQueryStop(int qstp)           { queryStop   = qstp; } 
    void setOrient(char ori)              { orientation = ori;  }
    void setTargetChrom(const string& tc) { targetChrom = tc;   }
    void setQueryChrom(const string& qc)  { queryChrom  = qc;   }


  private:
    /** Helper function to check whether num1 is within a +/- mergeBoundary of num2 */ 
    bool isWithinBound(int num1, int num2, int mergeBoundary);


    int    targetStart;  /// Basepair index on the target chromoseome where alignment starts
    int    targetStop;   /// Basepair index on the target chromoseome where alignment ends 
    int    queryStart;   /// Basepair index on the query chromoseome where alignment starts
    int    queryStop;    /// Basepair index on the query chromoseome where alignment ends
    char   orientation;  /// Orientation of the alignment (can be + or - for reverse)
    string targetChrom;  /// Identifying name for the target chromosome
    string queryChrom;   /// Identifying name for the query chromosome
};


#endif //_ALIGNMENTCOORDINATES_H_
